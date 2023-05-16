# *****************************************************************************
# R Script implementing conventional random f-fold cross-validation.  
# Related to the paper "Dealing with clustered samples for assessing map 
# accuracy by cross-validation".
# Contact: Sytze de Bruin, Wageningen University, Laboratory of Geo-information
# Science and Remote Sensing, email: sytze.debruin@wur.nl
# May 3, 2022
# *****************************************************************************

# ****** load required library *******
#Sys.sleep(round(runif(1, min = 1, max = 240)))

.libPaths("/home/j/jlinnenb/r_packages/")
library(ranger)
library(CAST)
library(sf)
library(raster)
library(caret)
library(parallel)
library(terra)

mecfu <- function(ref, pred){
  muref <- global(ref, "mean", na.rm=T)[[1]]
  residsq <- (ref - pred)^2
  SSR <- global(residsq, "sum", na.rm=T)[[1]]
  rm(residsq)
  residsq <- (ref - muref)^2
  SST <- global(residsq, "sum", na.rm=T)[[1]]
  1 - SSR/SST
}

rmsefu <- function(ref, pred){
  residsq <- (ref - pred)^2
  sqrt(global(residsq, "mean", na.rm=T)[[1]])
}

mefu <- function(ref, pred){
  resmap <- ref - pred
  global(resmap, "mean", na.rm=T)[[1]]
}

AOA_perc <- function(newdata, model, indexTrain=NULL, indexTest=NULL) {
  AOA <- aoa(newdata, model=model, CVtrain = indexTrain, CVtest = indexTest)$AOA
  AOAn <- terra::freq(AOA)[,2:3]
  overall_cells <- sum(AOAn$count)
  AOAp <- AOAn$count / overall_cells * 100
  list(data.frame(inside_AOA = AOAp[[2]], outside_AOA = AOAp[[1]]), AOA)
}

# ************ GLOBALS ***************
setwd("/scratch/tmp/jlinnenb/deBruin_add_nndm/")

samples   <- c("clusterMedium", "clusterStrong", "clusterGapped", "regular", 
               "simpleRandom")
infolder1 <- "./data"
infolder2 <- "./samples"
outfolder <- "./CVresults"
startseed <- 1234567

# create outfolders if they don't exist
if(!dir.exists(outfolder))
  dir.create(outfolder)

if(!dir.exists(paste0(outfolder, "/ffs_knndm")))
  dir.create(paste0(outfolder, "/ffs_knndm"))

if(!dir.exists(paste0(outfolder, "/tmp_rast")))
  dir.create(paste0(outfolder, "/tmp_rast"))

if(!file.exists(file.path(outfolder, "ffs_knndm", "runs.csv"))) {
  write.csv(data.frame("runs"=1), file.path(outfolder, "ffs_knndm", "runs.csv"), row.names = FALSE)
}


csv_file <- file.path(outfolder, "ffs_knndm", "runs.csv")
runs <- read.csv(csv_file)
lastIndex <- runs[nrow(runs),1]
thisIndex <- lastIndex + 1
print(paste0("this Index is: ", thisIndex))
runs[thisIndex,1] <- thisIndex
write.csv(runs, file = csv_file, row.names = FALSE)

# download data from https://doi.org/10.5281/zenodo.6513429
# ****** load input raster data ******
msk <- rast(file.path(infolder1, "TOTmask.tif"))
AGBstack <- rast(file.path(infolder1, "AGBstack.tif"))
OCSstack <- rast(file.path(infolder1, "OCSstack.tif"))

# test
#msk <- msk[2000:2500, 2000:2500, drop=FALSE] 
#AGBstack <- terra::crop(AGBstack, msk)
#OCSstack <- terra::crop(OCSstack, msk)
#

AGB_glc <- terra::segregate(AGBstack$glc2017)
names(AGB_glc) <- paste0("glc2017", c(1:3, 5:8))

AGBstack$glc2017 <- NULL
AGBstack <- c(AGBstack, AGB_glc)
OCSstack$glc2017 <- NULL
OCSstack <- c(OCSstack, AGB_glc)

OCt <- rast(file.path(infolder1, "ocs.tif")) #|> crop(msk)
AGt <- rast(file.path(infolder1, "agb.tif")) #|> crop(msk)

OCS <- mask(OCt, msk, filename="tmpocs.tif", overwrite=T)
AGB <- mask(AGt, msk, filename="tmpagb.tif", overwrite=T)
rm(AGt, OCt)

# ************ FUNCTIONS ***************


knndmCV <- function(smpl, number, variate, seed){

  #smpl="clusterStrong"
  #number=5
  #variate="AGB"
  #seed=1
  
  message(paste0("number ", seed))
  
  fname  <-  paste0(variate, "_", smpl, sprintf("%03d", number), ".csv")
  f_out  <- file.path(outfolder,"ffs_knndm", fname)
  
  if(!file.exists(f_out)) {
    fname <- paste0(variate, "data", sprintf("%03d", number), ".Rdata")
    f_in <- file.path(infolder2,smpl,fname)
    load(f_in)
    
    # load ppoints
    load(file.path(infolder2, "ppoints.Rdata"))
    
    if(variate == "AGB"){
      pts_df <- data.frame(x=AGBdata$xcoord * 1000, y=AGBdata$ycoord * 1000)
      pts_sf <- st_as_sf(pts_df, coords = c("x","y"), crs = st_crs(ppoints))
      
      pts_sf <- extract(AGB, vect(pts_sf), bind=TRUE) |> 
        st_as_sf() |> 
        na.omit() |> 
        st_as_sfc()
      
      pts_sf <- st_as_sf(st_set_crs(pts_sf, st_crs(ppoints)))
      
    } else {
      pts_df <- data.frame(x=OCSdata$xcoord * 1000, y=OCSdata$ycoord * 1000)
      pts_sf <- st_as_sf(pts_df, coords = c("x","y"), crs = st_crs(ppoints))
      
      pts_sf <- extract(OCS, vect(pts_sf), bind=TRUE) |> 
        st_as_sf() |> 
        na.omit() |> 
        st_as_sfc()
      
      pts_sf <- st_as_sf(st_set_crs(pts_sf, st_crs(ppoints)))
    }
    
    set.seed(seed)
    knndm <- knndm(pts_sf, ppoints = ppoints, k = 10, maxp = 0.8)
    trControl = trainControl(method = "cv", savePredictions = "final",
                             index=knndm$indx_train)
    
    f_out_p <- paste0(outfolder, "/tmp_rast/tmp_pred.tif")
    
    set.seed(seed)
    if(variate == "AGB"){
      
      mtry <- 4
      pgrid <- expand.grid(splitrule="variance",min.node.size=5,mtry=mtry)
      
      AGBdata <- extract(AGBstack, vect(pts_sf), ID = FALSE)
      
      var_l <- apply(AGBdata, 2, function(x) length(unique(x)))
      if(any(var_l==1)) {
        pred_zero_var <- names(which(var_l==1))
        AGBstack <- AGBstack[[!names(AGBstack) %in% pred_zero_var]]
        AGBdata <- AGBdata[,!names(AGBdata) %in% pred_zero_var]
        message(paste("Predictors", pred_zero_var[1], "and", pred_zero_var[2],
                      "had zero variance and were removed"))
      }
      
      
      f_mod <- ffs(AGBdata[,names(AGBdata) != "agb"], response = AGBdata[,names(AGBdata) == "agb"],
                   respect.unordered.factors=TRUE, method = "ranger",importance="permutation",
                   tuneGrid=pgrid, num.trees=500, trControl = trControl, globalval = TRUE,
                   verbose=FALSE)
      
      f_sel_var <- rownames(varImp(f_mod)$importance)
      if(any(c("xcoord", "ycoord") %in% f_sel_var)) {
        xy_as_predictor <- "1"
        xy_imp <- varImp(f_mod)$importance
        xy_imp <- sum(xy_imp$importance[rownames(xy_imp$importance) %in% c("xcoord", "ycoord"),], na.rm=TRUE)
        
      } else xy_as_predictor <- "0"; xy_imp <- 0
      
      n_sel_vars <- length(f_sel_var)
      
      message("predicting ffs model")
      map_ffs  <- terra::predict(AGBstack[[f_sel_var]], f_mod, 
                                 filename=f_out_p, overwrite=TRUE, na.rm=TRUE)
      
      true_ME_ffs   <- mefu(AGB, map_ffs)
      true_RMSE_ffs <- rmsefu(AGB, map_ffs)
      true_MEC_ffs  <- mecfu(AGB, map_ffs)
      
      CV_RMSE_ffs <- global_validation(f_mod)[[1]]
      
      diff_RMSE_ffs <- (CV_RMSE_ffs-true_RMSE_ffs)/true_RMSE_ffs*100
      
      message("calculating AOA for FFS model")
      AOA_AGB_ffs <- AOA_perc(AGBstack[[f_sel_var]], f_mod,
                              indexTrain = knndm$indx_train)
      
      message("setting values outside AOA to NA")
      map_ffs[AOA_AGB_ffs[[2]] == 0] <- NA
      true_ME_ffs_inAOA   <- mefu(AGB, map_ffs)
      true_RMSE_ffs_inAOA <- rmsefu(AGB, map_ffs)
      true_MEC_ffs_inAOA  <- mecfu(AGB, map_ffs)
      
      diff_RMSE_ffs_inAOA <- (CV_RMSE_ffs-true_RMSE_ffs_inAOA)/true_RMSE_ffs_inAOA*100
      
      AOA_stats_ffs <- AOA_AGB_ffs[[1]]
      
      message("training model without FFS")
      mod <- train(AGBdata[,names(AGBdata) != "agb"], y = AGBdata[,names(AGBdata) == "agb"],
                   respect.unordered.factors=TRUE, method = "ranger", importance="permutation",
                   tuneGrid=pgrid, num.trees=500, trControl = trControl)
      message("predicting normal model")
      map  <- terra::predict(AGBstack, mod, filename=f_out_p, overwrite=TRUE, na.rm=TRUE)
      true_ME  <- mefu(AGB, map)
      true_RMSE <- rmsefu(AGB, map)
      true_MEC  <- mecfu(AGB, map)
      
      CV_RMSE <- global_validation(mod)[[1]]
      diff_RMSE <- (CV_RMSE-true_RMSE)/true_RMSE * 100
      message("calculating AOA without FFS")
      AOA_AGB <- AOA_perc(AGBstack, mod,
                          indexTrain = knndm$indx_train)
      message("setting values outside AOA to NA (without FFS)")
      map[AOA_AGB[[2]] == 0] <- NA
      true_ME_inAOA   <- mefu(AGB, map)
      true_RMSE_inAOA <- rmsefu(AGB, map)
      true_MEC_inAOA  <- mecfu(AGB, map)
      
      diff_RMSE_inAOA <- (CV_RMSE-true_RMSE_inAOA)/true_RMSE_inAOA*100
      
      AOA_stats <- AOA_AGB[[1]]
      
    } else{
      
      message("starting OCS")
      mtry <- 4
      pgrid <- expand.grid(splitrule="variance",min.node.size=5,mtry=mtry)
      
      OCSdata <- extract(OCSstack, vect(pts_sf), ID = FALSE)
      
      var_l <- apply(OCSdata, 2, function(x) length(unique(x)))
      if(any(var_l==1)) {
        pred_zero_var <- names(which(var_l==1))
        OCSstack <- OCSstack[[!names(OCSstack) %in% pred_zero_var]]
        OCSdata <- OCSdata[,!names(OCSdata) %in% pred_zero_var]
        message(paste("Predictors", pred_zero_var[1], "and", pred_zero_var[2],
                      "had zero variance and were removed"))
      }
      
      
      f_mod <- ffs(OCSdata[,names(OCSdata) != "ocs"], response = OCSdata[,names(OCSdata) == "ocs"],
                   respect.unordered.factors=TRUE, method = "ranger",importance="permutation",
                   tuneGrid=pgrid, num.trees=500, trControl = trControl, globalval = TRUE)
      
      f_sel_var <- rownames(varImp(f_mod)$importance)
      if(any(c("xcoord", "ycoord") %in% f_sel_var)) {
        xy_as_predictor <- "1"
        xy_imp <- varImp(f_mod)$importance
        xy_imp <- sum(xy_imp$importance[rownames(xy_imp$importance) %in% c("xcoord", "ycoord"),], na.rm=TRUE)
        
      } else xy_as_predictor <- "0"; xy_imp <- 0
      
      map_ffs  <- terra::predict(OCSstack[[f_sel_var]], f_mod, filename=f_out_p, overwrite=TRUE, na.rm=TRUE)
      true_ME_ffs   <- mefu(OCS, map_ffs)
      true_RMSE_ffs <- rmsefu(OCS, map_ffs)
      true_MEC_ffs  <- mecfu(OCS, map_ffs)
      
      CV_RMSE_ffs <- global_validation(f_mod)[[1]]
      diff_RMSE_ffs <- (CV_RMSE_ffs-true_RMSE_ffs)/true_RMSE_ffs*100
      
      AOA_OCS_ffs <- AOA_perc(OCSstack[[f_sel_var]], f_mod)
      map_ffs[AOA_OCS_ffs[[2]] == 0] <- NA
      true_ME_ffs_inAOA   <- mefu(OCS, map_ffs)
      true_RMSE_ffs_inAOA <- rmsefu(OCS, map_ffs)
      true_MEC_ffs_inAOA  <- mecfu(OCS, map_ffs)
      
      diff_RMSE_ffs_inAOA <- (CV_RMSE_ffs-true_RMSE_ffs_inAOA)/true_RMSE_ffs_inAOA*100
      
      AOA_stats_ffs <- AOA_OCS_ffs[[1]]
      
      
      mod <- train(OCSdata[,names(OCSdata) != "ocs"], y = OCSdata[,names(OCSdata) == "ocs"],
                   respect.unordered.factors=TRUE, method = "ranger", importance="permutation",
                   tuneGrid=pgrid, num.trees=500, trControl = trControl)
      map  <- terra::predict(OCSstack, mod, filename=f_out_p, overwrite=TRUE, na.rm=TRUE)
      true_ME  <- mefu(OCS, map)
      true_RMSE <- rmsefu(OCS, map)
      true_MEC  <- mecfu(OCS, map)
      
      CV_RMSE <- global_validation(mod)[[1]]
      diff_RMSE <- (CV_RMSE-true_RMSE)/true_RMSE * 100
      
      AOA_OCS <- AOA_perc(OCSstack, mod)
      map[AOA_OCS[[2]] == 0] <- NA
      true_ME_inAOA   <- mefu(OCS, map)
      true_RMSE_inAOA <- rmsefu(OCS, map)
      true_MEC_inAOA  <- mecfu(OCS, map)
      
      diff_RMSE_inAOA <- (CV_RMSE-true_RMSE_inAOA)/true_RMSE_inAOA*100
      
      AOA_stats <- AOA_OCS[[1]]
      message(paste0("OCS stats:", AOA_stats))
      
    }
    
    write.csv(data.frame(
      run=thisIndex,
      variate=variate,
      smpl=smpl,
      xy_as_predictor = xy_as_predictor,
      xy_imp = xy_imp,
      n_sel_vars = n_sel_vars,
      inAOA_ffs = AOA_stats_ffs$inside_AOA,
      inAOA = AOA_stats$inside_AOA,
      AOA_gains_ffs = AOA_stats_ffs$inside_AOA - AOA_stats$inside_AOA,
      RMSE_gains_ffs = (true_RMSE - true_RMSE_ffs)/true_RMSE*100,
      RMSE_diff_ffs = diff_RMSE_ffs,
      RMSE_diff_ffs_aoa = diff_RMSE_ffs_inAOA,
      RMSE_diff = diff_RMSE,
      RMSE_diff_aoa = diff_RMSE_inAOA), file=f_out)
    } 
}


# ************ CALL THE FUNCTIONS ************ 
mclapply(list("AGB", "OCS"), function(x) {
  for(smpl in samples) {
    i <- thisIndex
    knndmCV(smpl, i, x, startseed)
  }
}, mc.cores = 2)
