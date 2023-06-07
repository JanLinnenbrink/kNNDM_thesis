# *****************************************************************************
# R Script implementing conventional random f-fold cross-validation.  
# Related to the paper "Dealing with clustered samples for assessing map 
# accuracy by cross-validation".
# Contact: Sytze de Bruin, Wageningen University, Laboratory of Geo-information
# Science and Remote Sensing, email: sytze.debruin@wur.nl
# May 3, 2022
# *****************************************************************************

# ****** load required library *******
.libPaths("/home/j/jlinnenb/r_packages/")
library(ranger)
library(CAST)
library(sf)
library(raster)
library(caret)
library(parallel)

# ************ GLOBALS ***************
setwd("/scratch/tmp/jlinnenb/deBruin_add_nndm/")
#setwd("C:/0_Msc_Loek/Z_Palma/deBruin_add_nndm/")

samples   <- c("clusterMedium", "clusterStrong", "clusterGapped", "regular", 
               "simpleRandom")

infolder <- "./samples"
outfolder <- "./CVresults"
startseed <- 1234567
n_CV   <- 3  # number of cross validation replications
n_samp <- 100  # number of sample replicates (for each design)
cores <- 30

# create outfolders if they don't exist
if(!dir.exists(outfolder))
  dir.create(outfolder)

if(!dir.exists(paste0(outfolder, "/nndm_caret")))
  dir.create(paste0(outfolder, "/nndm_caret"))


# ************ FUNCTIONS ***************

sumSquares <- function(ref, pred){
  SSR <- sum((ref - pred)^2)
  return(SSR)
  # cannot calculate muref, since only 1 prediction (loocv)
}

nndmCV <- function(smpl, number, variate, seed){
  fname  <-  paste0(variate, "_", smpl, sprintf("%03d", number), ".Rdata")
  f_out  <- file.path(outfolder,"nndm_caret", fname)
  
  if(!file.exists(f_out)) {
    fname <- paste0(variate, "data", sprintf("%03d", number), ".Rdata")
    f_in <- file.path(infolder,smpl,fname)
    load(f_in)
    
    # load ppoints
    load(file.path(infolder, "ppoints.Rdata"))
    
    RMSE=time=WS=time_mod <- numeric(n_CV)
    
    if(variate == "AGB"){
      pts_df <- data.frame(x=AGBdata$xcoord * 1000, y=AGBdata$ycoord * 1000)
      pts_sf <- st_as_sf(pts_df, coords = c("x","y"), crs = st_crs(ppoints)) #### only first 20 rows
    } else {
      pts_df <- data.frame(x=OCSdata$xcoord * 1000, y=OCSdata$ycoord * 1000)
      pts_sf <- st_as_sf(pts_df, coords = c("x","y"), crs = st_crs(ppoints))
    }
    
    
    n <- length(pts_df$x)
    
    for(i_CV in 1:n_CV){
      
      set.seed(seed)
      time[i_CV] <- system.time(nndm <- nndm(pts_sf, ppoints = ppoints))[[3]]
      WS[i_CV] <- twosamples::wass_stat(nndm$Gjstar, nndm$Gij)
      
      trControl <- trainControl(method = "cv", savePredictions = "final",
                                index=nndm$indx_train, indexOut = nndm$indx_test)
      
      set.seed(seed)
      if(variate == "AGB"){
        AGBdata$glc2017 <- as.factor(AGBdata$glc2017)
        pgrid <- expand.grid(splitrule="variance",min.node.size=5,mtry=4)
        time_mod[i_CV] <- system.time(RFmodel <- caret::train(agb~., AGBdata,
                                                              respect.unordered.factors=TRUE, 
                                                              method = "ranger",
                                                              tuneGrid=pgrid, num.trees=500,
                                                              trControl = trControl))[[3]]
      } else{
        OCSdata$glc2017 <- as.factor(OCSdata$glc2017)
        pgrid <- expand.grid(splitrule="variance",min.node.size=5,mtry=4)
        time_mod[i_CV] <- system.time(RFmodel <- caret::train(ocs~., OCSdata,
                                                              respect.unordered.factors=TRUE, 
                                                              method = "ranger",
                                                              tuneGrid=pgrid, num.trees=500,
                                                              trControl = trControl))[[3]]
      }
      
      RMSE[i_CV] <- global_validation(RFmodel)[[1]]
      
      seed <- seed + 1
      print(RMSE)
    } # loop over i_CV
    
    save(RMSE, time, time_mod, WS, file=f_out)
  }
  
}

# ************ CALL THE FUNCTIONS ************ 
mclapply(seq(n_samp), function(i) {
  for(smpl in samples) {
    nndmCV(smpl, i, "AGB", startseed)
    nndmCV(smpl, i, "OCS", startseed)
  }
}, mc.cores = cores)

