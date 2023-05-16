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

samples   <- c("clusterMedium", "clusterStrong", "clusterGapped", "regular", 
               "simpleRandom")
infolder <- "./samples"
outfolder <- "./CVresults"
startseed <- 1234567
n_CV   <- 3  # number of cross validation replications
n_samp <- 100  # number of sample replicates (for each design)
cores <- 10

# create outfolders if they don't exist
if(!dir.exists(outfolder))
  dir.create(outfolder)

if(!dir.exists(paste0(outfolder, "/knndm_caret_new")))
  dir.create(paste0(outfolder, "/knndm_caret_new"))

# ************ FUNCTIONS ***************


knndmCV <- function(smpl, number, variate, seed){
  
  fname  <-  paste0(variate, "_", smpl, sprintf("%03d", number), ".Rdata")
  f_out  <- file.path(outfolder,"knndm_caret_new", fname)
  
  if(!file.exists(f_out)) {
    fname <- paste0(variate, "data", sprintf("%03d", number), ".Rdata")
    f_in <- file.path(infolder,smpl,fname)
    load(f_in)
    
    # load ppoints
    load(file.path(infolder, "ppoints.Rdata"))
    
    if(variate == "AGB"){
      pts_df <- data.frame(x=AGBdata$xcoord * 1000, y=AGBdata$ycoord * 1000)
      pts_sf <- st_as_sf(pts_df, coords = c("x","y"), crs = st_crs(ppoints))
    } else {
      pts_df <- data.frame(x=OCSdata$xcoord * 1000, y=OCSdata$ycoord * 1000)
      pts_sf <- st_as_sf(pts_df, coords = c("x","y"), crs = st_crs(ppoints))
    }
    
    n <- length(pts_df$x)
    RMSE=time=time_mod <- numeric(n_CV)
    
    for(i_CV in 1:n_CV){
      
      set.seed(seed)
      tune_grid <- data.frame(k=2:10, maxp=seq(0.8,0.5,length.out=9))
      kndm <- apply(tune_grid, 1, function(x) knndm(pts_sf, ppoints = ppoints, k=x[[1]], maxp=x[[2]]))
      W_min <- lapply(kndm, function(x) x$W) |> 
        which.min()
      knndm <- kndm[[W_min]]
      
      trControl = trainControl(method = "cv", savePredictions = "final",
                               index=knndm$indx_train)
      
      set.seed(seed)
      if(variate == "AGB"){
        
        AGBdata$glc2017 <- as.factor(AGBdata$glc2017)
        mtry <- 4
        pgrid <- expand.grid(splitrule="variance",min.node.size=5,mtry=mtry)
        time_mod[i_CV] <- system.time(RFmodel <- caret::train(agb~., AGBdata,
                                                        respect.unordered.factors=TRUE, 
                                                        method = "ranger",
                                tuneGrid=pgrid, num.trees=500,
                                trControl = trControl))[[3]]
        
      } else{
        
        
        mtry <- 4
        pgrid <- expand.grid(splitrule="variance",min.node.size=5,mtry=mtry)
        time_mod[i_CV] <- system.time(RFmodel <- caret::train(ocs~., OCSdata,
                                                              respect.unordered.factors=TRUE,
                                                              method = "ranger",
                                tuneGrid=pgrid, num.trees=500, 
                                trControl = trControl))[[3]]
      }
      
      
      RMSE[i_CV] <- global_validation(RFmodel)[[1]]
      
      seed <- seed + 1
    } # loop over i_CV
    
    save(RMSE, time, time_mod, file=f_out)
  } 
}


# ************ CALL THE FUNCTIONS ************ 
mclapply(seq(n_samp), function(i) {
  for(smpl in samples) {
    knndmCV(smpl, i, "AGB", startseed)
    knndmCV(smpl, i, "OCS", startseed)
  }
}, mc.cores = cores)

