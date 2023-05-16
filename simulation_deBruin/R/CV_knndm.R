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
n_CV   <- 100  # number of cross validation replications
n_samp <- 100  # number of sample replicates (for each design)
cores <- 20

# create outfolders if they don't exist
if(!dir.exists(outfolder))
  dir.create(outfolder)

if(!dir.exists(paste0(outfolder, "/knndm_5k")))
  dir.create(paste0(outfolder, "/knndm_5k"))


# ************ FUNCTIONS ***************

sumSquares <- function(ref, pred){
  muref <- mean(ref, na.rm=T)
  SSR <- sum((ref - pred)^2)
  SST <- sum((ref - muref)^2)
  return(c(SSR, SST))
}

knndmCV <- function(smpl, number, variate, seed){
  
  fname  <-  paste0(variate, "_", smpl, sprintf("%03d", number), ".Rdata")
  f_out  <- file.path(outfolder,"knndm_5k", fname)
  
  if(!file.exists(f_out)) {
  fname <- paste0(variate, "data", sprintf("%03d", number), ".Rdata")
  f_in <- file.path(infolder,smpl,fname)
  load(f_in)
  
  # load ppoints
  load(file.path(infolder, "ppoints.Rdata"))
  
  MEC=RMSE=time=WS=time_mod <- numeric(n_CV)
  
  if(variate == "AGB"){
    pts_df <- data.frame(x=AGBdata$xcoord * 1000, y=AGBdata$ycoord * 1000)
    pts_sf <- st_as_sf(pts_df, coords = c("x","y"), crs = st_crs(ppoints))
  } else {
    pts_df <- data.frame(x=OCSdata$xcoord * 1000, y=OCSdata$ycoord * 1000)
    pts_sf <- st_as_sf(pts_df, coords = c("x","y"), crs = st_crs(ppoints))
  }
  
  n <- length(pts_df$x)
  
  for(i_CV in 1:n_CV){
    
    SSR <- 0
    SST <- 0
    
    set.seed(seed)
    time[i_CV] <- system.time(knndm <- knndm(pts_sf, ppoints = ppoints, k = 4, maxp = 0.8))[[3]]
    WS[i_CV] <- knndm$W
    fold <- knndm$clusters
    
    set.seed(seed)
    for(k in 1:length(unique(fold))){
      if(variate == "AGB"){
        RFmodel <- ranger(agb~., AGBdata[fold != k,], 
                          respect.unordered.factors=TRUE)
        refs <- AGBdata$agb[fold == k] 
        preds  <- predict(RFmodel, AGBdata[fold == k,])$predictions
      } else{
        RFmodel <- ranger(ocs~., OCSdata[fold != k,], 
                          respect.unordered.factors=TRUE)
        refs <- OCSdata$ocs[fold == k] 
        preds  <- predict(RFmodel, OCSdata[fold == k,])$predictions
      }
      
      
      # rmse function (global_Val)
      squares <- sumSquares(refs, preds)
      SSR <- SSR + squares[1]
      SST <- SST + squares[2]
    }
    
    MEC[i_CV]  <- 1 - SSR/SST
    RMSE[i_CV] <- sqrt(SSR/n)
    seed <- seed + 1
  } # loop over i_CV
  
  save(MEC, RMSE, time, time_mod, WS, file=f_out)
} 
}


# ************ CALL THE FUNCTIONS ************ 
mclapply(seq(n_samp), function(i) {
  for(smpl in samples) {
    knndmCV(smpl, i, "AGB", startseed)
    knndmCV(smpl, i, "OCS", startseed)
  }
}, mc.cores = cores)

