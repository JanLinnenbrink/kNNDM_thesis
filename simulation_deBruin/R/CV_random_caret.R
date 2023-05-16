# *****************************************************************************
# R Script implementing conventional random f-fold cross-validation.  
# Related to the paper "Dealing with clustered samples for assessing map 
# accuracy by cross-validation".
# Contact: Sytze de Bruin, Wageningen University, Laboratory of Geo-information
# Science and Remote Sensing, email: sytze.debruin@wur.nl
# May 3, 2022
# *****************************************************************************

# ****** load required library *******
#.libPaths("~/r_packages/")
library(ranger)
library(CAST)
library(sf)
library(raster)
library(caret)
library(parallel)

#source("~/deBruin_add_nndm/R/fold2index.R")
#source("~/deBruin_add_nndm/R/global_validation.R")

#source("./R/fold2index.R")
#source("./R/global_validation.R")

# ************ GLOBALS ***************
#setwd("/scratch/tmp/jlinnenb/deBruin_add_nndm/")

samples   <- c("clusterMedium", "clusterStrong", "clusterGapped", "regular", 
               "simpleRandom")

infolder <- "./samples"
outfolder <- "./CVresults"
startseed <- 1234567
n_CV   <- 3  # number of cross validation replications
n_samp <- 100 # number of sample replicates (for each design)
cores <- 10

# create outfolders if they don't exist
if(!dir.exists(outfolder))
  dir.create(outfolder)

if(!dir.exists(paste0(outfolder, "/random_caret_new")))
  dir.create(paste0(outfolder, "/random_caret_new"))


# ************ FUNCTIONS ***************


knndmCV <- function(smpl, number, variate, seed){
  
  #smpl="clusterMedium"
  #number=1
  #variate="AGB"
  #seed=1
  #i_CV=1
  
  fname  <-  paste0(variate, "_", smpl, sprintf("%03d", number), ".Rdata")
  f_out  <- file.path(outfolder,"random_caret_new", fname)
  
  
  fname <- paste0(variate, "data", sprintf("%03d", number), ".Rdata")
  f_in <- file.path(infolder,smpl,fname)
  load(f_in)
  
  # load ppoints
  load(file.path(infolder, "ppoints.Rdata"))
  
  RMSE <- numeric(n_CV)
  
  if(variate == "AGB"){
    pts_df <- data.frame(x=AGBdata$xcoord * 1000, y=AGBdata$ycoord * 1000)
    pts_sf <- st_as_sf(pts_df, coords = c("x","y"), crs = st_crs(ppoints))
  } else {
    pts_df <- data.frame(x=OCSdata$xcoord * 1000, y=OCSdata$ycoord * 1000)
    pts_sf <- st_as_sf(pts_df, coords = c("x","y"), crs = st_crs(ppoints))
  }
  
  for(i_CV in 1:n_CV){
    
    # RMSE random caret
    set.seed(seed)
    if(variate == "AGB"){
      
      mtry <- floor(sqrt(ncol(AGBdata[,-1])))
      pgrid <- expand.grid(mtry=mtry, splitrule="variance", min.node.size=5)
      RFmodel_rand <- caret::train(agb~., AGBdata,respect.unordered.factors=TRUE, method = "ranger",
                                   tuneGrid=pgrid, num.trees=500,
                                   trControl = trainControl(method = "cv", savePredictions = "final"))
      
    } else{
      
      
      mtry <- floor(sqrt(ncol(OCSdata[,-1])))
      pgrid <- expand.grid(mtry=mtry, splitrule="variance", min.node.size=5)
      RFmodel_rand <- caret::train(ocs~., OCSdata,respect.unordered.factors=TRUE, method = "ranger",
                                   tuneGrid=pgrid, num.trees=500, 
                                   trControl = trainControl(method = "cv", savePredictions = "final"))
    }
    
    err_rand_caret <- global_validation(RFmodel_rand)
    RMSE[i_CV] <- err_rand_caret[[1]]
    
    
    } # loop over i_CV
  save(RMSE, file=f_out)
  
}


# ************ CALL THE FUNCTIONS ************ 
mclapply(seq(n_samp), function(i) {
  for(smpl in samples) {
    knndmCV(smpl, i, "AGB", startseed)
    knndmCV(smpl, i, "OCS", startseed)
  }
}, mc.cores = cores)

