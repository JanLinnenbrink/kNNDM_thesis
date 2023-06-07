# *****************************************************************************
# R Script implementing exhaustive validation, which is used as the reference.  
# Related to the paper "Dealing with clustered samples for assessing map 
# accuracy by cross-validation".
# Contact: Sytze de Bruin, Wageningen University, Laboratory of Geo-information
# Science and Remote Sensing, email: sytze.debruin@wur.nl
# May 3, 2022
# *****************************************************************************
#Sys.sleep(round(runif(1,min=1,max=240)))
# ****** load required libraries *******
#.libPaths("/home/j/jlinnenb/r_packages")
library(ranger)
library(terra)
library(parallel)
library(caret)
library(CAST)

# ************ GLOBALS ***************
samples   <- c("clusterMedium", "clusterStrong", "clusterGapped", "regular", 
               "simpleRandom")
infolder1 <- "~/CrossValidation/deBruin_add_nndm/data"
infolder2 <- "~/CrossValidation/deBruin_add_nndm/samples"
outfolder <- "~/CrossValidation/deBruin_add_nndm/CVresults"
startseed <- 1234567
n_samp    <- 100  # number of sample replicates (for each design)


# create outfolders if they don't exist
if(!dir.exists(outfolder))
  dir.create(outfolder)

if(!dir.exists(paste0(outfolder, "/exhaustive")))
  dir.create(paste0(outfolder, "/exhaustive"))

if(!file.exists(file.path(outfolder, "exhaustive", "runs.csv"))) {
  write.csv(data.frame("runs"=1), file.path(outfolder, "exhaustive", "runs.csv"), row.names = FALSE)
}

csv_file <- file.path(outfolder, "exhaustive", "runs.csv")
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

OCt <- rast(file.path(infolder1, "ocs.tif"))
AGt <- rast(file.path(infolder1, "agb.tif"))
OCS <- mask(OCt, msk, filename="tmpocs.tif", overwrite=T)
AGB <- mask(AGt, msk, filename="tmpagb.tif", overwrite=T)
rm(AGt, OCt)

# thisIndex <- 100
# ************ FUNCTIONS ***************
predfun <- function(object, newdata){
  pred <- predict(object, newdata)
  pred[[1]]
}

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

exhaustive <- function(smpl, number, variate, seed){
  
  fname1 <- paste0(variate, "data", sprintf("%03d", number), ".Rdata")
  fname2 <- paste0(variate, smpl, sprintf("%03d", number), ".tif")
  
  f_in  <- file.path(infolder2,smpl,fname1)
  f_out <- file.path(outfolder, "exhaustive", fname2)
  
  load(f_in)
  
  set.seed(seed)
  if(variate == "AGB"){
    
    mtry <- floor(sqrt(ncol(AGBdata[,-1])))
    pgrid <- expand.grid(splitrule="variance",min.node.size=5,mtry=mtry)
    RFmodel <- caret::train(agb~., AGBdata,respect.unordered.factors=TRUE, method = "ranger",
                                 tuneGrid=pgrid, num.trees=500,
                                 trControl = trainControl(method = "cv", savePredictions = "final"))
    
  } else{
    
    
    mtry <- floor(sqrt(ncol(OCSdata[,-1])))
    pgrid <- expand.grid(splitrule="variance",min.node.size=5,mtry=mtry)
    RFmodel <- caret::train(ocs~., OCSdata,respect.unordered.factors=TRUE, method = "ranger",
                                 tuneGrid=pgrid, ntree=500, 
                                 trControl = trainControl(method = "cv", savePredictions = "final"))
  }
  
  err_rand_caret <- global_validation(RFmodel)
  RMSE <- err_rand_caret[[1]]
  
  fname <-  paste0(variate, "_", smpl, sprintf("%03d", number), ".Rdata")
  f_out2 <- file.path(outfolder,"exhaustive", fname)
  save(RMSE, file=f_out2)
  file.remove(f_out)
}


# ************ CALL THE FUNCTIONS ************ 
mclapply(list("AGB", "OCS"), function(x) {
  for(smpl in samples) {
    i <- thisIndex
    exhaustive(smpl, i, x, startseed)
  }
}, mc.cores = 10)
