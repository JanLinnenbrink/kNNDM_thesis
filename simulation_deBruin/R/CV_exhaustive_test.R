# *****************************************************************************
# R Script implementing exhaustive validation, which is used as the reference.  
# Related to the paper "Dealing with clustered samples for assessing map 
# accuracy by cross-validation".
# Contact: Sytze de Bruin, Wageningen University, Laboratory of Geo-information
# Science and Remote Sensing, email: sytze.debruin@wur.nl
# May 3, 2022
# *****************************************************************************
# ****** load required libraries *******
Sys.sleep(round(runif(1, min = 1, max = 240)))

.libPaths("/home/j/jlinnenb/r_packages/")
setwd("/scratch/tmp/jlinnenb/deBruin_add_nndm/")

library(caret)
library(CAST)
library(terra)
library(parallel)

# ************ GLOBALS ***************
samples   <- c("clusterMedium", "clusterStrong", "clusterGapped", "regular", 
               "simpleRandom")
# infolder1 <- "../data"
# infolder2 <- "../samples"
# outfolder <- "../CVresults"
infolder1 <- "./data"
infolder2 <- "./samples"
outfolder <- "./CVresults"
startseed <- 1234567
n_samp    <- 100  # number of sample replicates (for each design)


# create outfolders if they don't exist
if(!dir.exists(outfolder))
  dir.create(outfolder)

if(!dir.exists(paste0(outfolder, "/exhaustive_new")))
  dir.create(paste0(outfolder, "/exhaustive_new"))

if(!file.exists(file.path(outfolder, "exhaustive_new", "runs.csv"))) {
  write.csv(data.frame("runs"=1), file.path(outfolder, "exhaustive_new", "runs.csv"), row.names = FALSE)
}

csv_file <- file.path(outfolder, "exhaustive_new", "runs.csv")
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

AGBstack$glc2017 <- terra::as.factor(AGBstack$glc2017)
OCSstack$glc2017 <- terra::as.factor(OCSstack$glc2017)

OCt <- rast(file.path(infolder1, "ocs.tif"))
AGt <- rast(file.path(infolder1, "agb.tif"))
OCS <- mask(OCt, msk, filename="tmpocs.tif", overwrite=T)
AGB <- mask(AGt, msk, filename="tmpagb.tif", overwrite=T)
rm(AGt, OCt)

# thisIndex <- 100
# ************ FUNCTIONS ***************

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
  f_out <- file.path(outfolder, "exhaustive_new", fname2)
  
  load(f_in)
  
  set.seed(seed)
  
  if(variate == "AGB"){
    mtry <- 4
    pgrid <- expand.grid(mtry=mtry, splitrule="variance", min.node.size=5)
    RFmodel <- caret::train(agb~., AGBdata,respect.unordered.factors=TRUE, method = "ranger",
                            tuneGrid=pgrid, num.trees=500)
    map  <- terra::predict(AGBstack, RFmodel, filename=f_out, overwrite=TRUE, na.rm=TRUE)
    ME   <- mefu(AGB, map)
    RMSE <- rmsefu(AGB, map)
    MEC  <- mecfu(AGB, map)
  } else {
    mtry <- 4
    pgrid <- expand.grid(mtry=mtry, splitrule="variance", min.node.size=5)
    RFmodel <- caret::train(ocs~., OCSdata,respect.unordered.factors=TRUE, method = "ranger",
                                 tuneGrid=pgrid, num.trees=500)
    map  <- terra::predict(OCSstack, RFmodel, filename=f_out, overwrite=TRUE, na.rm=TRUE)
    ME   <- mefu(OCS, map)
    RMSE <- rmsefu(OCS, map)
    MEC  <- mecfu(OCS, map)
  }
  
  fname <-  paste0(variate, "_", smpl, sprintf("%03d", number), ".Rdata")
  f_out2 <- file.path(outfolder,"exhaustive_new", fname)
  save(MEC, ME, RMSE, file=f_out2)
  file.remove(f_out)
}


# ************ CALL THE FUNCTIONS ************ 
mclapply(list("AGB", "OCS"), function(x) {
  for(smpl in samples) {
    i <- 1
    exhaustive(smpl, i, x, startseed)
  }
}, mc.cores = 2)

