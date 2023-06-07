# *****************************************************************************
# R Script implementing conventional random cross-validation.  
# Related to the manuscript "Dealing with clustered samples for assessing map 
# accuracy by cross-validation".
# Contact: Sytze de Bruin, Wageningen University, Laboratory of Geo-information
# Science and Remote Sensing, email: sytze.debruin@wur.nl
# May 3, 2022
# *****************************************************************************

# ****** load required libraries *******
.libPaths("/home/j/jlinnenb/R")
library(ranger)
library(sperrorest)
library(parallel)


# ************ GLOBALS ***************
samples   <- c("clusterMedium", "clusterStrong", "clusterGapped", "regular", 
               "simpleRandom")
infolder <- "/scratch/tmp/jlinnenb/deBruin_add_nndm/samples"
outfolder <- "/scratch/tmp/jlinnenb/deBruin_add_nndm/CVresults"
startseed <- 1234567
n_CV      <- 3  # number of cross validation replications
n_samp    <- 100  # # number of sample replicates (for each design)
cores <- 20

# create outfolders if they don't exist
if(!dir.exists(outfolder))
  dir.create(outfolder)

if(!dir.exists(paste0(outfolder, "/spatial")))
  dir.create(paste0(outfolder, "/spatial"))


# ************ FUNCTIONS ***************

predfun <- function(object, newdata){
  pred <- predict(object, newdata)
  pred[[1]]
}


spatialCV <- function(smpl, number, variate, seed){
  
  fname <- paste0(variate, "data", sprintf("%03d", number), ".Rdata")
  f_in <- file.path(infolder,smpl,fname)
  load(f_in)
  
  if(variate == "AGB"){
    fo <- as.formula(paste0("agb~", paste(names(AGBdata)[-1], collapse = "+")))
    tst <- sperrorest(fo, data=AGBdata, model_fun=ranger, 
                      model_args=list(respect.unordered.factors=TRUE,
                                      mtry=4, splitrule="variance", min.node.size=5),
                      pred_fun=predfun, 
                      smp_fun=partition_kmeans, coords=c("xcoord", "ycoord"),
                      smp_args=list(balancing_steps = 1, seed1=seed, 
                                    repetition=1:n_CV, iter.max = 50), 
                      err_fun = err_fu)
  } else{
    fo <- as.formula(paste0("ocs~", paste(names(OCSdata)[-1], collapse = "+")))
    tst <- sperrorest(fo, data=OCSdata, model_fun=ranger, 
                      model_args=list(respect.unordered.factors=TRUE,
                                      mtry=4, splitrule="variance", min.node.size=5),
                      pred_fun=predfun, 
                      smp_fun=partition_kmeans, coords=c("xcoord", "ycoord"),
                      smp_args=list(balancing_steps = 1, seed1=seed, 
                                    repetition=1:n_CV, iter.max = 50), 
                      err_fun = err_fu)
  }

  ME   <- tst$error_rep$test_me
  RMSE <- tst$error_rep$test_rmse
  MEC  <- tst$error_rep$test_mec
  rm(tst)
  
  fname <-  paste0(variate, "_", smpl, sprintf("%03d", number), ".Rdata")
  f_out <- file.path(outfolder, "spatial", fname)
  save(MEC, ME, RMSE, file=f_out)
  
}


# ************ CALL THE FUNCTIONS ************ 
# change order to get it done with express partition
mclapply(order(seq(n_samp), decreasing=TRUE), function(i) {
  for(smpl in samples) {
    spatialCV(smpl, i, "AGB", startseed)
    spatialCV(smpl, i, "OCS", startseed)
  }
}, mc.cores = cores)
