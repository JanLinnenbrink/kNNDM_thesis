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
library(parallel)

# ************ GLOBALS ***************
samples   <- c("clusterMedium", "clusterStrong", "clusterGapped", "regular", 
               "simpleRandom")

infolder <- "/scratch/tmp/jlinnenb/deBruin_add_nndm/samples"
outfolder <- "/scratch/tmp/jlinnenb/deBruin_add_nndm/CVresults"
startseed <- 1234567
n_CV   <- 3  # number of cross validation replications
n_samp <- 100  # number of sample replicates (for each design)
cores <- 20

# create outfolders if they don't exist
if(!dir.exists(outfolder))
  dir.create(outfolder)

if(!dir.exists(paste0(outfolder, "/random")))
  dir.create(paste0(outfolder, "/random"))


# ************ FUNCTIONS ***************

sumSquares <- function(ref, pred){
  muref <- mean(ref, na.rm=T)
  SSR <- sum((ref - pred)^2)
  SST <- sum((ref - muref)^2)
  return(c(SSR, SST))
}

randomCV <- function(smpl, number, variate, seed){
  
  fname <- paste0(variate, "data", sprintf("%03d", number), ".Rdata")
  f_in <- file.path(infolder,smpl,fname)
  load(f_in)
  
  MEC  <- numeric(n_CV)
  RMSE <- numeric(n_CV)
  
  if(variate == "AGB"){
    pts_df <- data.frame(x=AGBdata$xcoord, y=AGBdata$ycoord)
  } else {
    pts_df <- data.frame(x=OCSdata$xcoord, y=OCSdata$ycoord)
  }
  
  n <- length(pts_df$x)
  
  for(i_CV in 1:n_CV){
    
    SSR <- 0
    SST <- 0
    resids <- numeric(n)
    
    # ************************
    # ****** 10 fold CV ******
    # ************************
    set.seed(seed)
    fold <- sample(rep(1:10, ceiling(n/10)))[1:n]
    
    set.seed(seed)
    for(k in 1:10){
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
      
      resids[fold == k] <- refs - preds
      squares <- sumSquares(refs, preds)
      SSR <- SSR + squares[1]
      SST <- SST + squares[2]
    }
    
    MEC[i_CV]  <- 1 - SSR/SST
    RMSE[i_CV] <- CAST::global_validation(RFmodel)["RMSE"][[1]]
    seed <- seed + 1
    pts_df  <- cbind(pts_df, tmp=resids)
    names(pts_df)[i_CV + 2] <- paste0("res", sprintf("%03d", i_CV))
    
  } # loop over i_CV
  
  fname  <-  paste0(variate, "_", smpl, sprintf("%03d", number), ".Rdata")
  fname2 <-  paste0("pts", variate, "_", smpl, 
                    sprintf("%03d", number), ".Rdata")
  f_out  <- file.path(outfolder,"random", fname)
  f_out2 <- file.path(outfolder,"random", fname2)
  save(MEC, RMSE, file=f_out)
  save(pts_df, file=f_out2)
}


# ************ CALL THE FUNCTIONS ************ 
mclapply(seq(n_samp), function(i) {
  for(smpl in samples) {
    randomCV(smpl, i, "AGB", startseed)
    randomCV(smpl, i, "OCS", startseed)
  }
}, mc.cores = cores)
