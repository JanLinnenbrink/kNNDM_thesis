library(ranger)
library(CAST)
library(sf)
library(raster)
library(caret)
library(parallel)


samples   <- c("clusterMedium", "clusterStrong", "clusterGapped", "regular", 
               "simpleRandom")

infolder <- "~/CrossValidation/deBruin_add_nndm/samples"
outfolder <- "~/CrossValidation/deBruin_add_nndm/CVresults"
startseed <- 1234567
n_samp <- 100  # number of sample replicates (for each design)


sumSquares <- function(ref, pred){
  muref <- mean(ref, na.rm=T)
  SSR <- sum((ref - pred)^2)
  SST <- sum((ref - muref)^2)
  return(c(SSR, SST))
}


randCV <- function(smpl, number, seed){
  
  # load training data
  fname <- paste0("AGB", "data", sprintf("%03d", number), ".Rdata")
  f_in <- file.path(infolder,smpl,fname)
  load(f_in)
  load(file.path(infolder, "ppoints.Rdata"))
  
  pts_df <- data.frame(x=AGBdata$xcoord * 1000, y=AGBdata$ycoord * 1000)
  pts_sf <- st_as_sf(pts_df, coords = c("x","y"), crs = st_crs(ppoints))
    
  # manual CV
  n <- length(pts_df$x)
  SSR <- 0
  SST <- 0
  
  fold <- sample(1:10, n, replace=TRUE)
  
  set.seed(seed)
  mtry <- round(sqrt(ncol(AGBdata[,-1])))

  pgrid <- expand.grid(mtry=mtry, splitrule="variance", min.node.size=5)
  
  RFmodel_rand <- caret::train(agb~., AGBdata, respect.unordered.factors=TRUE, method = "ranger",
                               tuneGrid=pgrid, num.trees=500,
                               trControl = trainControl(method = "cv", savePredictions = "final"))
    

  err_rand_caret <- global_validation(RFmodel_rand)
  RMSE_caret <- err_rand_caret[[1]]
  
  # load truth
  load(paste0("~/CrossValidation/deBruin_add_nndm/CVresults/exhaustive/AGB_simpleRandom", 
              stringr::str_pad(number, 3, pad = "0"), ".Rdata"))
  RMSE_diff_caret <- (RMSE_caret-RMSE)/RMSE * 100
  
  return(c(RMSE_diff_caret))
  seed <- seed + 1
  
}

i=1:100
smpl="simpleRandom"
diff <- lapply(i, function(j) {
   randCV(smpl, j, startseed)
})

#diff4 <- lapply(i, function(j) {
#  randCV(samples, j, startseed, mtry=4)
#})

diff=unlist(diff)
diff=data.frame(diff=diff, method="caret")

ggplot(data=diff, aes(x=method, y=diff)) +
  geom_boxplot() +
  geom_hline(yintercept=0)

# exhaustive sample: 1) mit mtry = 4 2) mit caret