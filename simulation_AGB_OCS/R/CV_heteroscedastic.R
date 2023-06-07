# *****************************************************************************
# R Script implementing heteroscedastic model-based based cross-validation.  
# Related to the paper "Dealing with clustered samples for assessing map 
# accuracy by cross-validation".
# Contact: Sytze de Bruin, Wageningen University, Laboratory of Geo-information
# Science and Remote Sensing, email: sytze.debruin@wur.nl
# May 3, 2022
# *****************************************************************************
Sys.sleep(round(runif(1, min = 1, max = 240)))
# ****** load required libraries *******
.libPaths("/home/j/jlinnenb/r_packages")
library(sf)
library(gstat)
library(terra)
library(ranger)
library(spatstat)
library(parallel)

# ************ GLOBALS ***************

infolder1 <- "/scratch/tmp/jlinnenb/deBruin_add_nndm/CVresults/modelbased"
infolder2 <- "/scratch/tmp/jlinnenb/deBruin_add_nndm/CVresults/random"
infolder3 <- "/scratch/tmp/jlinnenb/deBruin_add_nndm/data"
outfolder <- "/scratch/tmp/jlinnenb/deBruin_add_nndm/CVresults/heteroscedastic"
nsim <- 200        # number of sequential Gaussian simulations
i_CV <- 1:3  # cross validation replications analysed
startseed <- 1234567
cores=1

csv_file <- file.path(outfolder, "runs.csv")
runs <- read.csv(csv_file)
lastIndex <- runs[nrow(runs),1]
thisIndex <- lastIndex + 1
print(paste0("this Index is: ", thisIndex))
runs[thisIndex,1] <- thisIndex
write.csv(runs, file = csv_file, row.names = FALSE)


# coordinate reference system
CRSlaea   <- paste0("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 ",
                    "+ellps=GRS80 +units=m +no_defs")

# check whether infolders exist; if not, stop
if(!dir.exists(infolder1)){
  cat('First run "CV_model_based.R"\n\n')
  stop(paste("directory", infolder1, "does not exist")) 
}

if(!dir.exists(infolder2)){
  cat('First run "CV_random.R"\n\n')
  stop(paste("directory", infolder2, "does not exist")) 
}

# create outfolder if needed
if(!dir.exists(outfolder))
  dir.create(outfolder, recursive=T)

# ************ FUNCTIONS ***************

getDensity <- function(x, y, win, rsl){
  spp <- ppp(x, y, win)
  s <- bw.CvL(spp)
  den <- density.ppp(spp, eps=rsl, sigma=s, positive=T)
  denmat <- as.matrix.im(den)
  denrot <- transmat(denmat, from="spatstat", to="Europe")
  denext <- ext(round(c(den$xrange, den$yrange),0))
  denrast<- rast(denrot, crs=CRSlaea)
  ext(denrast) <- denext
  return(denrast)
}

getVgm <- function(pts, variate){
  bnds <- c(1:5 * 5000, 3:20 * 10000)
  if(variate=="AGB"){
    vgmod <- vgm(0.15, "Exp", 40000, add.to = vgm(0.25, "Sph", 15000, 0.6))
  } else {
    vgmod <- vgm(0.1, "Sph", 20000, 0.6, add.to=vgm(0.2, "Gau", 170000))
  }
  vg  <- variogram(resprime~1, pts, boundaries=bnds)
  vgmod <- tryCatch(fit.variogram(vg, vgmod), warning = function(w) w)
  if (is(vgmod, "warning")){
    vgmod <- vgm(0.4, "Exp", 40000, 0.6)
    vgmod <- fit.variogram(vg, vgmod)
  }
  return(list(vg, vgmod))
}


# ************ MAIN ************ 

studarea <- st_union(st_read(file.path(infolder3, "strata.shp")))
polbuf <- st_buffer(studarea, 5000)
polbuf <- st_buffer(polbuf,  -2500)
win    <- as.owin(polbuf)
rm(studarea, polbuf)

# simulation grid (fixed)
agg_area <- rast(file.path(infolder3, "aggArea.tif"))
agg_pnts <- as.points(agg_area)
sf_pnts <- st_as_sf(agg_pnts)
st_crs(sf_pnts) <- as.character(NA)
st_crs(sf_pnts) <- CRSlaea
rm(agg_area)

# dataframe with covariates on the grid
AGBstack <- AGBstack <- rast(file.path(infolder3, "AGBstack.tif"))
COVdata <- extract(AGBstack, agg_pnts)
COVdata$ID <- NULL
COVdata$agb <- NULL
rm(agg_pnts, AGBstack)

nanID <- which(apply(COVdata, 1, function(x) any(is.na(x))))
COVdata <- COVdata[-nanID,]
sf_pnts <- sf_pnts[-nanID,]

# find files with all design realizations
f_ins <- list.files(infolder2, glob2rx("ptsAGB_clusterMedium???.Rdata")) # formerly pts???
f_ins <- list(f_ins, list.files(infolder2, glob2rx("ptsAGB_clusterStrong???.Rdata")))
f_ins <- list(f_ins, list.files(infolder2, glob2rx("ptsAGB_regular???.Rdata")))

f_ins <- f_ins[c(thisIndex, thisIndex+100, thisIndex+200)]
# loop over all files

mclapply(f_ins, function(f_in) {
  # retrieve strings for naming purposes
  lchar <- nchar(f_in)
  variate <- substr(f_in, 4, 6)
  design <- substr(f_in, 8, lchar-9)
  number <- substr(f_in, lchar-8, lchar-6)
  f_out  <- substr(f_in, 4, lchar)
  ref_in <- file.path("../samples", design,
                      paste0(variate, "data", number, ".Rdata"))
  
  # load CV data
  load(file.path(infolder2, f_in))  #pts_df
  
  # load variogram from model-based method for computing muOK
  load(file.path(infolder1, f_out))
  rm(RMSEs, MECs)
  vgsOK <- vgs; rm(vgs)
  
  # load reference data
  load(ref_in)
  if(variate=="AGB"){
    refs <- AGBdata$agb
    rm(AGBdata)
  } else{
    refs <- OCSdata$ocs
    rm(OCSdata)
  }
  
  # compute sampling intensity and point weights
  pts_df$x <- pts_df$x*1000
  pts_df$y <- pts_df$y*1000
  dens <- getDensity(pts_df$x, pts_df$y, win, 2500)

  pts_df$dens <- extract(dens, pts_df[, 1:2])[,2]
  
  pts_df <- st_as_sf(x=pts_df, coords=c("x", "y"))
  st_crs(pts_df) <- CRSlaea
  
  sf_pnts$dens <- extract(dens, vect(sf_pnts))[,2]
  
  qs <- quantile(pts_df$dens, 1:99 * 0.01)
  
  pts_df$intvl <- findInterval(pts_df$dens, qs)
  avgdens <- aggregate(pts_df$dens, list(pts_df$intvl), mean)[,2]
  
  bnds <- c(1:5 * 5000, 3:20 * 10000)
  
  # Predict variate at prediction grid locations
  # get sample data
  fsamp <- paste0(variate, "data", number, ".Rdata")
  fsamp <- file.path("../samples", design, fsamp)
  load(fsamp)
  
  # fit RF on the entire sample
  set.seed(startseed)
  if(variate == "AGB"){
    RFmodel <- ranger(agb~., AGBdata, respect.unordered.factors=TRUE)
  } else{
    RFmodel <- ranger(ocs~., OCSdata, respect.unordered.factors=TRUE)
  }
  
  # predict on grid
  preds <- predict(RFmodel, COVdata)[[1]]
  
  # run the approach for a selection of folds set by i_CV
  RMSEs <- numeric()
  MECs <- numeric()
  vgs <- list()
  # loessmods <- list()
  # loessdata <- list()
  ilist <- 1
  
  for(icv in i_CV){
    
    sdresid <- aggregate(as.data.frame(pts_df)[,icv], list(pts_df$intvl), sd)[,2]
  
    loessmod <- loess(sdresid~avgdens, degree=0, span=0.5, surface="direct")
    # loessdata[[ilist]] <- data.frame(sdresid=sdresid, avgdens=avgdens)
    # loessmods[[ilist]] <- loessmod
    
    # compute muOK
    vgmOK <- vgsOK[[ilist]]$vgmod
    pts_df$z <- as.data.frame(pts_df)[,icv]
    muOK <- mean(krige(z~1, pts_df, newdata=sf_pnts, nmax=75,
                       model=vgsOK[[ilist]]$vgmod)$var1.pred, na.rm=T)
    
    # convert residual
    pts_df$resprime <- (as.data.frame(pts_df)[,icv] - muOK)/
      predict(loessmod, data.frame(avgdens=pts_df$dens))
    
    # predict transformed residual
    vg <- getVgm(pts_df, variate)
    
    set.seed(startseed)
    resmaps <- krige(resprime~1, pts_df, newdata=sf_pnts, nsim=nsim, 
                     model=vg[[2]], nmax=75)
    resmaps$geometry <- NULL  # just need a dataframe
    
    # back-transform
    idx <- which(is.na(sf_pnts$dens))
    if(length(idx) > 0){
      resmaps <- muOK + resmaps[-idx,] * 
        predict(loessmod, data.frame(avgdens=sf_pnts$dens[-idx]))
      
      refs <- resmaps + preds[-idx]
    } else{
      resmaps <- muOK + resmaps * 
        predict(loessmod, data.frame(avgdens=sf_pnts$dens))
      refs <- resmaps + preds
    }
    
    
    # accuracy metrics
    SST <- apply(refs, 2, function(x) sum((x - mean(x))^2))
    
    RMSEs <- rbind(RMSEs, apply(resmaps, 2, 
                                function(x) sqrt(mean(x^2))))
    MECs  <- rbind(MECs, 1 - apply(resmaps, 2, 
                                   function(x) sum(x^2, na.rm=T))/SST)
    vgs[ilist] <- list(list(vg = vg[[1]], vgmod = vg[[2]]))
    
    ilist <- ilist + 1
    
  }
  
  save(RMSEs, MECs, vgs, file = file.path(outfolder, f_out))

}, mc.cores=cores)

