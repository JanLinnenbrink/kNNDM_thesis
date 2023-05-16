# *****************************************************************************
# R Script implementing inverse sampling-intensity based cross-validation.  
# Related to the paper "Dealing with clustered samples for assessing map 
# accuracy by cross-validation".
# Contact: Sytze de Bruin, Wageningen University, Laboratory of Geo-information
# Science and Remote Sensing, email: sytze.debruin@wur.nl
# May 3, 2022
# *****************************************************************************

# ****** load required libraries *******
.libPaths("/home/j/jlinnenb/r_packages")
library(spatstat)
library(terra)
library(sf)
library(parallel)
library(caret)
library(CAST)

# ************ GLOBALS ***************

infolder1 <- "/scratch/tmp/jlinnenb/deBruin_add_nndm/data"
infolder2 <- "/scratch/tmp/jlinnenb/deBruin_add_nndm/CVresults/random"
outfolder <- "/scratch/tmp/jlinnenb/deBruin_add_nndm/CVresults/intensity_caret"

infolder1 <- "~/CrossValidation/deBruin_add_nndm/data"
infolder2 <- "~/CrossValidation/deBruin_add_nndm/CVresults/random"
outfolder <- "~/CrossValidation/deBruin_add_nndm/CVresults/intensity_caret"

CRSlaea   <- paste0("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 ",
                    "+ellps=GRS80 +units=m +no_defs")

n_CV      <- 3 # number of cross validation replications
cores     <- 10


# check whether infolder2 exists; if not, stop
if(!dir.exists(infolder2)){
  cat('First run "CV_random.R"\n\n')
  stop(paste("directory", infolder2, "does not exist")) 
}

# create outfolder if needed
if(!dir.exists(outfolder))
  dir.create(outfolder, recursive=T)


# ************ FUNCTIONS ***************

getDensity <- function(x, y, sf_pol, rsl){
  win <- as.owin(sf_pol)
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

weighted.mec <- function(ref, resid, w=rep(1, length(ref))){
  muref <- weighted.mean(ref, w, na.rm=T)
  SSR_w <- sum(w * resid^2, na.rm=T)
  SST_w <- sum(w * (ref - muref)^2, na.rm=T)
  1 - SSR_w/SST_w
}

weighted.rmse <- function(resid, w=rep(1, length(resid))){
  sqrt(sum(resid^2*w, na.rm=T)/sum(w))
}

# ************ CALL THE FUNCTION ************ 
studarea <- st_union(st_read(file.path(infolder1, "strata.shp")))
polbuf <- st_buffer(studarea, 5000)
polbuf <- st_buffer(polbuf,  -2500)
rm(studarea)

f_ins <- list.files(infolder2, glob2rx("pts*.Rdata"))

mclapply(f_ins, function(f_in) {
  # retrieve strings for naming purposes
  lchar <- nchar(f_in)
  f_out <- substr(f_in, 4, lchar)
  variate <- substr(f_in, 4, 6)
  design <- substr(f_in, 8, lchar-9)
  number <- substr(f_in, lchar-8, lchar-6)
  ref_in <- file.path("~/CrossValidation/deBruin_add_nndm/samples", design,
                      paste0(variate, "data", number, ".Rdata"))
  
  # load CV data
  load(file.path(infolder2, f_in))  #pts_df
  
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
  dens <- getDensity(pts_df$x*1000, pts_df$y*1000,
                     polbuf, 2500)
  w <- 1/(extract(dens, pts_df[, 1:2]*1000)[,2])
  
  # compute metrics
  MEC  <- numeric(n_CV)
  RMSE <- numeric(n_CV)
  
  for(i in 1:n_CV){
    MEC[i]  <- weighted.mec(refs, pts_df[,i+2], w) # pts_df[,i+2]
    RMSE[i] <- weighted.rmse(pts_df[,i+2], w) # same
  }
  
  # save output
  save(MEC, RMSE, file=file.path(outfolder, f_out))
  return()
}, mc.cores = cores)
