# Libraries and utils ----
library("tidyverse")
library("raster")
library("terra")
library("CAST")
library("sf")
library("caret")
library("gstat")
library("virtualspecies")
library("parallel")
library("doParallel")
library("pbapply")
source("code/sim_utils.R")

# No need for proj4 warnings
options("rgdal_show_exportToProj4_warnings"="none")

#' Sample simulation: virtual species.
#' @details
#' Simulates a series of sampling points for the simulation.
#' @param nsamples Integer. Number of samples to simulate.
#' @param dsamples Character. Spatial distribution of the samples. 5 are
#' possible: sregular, wregular, random, wclust, sclust.
#' @param sarea sf/sfc polygon where samples will be simulated.
sim2_samples <- function(nsamples, dsamples, sarea){


  if(dsamples=="sregular"){
    simpoints <- jitterreg_sample(sarea, nsamples, 40000)
  }else if(dsamples=="wregular"){
    simpoints <- jitterreg_sample(sarea, nsamples, 80000)
  }else if(dsamples=="random"){
    simpoints <- st_sample(sarea, nsamples)
  }else if(dsamples=="wclust"){
    simpoints <- clustered_sample(sarea, nsamples, 25, 80000)
  }else if(dsamples=="sclust"){
    simpoints <- clustered_sample(sarea, nsamples, 10, 80000)
  }

  simpoints <- st_sf(geometry=simpoints)
  simpoints
}


#' Simulation function.
#' @details
#' The function takes a virtual species simulated landscape for the Iberian
#' peninsula using bioclim data, simulates sampling points
#' and fits a RF and evaluates it using spatial 10-fold CV, random 10-fold CV,
#' kNNDM 10-fold CV and NNDM LOO CV; as well as the true error.
#' @param rgrid sf or sfc point object. Prediction grid.
#' @param rstack SpatRaster object. Contains the landscape data and the
#' simulated outcome.
#' @param sampling_area sf or sfc polygon object. Sampling area.
#' @param sample_dist String or vector or string. Distribution of the sampling
#' points. 5 are possible: "sregular", "wregular", "random", "wclust","sclust".
sim_species <- function(rgrid, rstack, sampling_area,
                        dist_it=c("sregular", "wregular", "random",
                                      "wclust","sclust")){
  
  # sample prediction points from the prediction area
  ppoints <- st_sample(sampling_area, 1000, type="regular")

  # Initiate results object and fixed information for all models
  res <- data.frame()
  grid_data <- as.data.frame(terra::extract(rstack, terra::vect(rgrid)))
  form <- as.formula(paste0("outcome~", paste0("bio", 1:19, collapse="+")))
  pgrid<- data.frame(mtry=6)
  res <- data.frame()
  
  # Simulate sampling points according to parameters and constraints
  train_points <- sim2_samples(100, dist_it, sampling_area) |> st_transform(st_crs(ppoints))
  
  # Get training and surface data for modelling and validation
  train_data <- terra::extract(rstack, train_points)
  
  # cluster points based on kmeans
  coords <- sf::st_coordinates(train_points) |> as.matrix()
  clust <- stats::kmeans(x=coords, centers = 2)
  pts_id <- cbind(train_points, clust$cluster)
  pts_id$rand <- sample(c(1:2),nrow(pts_id), replace = TRUE)
  
  #### Folds based on outcome range
  # Estimate outcome range
  train_points$outcome <- train_data$outcome
  # Define folds based on kNNDM
  folds_kndm <- CAST::knndm(train_points, ppoints = ppoints, k=2, maxp=0.8)
  # Define folds based on NNDM
  folds_ndm <- CAST::nndm(train_points, ppoints = ppoints, min_train =  0.5)
  
  
  # Store results of the iteration
  res_it <- list(ppoints=ppoints, train_points=train_points, 
                 folds_random=pts_id$rand, folds_kndm=folds_kndm$clusters, 
                 folds_ndm=folds_ndm, folds_spat=clust$cluster,
                 kout=folds_kndm)
  

  res_it
}
