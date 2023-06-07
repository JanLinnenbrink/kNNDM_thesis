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


#' Fits a RF model and evaluates it using several methods (species simulation).
#' @details
#' Fits a RF model and evaluates it using random 10-fold CV, spatial 10-fold CV,
#' NNDM LOO CV, kNNDM 10-fold CV and true errors.
#' @param form String. Model formula.
#' @param spatial_index list. Indices for spatial CV.
#' @param kndm_indexTrain list. Indices for kNNDM 10-fold CV training data.
#' @param kndm_indexTest list. Indices for kNNDM 10-fold CV test data.
#' @param nndm_indexTrain list. Indices for NNDM LOO CV test data.
#' @param nndm_indexTest list. Indices for NNDM LOO CV training data.
#' @param pgrid Data frame. Parameter rand of the model.
#' @param traindf Data frame. Training data to fit the model.
#' @param surfdf Data frame. Surface data.
fitval_rf_species <- function(form,
                              spatial_index,
                              kndm_indexTrain,
                              kndm_indexTest,
                              nndm_indexTrain,
                              nndm_indexTest,
                              pgrid, traindf,
                              surfdf) {


  # Validate with random CV and compute metrics
  r_cntrl <- trainControl(method="CV", savePredictions=TRUE)
  rand_mod <- train(form, data=traindf, method="rf",
                   trControl=r_cntrl, tuneGrid=pgrid, ntree=100)
  err_stats_rand <- global_validation(rand_mod)
  rand_stats <- data.frame(RMSE = err_stats_rand[[1]],
              MAE = err_stats_rand[[3]],
              R2 = err_stats_rand[[2]])
  names(rand_stats) <- paste0(names(rand_stats), "_rand")

  # Compute CV statistics in surface
  surfdf$preds <- predict(rand_mod, newdata=surfdf)
  surf_stats <- surfdf %>%
    summarise(RMSE = sqrt(mean((outcome-preds)^2)),
              MAE = mean(abs(outcome-preds)),
              R2 = cor(outcome, preds)^2)
  names(surf_stats) <- paste0(names(surf_stats), "_surf")

  # Validate with 10-fold cluster CV and compute CV statistics
  spatial_cntrl <- trainControl(method="cv", index=spatial_index,
                                 savePredictions=TRUE)
  spatial_mod <- suppressWarnings( # train() can't compute R2
    train(form, data=traindf, method="rf",
          trControl=spatial_cntrl, tuneGrid=pgrid, ntree=100))
  err_stats_spatial <- global_validation(spatial_mod)
  spatial_stats <- data.frame(RMSE = err_stats_spatial[[1]],
              MAE = err_stats_spatial[[3]],
              R2 = err_stats_spatial[[2]])
  names(spatial_stats) <- paste0(names(spatial_stats), "_spatial")

  # Validate with knndm and compute CV statistics
  kndm_cntrl <- trainControl(method="cv",
                                index=kndm_indexTrain,
                                indexOut=kndm_indexTest,
                                savePredictions=TRUE)
    time_kndm_mod <- system.time(kndm_mod <- suppressWarnings( # train() can't compute R2
      train(form, data=traindf, method="rf",
          trControl=kndm_cntrl, tuneGrid=pgrid, ntree=100)))[[3]]
  err_stats_kndm <- global_validation(kndm_mod)
  kndm_stats <- data.frame(RMSE = err_stats_kndm[[1]],
              MAE = err_stats_kndm[[3]],
              R2 = err_stats_kndm[[2]])
  names(kndm_stats) <- paste0(names(kndm_stats), "_kndm")

  # Validate with NNDM LOO and compute CV statistics (outcome range)
  nndm_cntrl <- trainControl(method="cv",
                                index=nndm_indexTrain,
                                indexOut=nndm_indexTest,
                                savePredictions=TRUE)
  time_nndm_mod <- system.time(nndm_mod <- suppressWarnings( # train() can't compute R2
    train(form, data=traindf, method="rf",
          trControl=nndm_cntrl, tuneGrid=pgrid, ntree=100)))[[3]]
  err_stats_nndm <- global_validation(nndm_mod)
  nndm_stats <- data.frame(RMSE = err_stats_nndm[[1]],
              MAE = err_stats_nndm[[3]],
              R2 = err_stats_nndm[[2]])
  names(nndm_stats) <- paste0(names(nndm_stats), "_nndm")

  # Tidy and return results
  data.frame(rand_stats,surf_stats,spatial_stats,kndm_stats,nndm_stats, time_kndm_mod, time_nndm_mod)
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
                        sample_dist=c("sregular", "wregular", "random",
                                      "wclust","sclust")){
  
  # sample prediction points from the prediction area
  ppoints <- st_sample(sampling_area, 1000, type="regular")

  # Initiate results object and fixed information for all models
  res <- data.frame()
  grid_data <- as.data.frame(terra::extract(rstack, terra::vect(rgrid)))
  form <- as.formula(paste0("outcome~", paste0("bio", 1:19, collapse="+")))
  pgrid<- data.frame(mtry=6)
  res <- data.frame()
  
  i <- 0
  # Start sampling loop
  for(dist_it in sample_dist){
    i=i+1
    # Simulate sampling points according to parameters and constraints
    train_points <- sim2_samples(100, dist_it, sampling_area) |> st_transform(st_crs(ppoints))

    # Get training and surface data for modelling and validation
    train_data <- terra::extract(rstack, train_points)

    # cluster points based on kmeans
    coords <- sf::st_coordinates(train_points) |> as.matrix()
    clust <- stats::kmeans(x=coords, centers = 10)
    pts_id <- cbind(train_points, clust$cluster)
    pts_id$rand <- sample(c(1:10),nrow(pts_id), replace = TRUE)

    #### Folds based on outcome range
    # Estimate outcome range
    train_points$outcome <- train_data$outcome
    # Define random folds
    folds_random <- CAST::CreateSpacetimeFolds(pts_id, spacevar="rand", k=10)
    # Define folds based on 10-fold cluster CV
    folds_spatial <- CAST::CreateSpacetimeFolds(pts_id, spacevar="clust.cluster",k=10)
    # Define folds based on kNNDM
    time_kndm_alg <- system.time(folds_kndm <- CAST::knndm(train_points, ppoints = ppoints))[[3]]
    # Define folds based on NNDM
    time_nndm_alg <- system.time(folds_ndm <- CAST::nndm(train_points, 
                                                         ppoints = ppoints, min_train =  0.5))[[3]]
    #### Model fitting and validation
    mod <- fitval_rf_species(form,
                             folds_spatial$index,
                             folds_kndm$indx_train,
                             folds_kndm$indx_test,
                             folds_ndm$indx_train,
                             folds_ndm$indx_test,
                             pgrid, train_data,
                             grid_data)

    # Store results of the iteration
    res_it <- cbind(data.frame(dsample=dist_it, stringsAsFactors = FALSE), mod, 
                    time_kndm_alg, time_nndm_alg)
    res <- bind_rows(res, res_it)
  }

  row.names(res) <- NULL
  res
}
