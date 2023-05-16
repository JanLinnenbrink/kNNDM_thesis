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
source("code/knndm_W.R")

# No need for proj4 warnings
options("rgdal_show_exportToProj4_warnings"="none")


#' Sample simulation: virtual species.
#' @details
#' Simulates a series of sampling points for the simulation.
#' Two more extreme sampling strategies to increase the range of different WS.dist.
#' @param nsamples Integer. Number of samples to simulate.
#' @param dsamples Character. Spatial distribution of the samples. 7 are
#' possible: "sregular", "wregular", "random", "wclust","sclust","vclust","eclust".
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
  }else if(dsamples=="vclust"){
    simpoints <- clustered_sample(sarea, nsamples, 10, 40000)
  }else if(dsamples=="eclust"){
    simpoints <- clustered_sample(sarea, nsamples, 5, 40000)
  }
  
  simpoints <- st_sf(geometry=simpoints)
  simpoints
}


#' Fits a RF model and evaluates it using different kNNDM configurations and sample point distributions.
#' @details
#' Fits a RF model and evaluates it using kNNDM 10-fold CV in different configurations and true errors.
#' @param form String. Model formula.
#' @param kndm_folds list. Indices for kndm CV.
#' @param pgrid Data frame. Parameter grid of the model.
#' @param traindf Data frame. Training data to fit the model.
#' @param surfdf Data frame. Surface data.
fitval_rf_species <- function(form,
                              kndm_folds,
                              pgrid, traindf,
                              surfdf) {
  
  # Validate with random CV and compute metrics
  r_cntrl <- trainControl(method="CV", savePredictions=TRUE)
  rand_mod <- train(form, data=traindf, method="rf",
                    trControl=r_cntrl, tuneGrid=pgrid, ntree=100)
  rand_stats <- rand_mod$pred %>%
    summarise(RMSE = sqrt(mean((obs-pred)^2)),
              MAE = mean(abs(obs-pred)),
              R2 = cor(obs, pred)^2)
  names(rand_stats) <- paste0(names(rand_stats), "_rand")
  
  # Compute CV statistics in surface
  surfdf$preds <- predict(rand_mod, newdata=surfdf)
  surf_stats <- surfdf %>%
    summarise(RMSE = sqrt(mean((outcome-preds)^2)),
              MAE = mean(abs(outcome-preds)),
              R2 = cor(outcome, preds)^2)
  names(surf_stats) <- paste0(names(surf_stats), "_surf")
  
  # Validate with knndm
  if (class(kndm_folds)=="list") {
    kndm_stats <- lapply(kndm_folds, function(x) {
      fdf <- data.frame(f=x)
      ctrl <- CAST::CreateSpacetimeFolds(fdf, spacevar=1, k = max(fdf))
      kndm_ctrl <- trainControl(method="cv",
                                index=ctrl$index,
                                indexOut=ctrl$indexOut,
                                savePredictions=TRUE)
      kndm_out_mod <- suppressWarnings( # train() can't compute R2
        train(form, data=traindf, method="rf",
              trControl=kndm_ctrl, tuneGrid=pgrid, ntree=100))
      kndm_stats <- kndm_out_mod$pred %>%
        summarise(RMSE = sqrt(mean((obs-pred)^2)),
                  MAE = mean(abs(obs-pred)),
                  R2 = cor(obs, pred)^2)
      names(kndm_stats) <- paste0(names(kndm_stats), "_kndm")
      kndm_stats
    })
    kndm_stats <- do.call(rbind, kndm_stats)
  } else {
    fdf <- data.frame(f=kndm_folds)
    ctrl <- CAST::CreateSpacetimeFolds(fdf, spacevar=1, k = max(fdf))
    kndm_ctrl <- trainControl(method="cv",
                              index=ctrl$index,
                              indexOut=ctrl$indexOut,
                              savePredictions=TRUE)
    kndm_out_mod <- suppressWarnings( # train() can't compute R2
      train(form, data=traindf, method="rf",
            trControl=kndm_ctrl, tuneGrid=pgrid, ntree=100))
    kndm_stats <- kndm_out_mod$pred %>%
      summarise(RMSE = sqrt(mean((obs-pred)^2)),
                MAE = mean(abs(obs-pred)),
                R2 = cor(obs, pred)^2)
    names(kndm_stats) <- paste0(names(kndm_stats), "_kndm")
    kndm_stats
  }
  
  kndm_stats_diff <- data.frame(
    RMSE_surf = surf_stats$RMSE_surf, R2_surf = surf_stats$R2_surf, MAE_surf = surf_stats$MAE_surf,
    RMSE_kndm = kndm_stats$RMSE_kndm, R2_kndm = kndm_stats$R2_kndm, MAE_kndm = kndm_stats$MAE_kndm)
  kndm_stats_diff
}


#' Simulation function to generate many possible configurations of W and CV error.
#' @details
#' The function takes a virtual species simulated landscape for the Iberian
#' peninsula using bioclim data, simulates sampling points and fits a RF and 
#' evaluates it using kNNDM CV, as well as the true error.
#' @param rgrid sf or sfc point object. Prediction grid.
#' @param rstack spatRast object. Contains the landscape data and the
#' simulated outcome.
#' @param sampling_area sf or sfc polygon object. Sampling area.
#' @param sample_dist String or vector or string. Distribution of the sampling
#' points. 7 are possible: "sregular", "wregular", "random", "wclust","sclust","vclust","eclust".
sim_species <- function(rgrid, rstack, sampling_area,
                        sample_dist=c("sregular", "wregular", "random", "wclust","sclust","vclust","eclust")){
  
  ppoints <- st_sample(sampling_area, 1000, type="regular")
  
  # Initiate results object and fixed information for all models
  res <- data.frame()
  grid_data <- as.data.frame(terra::extract(rstack, terra::vect(rgrid)))
  form <- as.formula(paste0("outcome~", paste0("bio", 1:19, collapse="+")))
  pgrid <- data.frame(mtry=6)
  
  i <- 0
  # Start sampling loop
  for(dist_it in sample_dist){
    i=i+1
    
    # Simulate sampling points according to parameters and constraints
    train_points <- sim2_samples(100, dist_it, sampling_area)
    ppoints <- sf::st_transform(ppoints, sf::st_crs(train_points))
    
    # Get training and surface data for modelling and validation
    train_data <- terra::extract(rstack, train_points)
    
    # Estimate outcome range
    train_points$outcome <- train_data$outcome
    
    # kndm
    folds_kndm <- knndmW(train_points, ppoints = ppoints, clustering = "kmeans", k=10, maxp=0.5)
    
    #### Model fitting and validation
    mod <- fitval_rf_species(form,
                             folds_kndm$clusters,
                             pgrid, train_data,
                             grid_data)
    mod_all <- cbind(mod, data.frame(WS=folds_kndm$W))
    
    # Store results of the iteration
    res_it <- cbind(data.frame(dsample=dist_it, stringsAsFactors = FALSE),
                    mod_all)
    res <- bind_rows(res, res_it)
  }
  
  row.names(res) <- NULL
  res
}
