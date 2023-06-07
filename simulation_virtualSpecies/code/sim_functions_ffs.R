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
source("./code/sim_utils.R")

# No need for proj4 warnings
options("rgdal_show_exportToProj4_warnings"="none")


# helper function: calculate AOA percentages
AOA_perc <- function(model, indexTrain=NULL, indexTest=NULL, predictor_stack) {
  AOA_res <- aoa(predictor_stack, model=model, CVtrain = indexTrain, CVtest = indexTest)
  AOA <- AOA_res$AOA
  AOAn <- terra::freq(AOA)[,2:3]
  overall_cells <- sum(AOAn$count)
  AOAp <- AOAn$count / overall_cells * 100
  if(class(AOAp)!="list") list(data.frame(inside_AOA = 1, outside_AOA = 0), AOA, AOA_res)
  else list(data.frame(inside_AOA = AOAp[[2]], outside_AOA = AOAp[[1]]), AOA, AOA_res)
  
}

# helper function: calculate true error inside AOA
AOA_err <- function(AOA_raster, predictor_stack, outcome_rast, mod) {
  rstack <- c(predictor_stack, outcome_rast)
  rstack[AOA_raster<1] <- NA
  wgrid_curr <- st_as_sf(as.points(rstack))
  surfdf_curr <- terra::extract(rstack, terra::vect(wgrid_curr), ID=FALSE)
  surfdf_curr$preds <- predict(mod, newdata=surfdf_curr)
  surf_curr_stats <- surfdf_curr %>%
    summarise(RMSE_AOA = sqrt(mean((outcome-surfdf_curr$preds)^2)),
              MAE_AOA = mean(abs(outcome-surfdf_curr$preds)),
              R2_AOA = cor(outcome, surfdf_curr$preds)^2)
  names(surf_curr_stats) <- paste0(names(surf_curr_stats), "_surf")
  surf_curr_stats
}

# helper function: calculate true error for the whole AOI
surf_err <- function(mod, surfdf, what) {
  surfdf$preds <- predict(mod, newdata=surfdf)
  surf_stats <- surfdf %>%
    summarise(RMSE = sqrt(mean((outcome-preds)^2)),
              MAE = mean(abs(outcome-preds)),
              R2 = cor(outcome, preds)^2)
  names(surf_stats) <- paste0(names(surf_stats), "_surf_", what)
  surf_stats
}



calc_diff <- function(cntrl, what, traindf, predictor_names,outcome_rast,
                      prand, predictor_stack, rstack, surfdf) {
  
  quiet(mod <- ffs(traindf[,predictor_names], traindf$outcome, method="rf",
                     trControl=cntrl, tunerand=prand, ntree=100,
                   tuneLength=1))
  err_stats_spat <- global_validation(mod)
  err_stats_spat <- data.frame(
    RMSE = err_stats_spat[[1]],
    MAE = err_stats_spat[[3]],
    R2 = err_stats_spat[[2]])
  AOA_spat <- AOA_perc(mod, predictor_stack=predictor_stack)
  AOA_err <- AOA_err(AOA_spat[[2]], predictor_stack, outcome_rast= outcome_rast, mod=mod)

  whole_AOI_err <- surf_err(mod, surfdf, what)

  diff_CV_AOI <- err_stats_spat - whole_AOI_err
  diff_CV_AOA <- err_stats_spat - AOA_err
  names(diff_CV_AOI) <- paste0("diff_CV_AOI_", names(diff_CV_AOI))
  names(diff_CV_AOA) <- paste0("diff_CV_AOA_", names(diff_CV_AOA))

  AOA_cal <- calibrate_aoa(AOA_spat[[3]], mod, showPlot = FALSE)$AOA$expected_RMSE
  if(class(AOA_cal) != "SpatRaster") {
    AOA_cal <- rast(AOA_cal)
  }
  RMSE_caoa <- global(AOA_cal, mean, na.rm=TRUE)[[1]]
  diff_calAOA_RMSE <- RMSE_caoa - AOA_err$RMSE_AOA_surf
  
  form <- as.formula(paste0(paste0("outcome~", paste0("bio", 1:19, collapse="+")), " + lon", " + lat")) #1
  model_wo <- train(form, data=traindf, method="rf",
                    trControl=cntrl, tunespat=prand, ntree=100)
  
  true_rmse_ffs <- whole_AOI_err$RMSE
  true_rmse_wo <- surf_err(model_wo, surfdf, what)$RMSE
  improvement_rmse <- (true_rmse_wo-true_rmse_ffs)/true_rmse_wo * 100

  AOA_wo <- AOA_perc(model_wo, predictor_stack=predictor_stack)

  AOA_improvement <- (AOA_spat[[1]]$inside_AOA-AOA_wo[[1]]$inside_AOA)/AOA_wo[[1]]$inside_AOA * 100
  n_predictors_ffs <- nrow(mod$finalModel$importance)
  xy_as_predictor <- ifelse(any(rownames(mod$finalModel$importance) %in% c("lat", "lon")),
                            "yes","no")
  
  ffs_stats <- data.frame(n_predictors=n_predictors_ffs, 
                          xy_as_predictor=xy_as_predictor,
                          improvement_rmse=improvement_rmse,
                          improvement_aoa=AOA_improvement)
  
  out_stats <- cbind(diff_CV_AOI, diff_CV_AOA, diff_calAOA_RMSE, RMSE_caoa, ffs_stats)
  out_stats$in_AOA <- AOA_spat[[1]]$inside_AOA
  out_stats$in_AOA_wo <- AOA_wo[[1]]$inside_AOA
  names(out_stats) <- paste0(names(out_stats), "_", what)
  out_stats
}

#' Sample simulation 2: virtual species.
#' @details
#' Simulates a series of sampling points for simulation problem 2.
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
  }else if(dsamples=="vclust"){
    simpoints <- clustered_sample(sarea, nsamples, 10, 40000)
  }else if(dsamples=="eclust"){
    simpoints <- clustered_sample(sarea, nsamples, 5, 40000)
  }
  
  simpoints <- st_sf(geometry=simpoints)
  simpoints
}


#' Fits a RF model and evaluates it using several methods (species simulation).
#' @details
#' Fits a RF model and evaluates it using LOO CV, bLOO CV, NNDM LOO CV and true errors.
#' @param form String. Model formula.
#' @param folds_rand list. Indices for rand CV
#' training data.
#' @param folds_kndm list. Indices for kndm CV
#' data.
#' @param ndmout_indexTrain list. Indices for NNDM LOO CV (outcome range) training
#' data.
#' @param ndmout_indexTest list. Indices for NNDM LOO CV (outcome range) test data.
#' @param ndmres_indexTrain list. Indices for NNDM LOO CV (residual range) training
#' data.
#' @param ndmres_indexTest list. Indices for NNDM LOO CV (residual range) test data.
#' @param prand Data frame. Parameter rand of the model.
#' @param traindf Data frame. Training data to fit the model.
#' @param surfdf Data frame. Surface data.
fitval_rf_species <- function(form,
                              clust_indexTrain,
                              kndm_indexTrain,
                              kndm_indexTest,
                              nndm_indexTrain,
                              nndm_indexTest,
                              prand, traindf,
                              surfdf,train_points,
                              rstack, outcome) {
  
  predictor_names <- names(traindf)[names(traindf)!= "outcome"] 
  predictor_stack <- rstack[[predictor_names]]
  outcome_rast <- rstack$outcome
  
  
  # Validate with random CV and compute metrics
  r_cntrl <- trainControl(method="CV", savePredictions=TRUE)
  rand_stats <- calc_diff(r_cntrl, "rand",  traindf, predictor_names,
                          outcome_rast, prand, 
                          predictor_stack, rstack, surfdf)
  
  # Validate with 10-fold cluster CV and compute CV statistics
  spat_cntrl <- trainControl(method="cv", index=clust_indexTrain, savePredictions=TRUE)
  spat_stats <- calc_diff(spat_cntrl, "spat", traindf, predictor_names,
                          outcome_rast, prand, 
                          predictor_stack, rstack, surfdf)
  
  # Validate with knndm and compute CV statistics
  kndm_cntrl <- trainControl(method="cv",
                             index=kndm_indexTrain,
                             indexOut=kndm_indexTest,
                             savePredictions=TRUE)
  kndm_stats <- calc_diff(kndm_cntrl, "kndm", traindf, predictor_names,
                          outcome_rast, prand, 
                          predictor_stack, rstack, surfdf)
  
  
  # Validate with NNDM LOO and compute CV statistics (outcome range)
  nndm_cntrl <- trainControl(method="cv",
                             index=nndm_indexTrain,
                             indexOut=nndm_indexTest,
                             savePredictions=TRUE)
  nndm_stats <- calc_diff(nndm_cntrl, "nndm", traindf, predictor_names,
                          outcome_rast, prand, 
                          predictor_stack, rstack, surfdf)
  
  # Tidy and return results
  data.frame(rand_stats,spat_stats,kndm_stats,nndm_stats)
}


#' Simulation function 2: virtual species.
#' @details
#' The function takes a virtual species simulated landscape for the Iberian
#' peninsula using bioclim data, simulates sampling points, computes outcome
#' and residual autocorrelation range, and fits a RF and evaluates it using LOO,
#' bLOO, and NNDM LOO CV; as well as the true error.
#' @param wgrid sf or sfc point object. Prediction rand.
#' @param rstack stack raster object. Contains the landscape data and the
#' simulated outcome.
#' @param sampling_area sf or sfc polygon object. Sampling area.
#' @param sample_dist String or vector or string. Distribution of the sampling
#' points. 5 are possible: "sregular", "wregular", "random", "wclust","sclust".
sim_species <- function(wgrid, rstack, sampling_area,
                        sample_dist=c("sregular", "wregular", "random",
                                      "wclust","sclust")){

  
  ppoints <- st_sample(sampling_area, 1000, type="regular")
  
  # Initiate results object and fixed information for all models
  res <- data.frame()
  rand_data <- terra::extract(rstack, terra::vect(wgrid), ID=FALSE)
  rand_data <- rand_data[complete.cases(rand_data),]
  form <- as.formula(paste0("outcome~", paste0("bio", 1:19, collapse="+")))
  prand <- data.frame(mtry=6)
  
  i <- 0
  # Start sampling loop
  for(dist_it in sample_dist){
    
    # Simulate sampling points according to parameters and constraints
    train_points <- sim2_samples(100, dist_it, sampling_area) |> st_transform(st_crs(ppoints))
    
    # Get training and surface data for modelling and validation
    train_data <- terra::extract(rstack, train_points, ID=FALSE)
    
    # cluster points based on kmeans
    coords <- sf::st_coordinates(train_points) |> as.matrix()
    clust <- stats::kmeans(x=coords, centers = 10)
    pts_id <- cbind(train_points, clust$cluster)
    pts_id$rand <- sample(c(1:10),nrow(pts_id), replace = TRUE)
    
    #### Folds based on outcome range
    # Estimate outcome range
    train_points$outcome <- train_data$outcome
    # Define folds based on 5-fold cluster CV
    folds_clust <- CAST::CreateSpacetimeFolds(pts_id, spacevar="clust.cluster",k=10)
    # Define folds based on kNNDM
    folds_kndm <- CAST::knndm(train_points, ppoints = ppoints)
    # Define folds based on NNDM
    folds_ndm_out <- CAST::nndm(train_points, ppoints = ppoints, min_train =  0.5)
    #### Model fitting and validation
    mod <- fitval_rf_species(form,
                             folds_clust$index,
                             folds_kndm$indx_train,
                             folds_kndm$indx_test,
                             folds_ndm_out$indx_train,
                             folds_ndm_out$indx_test,
                             prand, train_data,
                             rand_data,train_points,rstack, train_data$outcome)
    
    # Store results of the iteration
    res_it <- cbind(data.frame(dsample=dist_it, stringsAsFactors = FALSE),
                    mod)
    res <- bind_rows(res, res_it)
  }
  
  row.names(res) <- NULL
  res
}
