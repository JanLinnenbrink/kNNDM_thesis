train_model = function(modelname, training_samples, predictors, response, folds, hyperparameter){
  
  if(!dir.exists(modelname)) dir.create(modelname)
  training_samples = training_samples %>% 
    dplyr::select(all_of(c(predictors, response))) %>% 
    st_drop_geometry()
  
  i = fold2index(folds)
  
  
  set.seed(7353)
  rfmodel = caret::train(x = training_samples %>% dplyr::select(all_of(predictors)),
                         y = training_samples %>% dplyr::pull(response),
                         method = "ranger",
                         tuneGrid = hyperparameter,
                         num.trees = 300,
                         trControl = trainControl(method = "cv", number = length(unique(folds)),
                                                  index = i$index, indexOut = i$indexOut,
                                                  savePredictions = "final"),
                         importance = "impurity")
  
  saveRDS(rfmodel, paste0(modelname, "/rfmodel_new.RDS"))
  return(rfmodel)
  
}

fold2index = function(fold){
  
  fold = data.frame(fold = fold)
  
  indOut = fold %>% dplyr::group_by(fold) %>%
    attr('groups') %>% dplyr::pull(.rows)
  
  ind = purrr::map(seq(length(indOut)), function(x){
    s = seq(nrow(fold))
    s = s[!s %in% indOut[[x]]]
    return(s)
  })
  return(
    list(
      index = ind,
      indexOut = indOut
    )
  )
  
}