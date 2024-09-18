##### Model fits using bootstrapping
##### Daniel Dempsey

### Read in libraries
library( xgboost ) # Fit xgboost model
library( randomForest ) # Fit random forest model
library( mice ) # Missing data imputation
library( SHAPforxgboost ) # For shapley value plots
library( treeshap ) # To compute shapley values
library( readr ) # Load in csv files
library( readxl ) # Load in excel files
library( dplyr ) # Data frame manipulation tools
library( tidyr ) # Data frame manipulation tools
library( pROC ) # Computing ROC curves and AUC
library( parallel ) # Parallel computation
library( tictoc ) # For timing code
library( caret ) # Cross validation and random forest tools
library( data.table ) # Used when computing Shapley values
library( ggplot2 ) # Plotting functionality
library( ggridges ) # For making joyplots/ridgeplots
mycol <- c( '#56B4E9', '#E69F00' )

### Function to make directories if they don't already exist
mkdir <- function( x ) {
  if( !dir.exists( x ) ) {
    dir.create( x )
  }
}

### Loss function for cross validation (with hacky fix to avoid bugs)
log_loss <- function( pred, true ) {
  pred <- ifelse( pred == 0, pred + .Machine$double.eps, pred )
  pred <- ifelse( pred == 1, pred - .Machine$double.eps, pred )
  -sum(true*log(pred) + (1-true)*log(1-pred))   
}


### Alternate version of randomForest.unify that fixes a bug that occurs for classification problems
randomForest2.unify <- function (rf_model, data) {
  if (!inherits(rf_model, "randomForest")) {
    stop("Object rf_model was not of class \"randomForest\"")
  }
  if (any(attr(rf_model$terms, "dataClasses") != "numeric")) {
    stop("Models built on data with categorical features are not supported - please encode them before training.")
  }
  n <- rf_model$ntree
  ret <- data.table()
  x <- lapply(1:n, function(tree) {
    tree_data <- as.data.table(randomForest::getTree(rf_model, 
                                                     k = tree, labelVar = TRUE))
    tree_data[, c("left daughter", "right daughter", "split var", 
                  "split point", "prediction")]
  })
  times_vec <- sapply(x, nrow)
  y <- rbindlist(x)
  y$prediction <- as.numeric( y$prediction ) # <-------- this is the fix
  y[, `:=`(Tree, rep(0:(n - 1), times = times_vec))]
  y[, `:=`(Node, unlist(lapply(times_vec, function(x) 0:(x - 1))))]
  setnames(y, c("Yes", "No", "Feature", "Split", "Prediction", 
                "Tree", "Node"))
  y[, `:=`(Feature, as.character(Feature))]
  y[, `:=`(Yes, Yes - 1)]
  y[, `:=`(No, No - 1)]
  y[y$Yes < 0, "Yes"] <- NA
  y[y$No < 0, "No"] <- NA
  y[, `:=`(Missing, NA)]
  y[, `:=`(Missing, as.integer(Missing))]
  ID <- paste0(y$Node, "-", y$Tree)
  y$Yes <- match(paste0(y$Yes, "-", y$Tree), ID)
  y$No <- match(paste0(y$No, "-", y$Tree), ID)
  y$Cover <- 0
  y$Decision.type <- factor(x = rep("<=", times = nrow(y)), 
                            levels = c("<=", "<"))
  y[is.na(Feature), `:=`(Decision.type, NA)]
  y[!is.na(Feature), `:=`(Prediction, NA)]
  y[is.na(Feature), `:=`(Prediction, Prediction/n)]
  setcolorder(y, c("Tree", "Node", "Feature", "Decision.type", 
                   "Split", "Yes", "No", "Missing", "Prediction", "Cover"))
  ret <- list(model = as.data.frame(y), data = as.data.frame(data))
  class(ret) <- "model_unified"
  attr(ret, "missing_support") <- FALSE
  attr(ret, "model") <- "randomForest"
  return(set_reference_dataset(ret, as.data.frame(data)))
}


### Fitting function
# Model mode = 1 corresponds to XGBoost, anything else is Random Forest
fit_fun <- function( X, y, param_grid, model_mode = 1, num_folds = 10, B = 1000, seed = 1 ) {
  
  set.seed( seed )
  
  # set working directory
  if ( model_mode == 1 ) { directory <- 'XGBoost' }
  else { directory <- 'RandomForest' }
  mkdir( directory )
  setwd( directory )
  
  # Save parameter grid
  grid_save <- as.data.frame( do.call( 'cbind', param_grid ) )
  write_csv( grid_save, 'parameter_grid.csv' )
  
  # Prepare xgboost sets
  if ( model_mode == 1 ) {
    xgb_dat <- xgb.DMatrix( data = X, label = y )
  }
  
  # Prepare parameter grid for cross-validation
  nrounds_ind <- which( names(param_grid) == 'nrounds' )
  all_param_grid <- expand.grid( param_grid )
  all_param_grid <- split( all_param_grid, 1:nrow(all_param_grid) )
  
  # Cross validation settings
  cv_fold_inds <- createFolds( y, k = num_folds )
  
  # Run cross validation
  ce_list <- auc_list <- list()
  for ( i in 1:num_folds ) {
    
    # Select cross-validation set
    X_cv_train <- X[-cv_fold_inds[[i]], ]
    y_cv_train <- y[-cv_fold_inds[[i]]]
    
    X_cv_test <- X[cv_fold_inds[[i]], ]
    y_cv_test <- y[cv_fold_inds[[i]]]
    
    if (model_mode == 1) {
      
      cv_train_xgb <- xgb.DMatrix( data = X_cv_train, label = y_cv_train )
      cv_test_xgb <- xgb.DMatrix( data = X_cv_test, label = y_cv_test )
      
      cv_fit <- lapply( all_param_grid, function(x) {
        xgb.train( params = as.list(x[-nrounds_ind]), data = cv_train_xgb, 
                   nrounds = x[[nrounds_ind]], verbose = 0 )
      } )
      
      test_preds <- lapply( cv_fit, predict, newdata = cv_test_xgb )
      
    }
    else {
      
      cv_fit <- lapply( all_param_grid, function(z) {
        randomForest( x = X_cv_train, y = factor(y_cv_train),
                      ntree = z$ntree, sampsize = z$sampsize,
                      nodesize = z$nodesize, mtry = z$mtry ) 
      } )
      
      test_preds <- lapply( cv_fit, function(x) { predict( x, newdata = X_cv_test, type = 'prob' )[, 2] } )
      
    }
    
    auc_list[[i]] <- sapply( test_preds, function(x){ roc( y_cv_test, x, direction = '<', quiet = TRUE )$auc } )
    ce_list[[i]] <- sapply( test_preds, function(x){ log_loss(x, y_cv_test) } )
    
  }
  
  # Compile results; cross entropys and AUCs
  ce_dat <- do.call( 'rbind', ce_list )
  auc_dat <- do.call( 'rbind', auc_list )
  
  ce_mean <- apply( ce_dat, 2, mean )
  ce_se <- apply( ce_dat, 2, sd ) / sqrt( num_folds ) 
  auc_mean <- apply( auc_dat, 2, mean )
  auc_se <- apply( auc_dat, 2, sd ) / sqrt( num_folds )
  
  # Select best parameters based on whichever set minimised the CV error on average
  best_ind <- which.min( ce_mean )
  best_params <- all_param_grid[[best_ind]]
  
  # Create some diagnostic plots
  png( 'CV_CE.png' )
  cv_ce_plot <- ggplot( mapping = aes(x = 1:length(ce_mean), y = ce_mean) ) + geom_line() + 
    geom_vline( xintercept = best_ind, linetype = 'dashed', col = 'red' ) +
    ggtitle( '10-fold Cross Validation Cross-Entropy' ) + xlab( 'Index' ) + 
    ylab( '' )
  print( cv_ce_plot )
  dev.off()
  
  png( 'CV_AUC.png' )
  cv_auc_plot <- ggplot( mapping = aes(x = 1:length(auc_mean), y = auc_mean) ) + geom_line() + 
    geom_vline( xintercept = best_ind, linetype = 'dashed', col = 'red' ) +
    ggtitle( '10-fold Cross Validation AUC' ) + xlab( 'Index' ) + 
    ylab( '' )
  print( cv_auc_plot )
  dev.off()
  
  # Save cross-validation results
  write_csv( best_params, 'CV_chosen_parameters.csv' )
  write_csv( data.frame(mean = ce_mean, se = ce_se), file = 'CV_CE.csv' )
  write_csv( data.frame(mean = auc_mean, se = auc_se), file = 'CV_AUC.csv' )
  
  cv_results <- list( cv_err_dat = ce_dat, cv_auc_dat = auc_dat, 
                      chosen_parameters = best_params, best_param_ind = best_ind,
                      param_grid = all_param_grid, X = X, y = y )
  
  # Bootstrap settings and initialization
  train_prop <- 0.8
  ncol_x <- ncol( X )
  nrow_x <- nrow( X )
  dim_x <- ncol_x * nrow_x
  nms_x <- colnames( X )
  avg_shaps <- matrix( 0, nrow = B, ncol = ncol_x )
  all_auc <- matrix( 0, nrow = B, ncol = 2 )
  all_train_inds <- all_fits <- train_rocs <- test_rocs <- train_preds <- test_preds <- list()
  colnames( avg_shaps ) <- nms_x
  colnames( all_auc ) <- c( 'Train', 'Test' )
  train_num <- ceiling( nrow_x * train_prop )
  test_num <- nrow_x - train_num
  all_shaps <- data.frame( boot_ind = rep(1:B, each = nrow_x * ncol_x ), 
                           variable = rep(rep(nms_x, each = nrow_x), B),
                           obs_ind = rep(1:nrow_x, B * ncol_x),
                           shap = 0 )
  all_preds <- data.frame( boot_ind = rep(1:B, each = nrow_x),
                           obs_ind = 0, train_ind = rep(c(rep(1, train_num), rep(0, test_num)), B), 
                           pred = 0, true_val = 0 )
  
  # Run bootstrap
  for ( i in 1:B ) {
    
    # Split into training / test sets
    train_inds <- all_train_inds[[i]] <- createDataPartition( y, p = train_prop ) %>% unlist
    test_inds <- setdiff( 1:nrow_x, train_inds )
    
    preds_inds <- (((i-1)*nrow_x)+1):(i*nrow_x)
    all_preds$obs_ind[preds_inds] <- c( train_inds, test_inds )
    
    X_train <- X[train_inds, ]
    X_test <- X[test_inds, ]
    y_train <- y[train_inds]
    y_test <- y[test_inds]
    
    all_preds$true_val[preds_inds] <- c( y_train, y_test )
    
    if ( model_mode == 1 ) {
      # XGBoost fit
      train_xgb <- xgb.DMatrix( data = X_train, label = y_train )
      test_xgb <- xgb.DMatrix( data = X_test, label = y_test )
      fit <- all_fits[[i]] <- xgb.train( params = as.list(best_params[-nrounds_ind]), 
                                         data = train_xgb, nrounds = best_params[[nrounds_ind]] )
      preds_train <- predict( fit, newdata = X_train )
      preds_test <- predict( fit, newdata = X_test )
      fit_unify <- xgboost.unify( fit, X_train )
    }
    else {
      # Random forest fit
      fit <- all_fits[[i]] <- randomForest( x = X_train, y = factor(y_train), 
                                            ntree = best_params$ntree, sampsize = best_params$sampsize,
                                            nodesize = best_params$nodesize, mtry = best_params$mtry )
      preds_train <- predict( fit, newdata = X_train, type = 'prob' )[, 2]
      preds_test <- predict( fit, newdata = X_test, type = 'prob' )[, 2]
      fit_unify <- randomForest2.unify( fit, X_train )
    }
    
    # Predictions
    all_preds$pred[preds_inds] <- c( preds_train, preds_test )
    
    # Shapley values
    shap_X <- treeshap( fit_unify, X, verbose = FALSE )$shaps
    shap_train <- treeshap( fit_unify, X_train, verbose = FALSE )$shaps
    all_shaps$shap[(((i-1)*dim_x)+1):(i*dim_x)] <- unlist( shap_X )
    avg_shaps[i, ] <- apply( shap_train, 2, function(x) { mean(abs(x)) } )
    
    # AUCs
    train_rocs[[i]] <- roc( y_train, preds_train, direction = '<', quiet = TRUE )
    test_rocs[[i]] <- roc( y_test, preds_test, direction = '<', quiet = TRUE )
    all_auc[i, ] <- c( train = train_rocs[[i]]$auc,
                       test = test_rocs[[i]]$auc )
    
  }
  
  # Compile results
  shap_ranks <- t( apply( -avg_shaps, 1, rank ) )
  rank_medians <- apply( shap_ranks, 2, quantile, probs = 0.5 )
  
  # Predictions plot
  tru_val_lab <- c( 'Non-Exposed', 'Exposed' )
  train_val_lab <- c( 'Train', 'Test' )
  all_preds$true_lab <- factor( tru_val_lab[all_preds$true_val + 1], tru_val_lab )
  all_preds$train_lab <- factor( rev(train_val_lab)[all_preds$train_ind + 1], train_val_lab )
  
  png( 'Bootstrapped_Predictions.png' )
  all_preds_plot <- ggplot( all_preds, aes( true_lab, pred, fill = train_lab ) ) + 
    geom_boxplot( ) + scale_fill_manual( values = mycol, name = '' ) +
    ylim(0, 1) + ggtitle( 'Bootstraped Predictions' ) + ylab( 'Prediction' ) + xlab( '' )
  print( all_preds_plot )
  dev.off()
  
  png( 'Bootstrapped_AUC.png' )
  all_auc_long <- pivot_longer( as.data.frame(all_auc), cols = all_of(train_val_lab) )
  all_auc_long$name <- factor( all_auc_long$name, levels = train_val_lab )
  boot_auc <- ggplot( all_auc_long, aes(name, value) ) + geom_boxplot( fill = mycol ) + ylim(0, 1) +
    geom_hline( yintercept = 0.5, linetype = 'dashed', col = 'darkgrey' ) + ylab( 'AUC' ) + xlab( '' ) +
    ggtitle( 'Bootstrapped AUC' )
  print( boot_auc )
  dev.off()
  
  # Shapley ranks plot
  png( 'Bootstrap_Shap_Ranks.png' )
  long_shap <- data.frame( var = rep(colnames(shap_ranks), each = B), 
                           rank = unlist(as.data.frame(shap_ranks)) )
  long_shap$var <- factor( long_shap$var, levels = colnames(shap_ranks)[order(rank_medians, decreasing = TRUE)] )
  shap_joy <- ggplot( long_shap, aes( x = rank, y = var, fill = factor(ifelse(var=="random_num","Highlighted","Normal")) ) ) + 
    scale_fill_manual(name = "var", values = rev(mycol)) + xlab( '' ) +
    geom_density_ridges( ) + xlim( c(0, ncol(shap_ranks)) ) + theme( legend.position = 'none' )
  print( shap_joy )
  dev.off()
  
  # Save and return results
  bootstrap_results <- list( all_fits = all_fits, auc = all_auc, avg_shaps = avg_shaps, 
                             all_shaps = all_shaps, shap_ranks = shap_ranks, 
                             all_preds = all_preds, all_train_inds = all_train_inds,
                             train_rocs = train_rocs, test_rocs = test_rocs )
  
  res <- list( cv_results = cv_results, bootstrap_results = bootstrap_results )
  save( res, file = 'res.Rdata' )
  setwd( '..' )
  
}

### Load and prepare data
bsi_dat <- read_csv( 'Data/bsi_dat.csv' )
X_full <- select( bsi_dat, -c(StudyNo, Exposure) ) %>% as.matrix
X <- complete( mice(X_full, printFlag = FALSE, seed = 42) ) %>% as.matrix
y <- bsi_dat$Exposure

### Prepare parameter grids
# Test parameter grid for debug purposes
#param_grid_xgb <- list( objective = 'binary:logistic', 
#                        nrounds = c(20, 50), max_depth = 1 )

# XGBoost
param_grid_xgb <- list( max_depth = c(1, 3, 5), eta = c(0.005, 0.01, 0.05, 0.1),
                        gamma = c(5, 10, 15), subsample = c(0.5, 0.75, 1), 
                        colsample_bytree = c(0.5, 0.75, 1),
                        lambda = c(0.1, 0.9), objective = 'binary:logistic', 
                        nrounds = c(10, 20, 50, 100) )

# Random forest
param_grid_rf <- list( ntree = c( 10, 50, 100 ), sampsize = c( 72, 92, 116 ),
                       mtry = c( 4, 6, 8 ), nodesize = c( 1, 3, 5 ) )

### Set working directory
home_dir <- 'Output/Model_Fits/'
mkdir( home_dir )
setwd( home_dir )

### Run model fit function
# XGBoost
tic()
fit_fun( X, y, model_mode = 1, param_grid = param_grid_xgb, seed = 42 ) # XGBoost
xgb_time <- toc() # 655.537 sec elapsed

# Random forest
tic()
fit_fun( X, y, model_mode = 2, param_grid = param_grid_rf, seed = 42 ) # Random forest
rf_time <- toc() # 268.588 sec elapsed

