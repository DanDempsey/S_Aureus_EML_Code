##### Code for Further Graphs and Analysis
##### Daniel Dempsey

### Load libraries
library( pROC ) # For computing ROC curves
library( Gmedian ) # For estimating the geometric median
library( tidyr ) # For data frame manipulation tools
library( dplyr ) # For data frame manipulation tools
library( SHAPforxgboost ) # For shapley value plots
library( data.table ) # For creating data table - required for the shapley computations
library( ggplot2 ) # Extra plotting functionality
library( ggridges ) # For creating ridge plots
library( ggforce ) # Used in the creation of violin plots
library( ggpubr ) # For placing multiple graphs on the same plot
library( readr ) # For reading in csv files
setwd( 'Output/Model_Fits/XGBoost/' )
mycol <- c( '#56B4E9', '#E69F00' )

### Function to make directories if they don't already exist
mkdir <- function( x ) {
  if( !dir.exists( x ) ) {
    dir.create( x )
  }
}

### Load in XGBoost results
load( 'res.Rdata' )

### Geometric Median of bootstrapped Shapley values
# Construct a 2x2 matrix 
shap_mat <- pivot_wider( res$bootstrap_results$all_shaps, id_cols = 'boot_ind', 
                         names_from = c('obs_ind', 'variable'), values_from = 'shap' )[, -1]

# Compute geometric median
med_shap <- Gmedian( shap_mat )

# Plot geometric median
X <- res$cv_results$X
med_shap_list <- split( t(med_shap), rep( 1:37, each = nrow(X) ) )
colnames(X)[which(colnames(X) == 'AetiologyofCKD.4')] <- 'Diabetes'
names( med_shap_list ) <- colnames( X )
shap_contrib <- do.call( 'cbind', med_shap_list ) %>% as.data.table

png( 'Geometric_Median_Shapley_Values_full.png' )
shap_long <- shap.prep( shap_contrib = shap_contrib, X_train = X )
shap_long_split <- split( shap_long, shap_long$variable )
shap_long$stdfvalue <- unlist( lapply( shap_long_split, function(x) { (rank(x$rfvalue)-1) / (nrow(x)-1) } ) )
s_plot <- shap.plot.summary( shap_long )
print( s_plot )
dev.off()

# Plot a smaller version
shap.plot.summary2 <- function (data_long, x_bound = NULL, dilute = FALSE, scientific = FALSE, 
          my_format = NULL) 
{
  if (scientific) {
    label_format = "%.1e"
  }
  else {
    label_format = "%.3f"
  }
  if (!is.null(my_format)) 
    label_format <- my_format
  N_features <- setDT(data_long)[, uniqueN(variable)]
  if (is.null(dilute)) 
    dilute = FALSE
  nrow_X <- nrow(data_long)/N_features
  if (dilute != 0) {
    dilute <- ceiling(min(nrow_X/10, abs(as.numeric(dilute))))
    set.seed(1234)
    data_long <- data_long[sample(nrow(data_long), min(nrow(data_long)/dilute, 
                                                       nrow(data_long)/2))]
  }
  x_bound <- if (is.null(x_bound)) 
    max(abs(data_long$value)) * 1.1
  else as.numeric(abs(x_bound))
  plot1 <- ggplot(data = data_long) + coord_flip(ylim = c(-x_bound, 
                                                          x_bound)) + geom_hline(yintercept = 0) + ggforce::geom_sina(aes(x = variable, 
                                                                                                                          y = value, color = stdfvalue), method = "counts", maxwidth = 0.7, 
                                                                                                                      alpha = 0.7) + geom_text(data = unique(data_long[, c("variable", 
                                                                                                                                                                           "mean_value")]), aes(x = variable, y = -Inf, label = sprintf(label_format, 
                                                                                                                                                                                                                                        mean_value)), size = 6, alpha = 0.7, hjust = -0.2, fontface = "bold") + 
    scale_color_gradient(low = "#FFCC33", high = "#6600CC", 
                         breaks = c(0, 1), labels = c(" Low", "High "), guide = guide_colorbar(barwidth = 12, 
                                                                                               barheight = 0.3)) + theme_bw() + theme(axis.line.y = element_blank(), 
                                                                                                                                      axis.ticks.y = element_blank(), legend.position = "bottom", 
                                                                                                                                      legend.title = element_text(size = 10), legend.text = element_text(size = 8), 
                                                                                                                                      axis.title.x = element_text(size = 10)) + scale_x_discrete(limits = rev(levels(data_long$variable)), 
                                                                                                                                                                                                 labels = SHAPforxgboost:::label.feature(rev(levels(data_long$variable)))) + 
    labs(y = "SHAP value (impact on model output)", x = "", 
         color = "Feature value  ")
  return(plot1)
}
cutoff <- min( which( shap_long$variable == 'random_num' ) ) - 1
shap_long_trunc <- shap_long[1:cutoff, ]
shap_long_trunc$variable <- factor( shap_long_trunc$variable )
s_plot_trunc <- shap.plot.summary2( shap_long_trunc ) + theme( text = element_text(size = 20),
                                                              legend.text = element_text(size = 12),
                                                              legend.title = element_text(size = 14),
                                                              axis.title.x = element_text(size = 20),
                                                              )
### Re-draw bootstrapped Shapley value ranks to better fit the paper
vars <- unique(shap_long_trunc$variable) %>% as.vector %>% rev
shap_ranks <- res$bootstrap_results$shap_ranks
B <- length( res$bootstrap_results$all_fits )
shap_ranks_long <- data.frame( var = rep(colnames(shap_ranks), each = B), 
                               rank = unlist(as.data.frame(shap_ranks)) ) %>%
  filter( var %in% vars )

shap_ranks_long$var <- factor( shap_ranks_long$var, levels = vars )
shap_joy <- ggplot( shap_ranks_long, aes( x = rank, y = var ) ) + xlab( '' ) +
  geom_density_ridges( fill = mycol[1] ) + xlim( c(0, ncol(shap_ranks)) ) + 
  theme( legend.position = 'none', text = element_text(size = 20) ) + xlim( c(-1, 28) ) + ylab( '' )

shap_box <- ggplot( shap_ranks_long, aes( x = rank, y = var ) ) + xlab( '' ) +
  geom_boxplot( fill = mycol[1] ) + xlim( c(0, ncol(shap_ranks)) ) + 
  theme( legend.position = 'none' ) + xlim( c(-1, 28) ) + ylab( '' )

png( 'Ridge_Box_Comparison.png', height = 500, width = 900 )
ggarrange( shap_joy, shap_box, nrow = 1, ncol = 2 )
dev.off()

png( 'Bootstrapped_Shaps.png', height = 600, width = 1000 )
ggarrange( shap_joy, s_plot_trunc, nrow = 1, ncol = 2 )
dev.off()

png( 'Bootstrapped_Quantile_Shaps_Trunc.png', height = 600, width = 1000 )
print( s_plot_trunc )
dev.off()

png( 'Bootstrapped_Rank_Distributions.png', height = 600, width = 1000 )
print( shap_joy )
dev.off()

# How often ranked in the top 10
top10 <- apply( res$bootstrap_results$shap_ranks, 2, function(x) { mean(x < 10) } ) %>% 
  sort( decreasing = TRUE ) %>% print
plot( top10, type = 'h', xaxt = 'n', xlab = '', ylab = '' )
axis( 1, at = 1:length(top10), las = 2, labels = names(top10) )
abline( h = 0.3, lty = 2 )

# How often ranked above random noise
rn_pos <- res$bootstrap_results$shap_ranks[, 'random_num']
above_noise_fun <- function( x ) {
  (x <= rn_pos) %>% mean
}
apply( res$bootstrap_results$shap_ranks, 2, above_noise_fun ) %>% sort( decreasing = TRUE ) %>% print

### Bootstrapped ROC curves
# Function for interpolating ROC coordinates
inter_coords <- function( y, len = 500 ) {
  
  inter_len <- 1:len
  quants <- ( inter_len - min(inter_len) ) / ( max(inter_len) - min(inter_len) )
  
  sens <- coords( roc = y, x = quants, ret = 'sensitivity', inp = 'threshold' ) %>% unlist
  spec <- coords( roc = y, x = quants, ret = 'specificity', inp = 'threshold' ) %>% unlist
  
  data.frame( FPR = 1 - spec, Sensitivity = sens )
  
}

# Utility function for computing mean and quantiles of a column in a dataframe 
col_stats <- function( y, x ) {
  
  col_extract <- lapply( x, '[', y )
  col_set <- do.call( 'cbind', col_extract ) %>% as.matrix
  col_median <- apply( col_set, 1, quantile, 0.5 )
  col_q1 <- apply( col_set, 1, quantile, 0.1 )
  col_q2 <- apply( col_set, 1, quantile, 0.25 )
  col_q3 <- apply( col_set, 1, quantile, 0.75 )
  col_q4 <- apply( col_set, 1, quantile, 0.9 )
  
  data.frame( median = col_median, q1 = col_q1, q2 = col_q2, q3 = col_q3, q4 = col_q4 )
  
}

# Utility function to extract the correct columns from dataframes in a list and combine them
plot_cols <- function( x, y ) {
  cols <- lapply( x, '[', y )
  do.call( 'cbind', cols )
}

# Function for computing bootstrapped ROC confidence intervals
roc_ci <- function( x, col_poly ) {
  
  roc_comps <- c( 'FPR', 'Sensitivity' )
  sets <- lapply( roc_comps, col_stats, x = x )
  median_curve <- plot_cols( sets, 'median' )
  ul_curve1 <- cbind( sets[[1]]$q1, sets[[2]]$q4 )
  ll_curve1 <- cbind( rev(sets[[1]]$q4), rev(sets[[2]]$q1) )
  ul_curve2 <- cbind( sets[[1]]$q2, sets[[2]]$q3 )
  ll_curve2 <- cbind( rev(sets[[1]]$q3), rev(sets[[2]]$q2) )
  
  colnames( median_curve ) <- colnames( ul_curve1 ) <- colnames( ll_curve1 ) <- 
    colnames( ul_curve2 ) <- colnames( ll_curve2 ) <- roc_comps
  poly_curve1 <- rbind( ul_curve1, ll_curve1 ) %>% as.data.frame
  poly_curve2 <- rbind( ul_curve2, ll_curve2 ) %>% as.data.frame
  
  list( median_curve = median_curve, poly_curve1 = poly_curve1, poly_curve2 = poly_curve2 )
  
}

# Create interpolated ROC curves
train_rocs <- lapply( res$bootstrap_results$train_rocs, inter_coords )
test_rocs <- lapply( res$bootstrap_results$test_rocs, inter_coords )

train_ci <- roc_ci( train_rocs )
test_ci <- roc_ci( test_rocs )

boot_roc <- ggplot( train_ci$median_curve, aes(FPR, Sensitivity) ) + geom_line( col = mycol[1] ) +
  #geom_polygon( data = train_ci$poly_curve1, aes(FPR, Sensitivity), fill = mycol[1], alpha = 0.1 ) +
  geom_polygon( data = train_ci$poly_curve2, aes(FPR, Sensitivity), fill = mycol[1], alpha = 0.3 ) +
  geom_line( data = test_ci$median_curve, aes(FPR, Sensitivity), col = mycol[2] ) +
  #geom_polygon( data = test_ci$poly_curve1, aes(FPR, Sensitivity), fill = mycol[2], alpha = 0.1 ) +
  geom_polygon( data = test_ci$poly_curve2, aes(FPR, Sensitivity), fill = mycol[2], alpha = 0.3 ) +
  ggtitle( 'Bootstrapped ROC Curves' ) + geom_text(aes(0.5, 1, label = 'Train'), col = mycol[1], size = 7) + 
  geom_text(aes(0.8, 0.5, label = 'Test'), col = mycol[2], size = 7) + theme( text = element_text(size = 15) )

name_lab <- c( 'Train', 'Test' )
all_auc_long <- pivot_longer( as.data.frame(res$bootstrap_results$auc), cols = all_of(name_lab) )
all_auc_long$name <- factor( all_auc_long$name, levels = name_lab )
boot_auc <- ggplot( all_auc_long, aes(name, value) ) + geom_boxplot( fill = mycol ) + ylim(0, 1) +
  geom_hline( yintercept = 0.5, linetype = 'dashed', col = 'darkgrey' ) + ylab( 'AUC' ) + xlab( '' ) +
  ggtitle( 'Bootstrapped AUC' ) + theme( text = element_text(size = 15) )

png( 'Bootstrapped_ROC.png', height = 600, width = 1000 )
ggarrange( boot_roc, boot_auc, ncol = 2, nrow = 1 )
dev.off()

png( 'Bootstrapped_ROC_band.png', height = 600, width = 1000 )
print( boot_roc )
dev.off()

### Cross Validation Analysis
mkdir( 'CV_plots' )

param_grid <- do.call( 'rbind', res$cv_results$param_grid )
cv_ce <- apply( res$cv_results$cv_err_dat, 2, mean )
cv_auc <- apply( res$cv_results$cv_auc_dat, 2, mean )

# Utility function for finding the first instance of a value in a vector
first_val <- function( x, y, dat ) {
  which( dat[[y]] == x ) %>% min
}

# Utility function for creating all two-way concatenations in a vector
concat_fun <- function( x ) {
  
  res <- list()
  for ( i in 1:length(x) ) {
    y <- x[i]
    z <- x[-i]
    res_i <- lapply( z, function(w) { c(y, w) } )
    res <- c( res, res_i )
  }
  res
  
}

# Function to plot out how CV results change over a given parameter
CV_var_rel <- function( x ) {
  
  setwd( 'CV_plots' )
  
  primary <- x[1]
  secondary <- x[2]
  
  if( !is.na(secondary) ) {
    var_ord <- order(param_grid[[primary]], param_grid[[secondary]])
  }
  else {
    var_ord <- order(param_grid[[primary]])
  }
  pargrid_var <- param_grid[ var_ord, ]
  
  cv_ce_var <- data.frame( ind = 1:length(var_ord), val = cv_ce[ var_ord ] )
  cv_auc_var <- data.frame( ind = 1:length(var_ord), val = cv_auc[ var_ord ] )
  
  var_unique <- unique( param_grid[[primary]] ) %>% sort
  var_inds <- sapply( var_unique, first_val, y = primary, dat = pargrid_var )
  
  if ( !is.na(secondary) ) {
    sec_var <- pargrid_var[[secondary]]
    sec_var_set <- unique( sec_var ) %>% sort( decreasing = TRUE )
    col_scale_fun <- colorRampPalette( c('purple', 'yellow') )
    col_scale <- col_scale_fun( length( sec_var_set ) )
    use_cols <- col_scale[ match( sec_var, sec_var_set ) ]
    filename <- paste0( primary, '_', secondary )
    gtitle <- paste0( primary, ' and ', secondary )
  }
  else { 
    use_cols <- 'black' 
    filename <- gtitle <- primary
  }
  
  ce_min_ind <- which.min( cv_ce_var$val )
  
  png( paste0('CE_', filename, '.png') )
  ce_plot <- ggplot( cv_ce_var, aes(ind, val) ) + geom_line() + geom_point( col = use_cols ) +
    ggtitle( gtitle ) + geom_vline( xintercept = var_inds[-1], linetype = 'dashed', col = 'red' ) +
    geom_vline( xintercept = ce_min_ind, col = 'blue' ) + xlab( 'Index' ) + ylab( 'Cross-Entropy Error' )
  print( ce_plot )
  dev.off()
  
  png( paste0('AUC_', filename, '.png') )
  auc_plot <- ggplot( cv_auc_var, aes(ind, val) ) + geom_line() + geom_point( col = use_cols ) +
    ggtitle( gtitle ) + geom_vline( xintercept = var_inds[-1], linetype = 'dashed', col = 'red' ) +
    geom_vline( xintercept = ce_min_ind, col = 'blue' ) + xlab( 'Index' ) + ylab( 'AUC' )
  print( auc_plot )
  dev.off()
  
  setwd( '..' )
  
  NULL
  
}

# Run function for all tuning parameters
parnames <- colnames( param_grid )[-7] # ignore objective
parnames_concat <- concat_fun( parnames )
lapply( parnames, CV_var_rel )
lapply( parnames_concat, CV_var_rel )

### Bootstrapped Heat Maps to Analyse Consistency
mkdir( 'Heatmaps' )

all_vars <- unique(shap_long$variable) %>% as.vector
all_shaps <- res$bootstrap_results$all_shaps
all_shaps_val_list <- split( all_shaps$shap, all_shaps$variable )
all_shaps_val_list <- all_shaps_val_list[all_vars]
all_shaps_val <- do.call( 'cbind', all_shaps_val_list ) %>% as.data.frame
all_shaps_val$boot_ind <- rep( 1:B, each = 180 )

rank_shaps <- function( x ) {
  dat <- filter( all_shaps_val, boot_ind == x ) %>% select( -boot_ind )
  lapply( dat, rank ) %>% as.data.frame
}

all_shaps_rank <- do.call( 'rbind', lapply( 1:B, rank_shaps ) )
N <- ncol( all_shaps_rank )
all_shaps_rank$boot_ind <- all_shaps_val$boot_ind

for ( i in 1:N ) {
  
  cn <- all_vars[i]
  dat <- select( all_shaps_rank, all_of(c(i, N+1) ) )
  colnames( dat )[1] <- 'val'
  obs_ind <- order( X[, cn] )
  dat_list <- split( dat, rep( 1:B, each = 180 ) )
  reorder_list <- lapply( dat_list, function(x) { x$val[obs_ind] } )
  hm_dat <- data.frame( val = unlist(reorder_list), obs_ind = rep(1:180, B), boot_ind = rep(1:B, each = 180) )
  
  if ( length(unique(hm_dat$val)) == 1 ) { next }
  
  hm <- ggplot( hm_dat, aes(obs_ind, boot_ind, fill = val) ) + geom_tile() +
    xlab('Observation Rank') + ylab( 'Bootstrap Index' ) + 
    ggtitle( paste0( cn, ' Shapley Value Heatmap' ) ) + 
    theme(legend.position="bottom") + 
    scale_fill_gradient('Shapley Value', low="#FFCC33", high="#6600CC", 
                        breaks = range(hm_dat$obs_ind), labels = c("low", "high"))
  
  png( paste0('Heatmaps/shap_heatmap_', i, '.png') )
  print( hm )
  dev.off()
  
}

### Violin plots of top 10 variables
num_top <- 10
top_vars <- unique( shap_long$variable )[1:num_top] %>% as.vector
top_dat <- X[, top_vars] %>% as.data.frame %>% cbind( Exposure = as.factor(res$cv_results$y) )

group_cols <- mycol
use_cols <- group_cols[ match( top_dat$Exposure, c( 1, 0 ) ) ]

exposed_labs <- c('Non-Exposed', 'Exposed')
top_dat$Exposure_lab <- factor( exposed_labs[as.numeric(top_dat$Exposure)], levels = exposed_labs )

vio_plot_fun <- function( z ) {
  
  if ( z == "MonthsSinceStartingDialysis" ) {
    nm <- "Months on Dialysis" # This is so the title fits on the page
  }
  else {
    nm <- z
  }
  
  ggplot(top_dat, aes(x = Exposure_lab, y = get(z))) +
    geom_violin(alpha = 0.3, size = 0.1, fill = 'grey', col = 'grey') +
    geom_sina( aes(x = Exposure_lab, y = get(z), color = use_cols), method = "counts", maxwidth = 0.7, 
               alpha = 0.7 ) + theme(legend.position = "none", text = element_text(size = 18)) + 
    ggtitle( nm ) + ylab( '' ) + xlab( '' )
  
}

set.seed( 1 ) # To ensure the random jitter is reproducible
vio_plots <- lapply( top_vars, vio_plot_fun )

png( paste0('Top_Violins.png'), height = 1000, width = 1400 )
ggarrange( plotlist = vio_plots, ncol = 5, nrow = 2 )
dev.off()

### Create tables for clinical variables
### Numerical Table
bsi_dat <- read_csv( '../../../Data/bsi_dat.csv' )
tab_names <- c( 'Variable', 'No Exposure Mean \u00B1 1se', 'Exposure Mean \u00B1 1se', 'p-value' )
pm_paste_fun <- function( x, y ) {
  paste0( round(x, 1), ' \u00B1 ', round(y, 1) )
}

# Age
age_cohort <- matrix( '', nrow = 5, ncol = 4 ) %>% as.data.frame
names( age_cohort ) <- tab_names
age_cohort$Variable <- c( 'Age (years)', '20-45', '46-70', '>70', '' )
age_split <- list()
age_split[[1]] <- split( bsi_dat$AgeatConsent, bsi_dat$Exposure )
age_split[[2]] <- lapply( age_split[[1]], function(x) { x[x <= 45] } )
age_split[[3]] <- lapply( age_split[[1]], function(x) { x[(x > 45) & (x <= 70)] } )
age_split[[4]] <- lapply( age_split[[1]], function(x) { x[x > 70] } )

age_mean <- lapply( age_split, function(x) { sapply(x, 'mean') } )
age_se <- lapply( age_split, function(x) { sapply(x, 'sd')/sqrt(sapply(x, length)) } )
age_p <- sapply( age_split, function(x) { round( t.test(x[[1]], x[[2]])$p.value, 2 ) } )

age_cohort[1:4, 2:4] <- cbind( do.call( 'rbind', Map( pm_paste_fun, age_mean, age_se ) ), age_p )

# Months On Dialysis
mod_cohort <- matrix( '', nrow = 4, ncol = 4 ) %>% as.data.frame
names( mod_cohort ) <- tab_names
mod_cohort$Variable <- c( 'Months on Dialysis', '1-48', '49-96', '>96' )
mod_split <- list()
mod_split[[1]] <- split( bsi_dat$MonthsSinceStartingDialysis, bsi_dat$Exposure )
mod_split[[2]] <- lapply( mod_split[[1]], function(x) { x[x <= 48] } )
mod_split[[3]] <- lapply( mod_split[[1]], function(x) { x[(x > 48) & (x <= 96)] } )
mod_split[[4]] <- lapply( mod_split[[1]], function(x) { x[x > 96] } )

mod_mean <- lapply( mod_split, function(x) { sapply(x, 'mean') } )
mod_se <- lapply( mod_split, function(x) { sapply(x, 'sd')/sqrt(sapply(x, length)) } )
mod_p <- sapply( mod_split, function(x) { round( t.test(x[[1]], x[[2]])$p.value, 2 ) } )

mod_cohort[1:4, 2:4] <- cbind( do.call( 'rbind', Map( pm_paste_fun, mod_mean, mod_se ) ), mod_p )

numerical_table <- rbind( age_cohort, mod_cohort )
write.csv( numerical_table, 'Numerical_table.csv' )

### Categorical variables
cat_names <- c( 'Variable', 'No Exposure Amount', 'Exposure Amount', 'p-value' )

# Sex
sex_cohort <- matrix( '', nrow = 4, ncol = 4 ) %>% as.data.frame
names( sex_cohort ) <- cat_names
sex_cohort$Variable <- c( 'Sex', 'Male', 'Female', '' )
sex <- bsi_dat$Gender - 1
sex_cohort[2:3, 2:3] <- sex_props <- table( sex, bsi_dat$Exposure ) %>% as.matrix
sex_cohort[1, 4] <- round( prop.test( sex_props )$p.value, 2 )

# CKD Aetiology
ckd_cohort <- matrix( '', nrow = 9, ncol = 4 ) %>% as.data.frame
names( ckd_cohort ) <- cat_names
ckd_cohort$Variable <- c( 'Aetiology of CKD', 'Chronic Glomerulonephritis', 'Ischaemic Nephrology',
                          'Polycystic Kidney Disease', 'Diabetes', 'Congenital', 'Other',
                          'Unknown', '' )
ckd_dat <- split( select( bsi_dat, all_of( paste0('AetiologyofCKD.', 2:7) ) ), bsi_dat$Exposure )
ckd_rowsums <- sapply( ckd_dat, rowSums )
ckd_cohort[2, 2:3] <- sapply( ckd_rowsums, function(x) { sum( x == 0 ) } )

ind_sum <- function( x, i ) {
  sum( x[[i]], na.rm = TRUE )
}

for ( i in 3:8 ) {
  ckd_cohort[i, 2:3] <- sapply( ckd_dat, ind_sum, i = i-2 )
}

ckd_props <- lapply( ckd_cohort[2:8, 2:3], as.numeric ) %>% as.data.frame %>% as.matrix
ckd_cohort[1, 4] <- round( prop.test( ckd_props )$p.value, 2 )

# Previous Transplant
pt_cohort <- matrix( '', nrow = 5, ncol = 4 ) %>% as.data.frame
names( pt_cohort ) <- cat_names
pt_cohort$Variable <- c( 'Previous Transplant', 'Yes', 'No', 'Unknown', '' )
pt_cohort[2:4, 2:3] <- pt_props <- table( bsi_dat$TransplantEver, bsi_dat$Exposure, useNA = 'ifany' ) %>% as.matrix
pt_cohort[1, 4] <- round( prop.test( pt_props[-3, ] )$p.value, 2 )

# Immunosuppressive Therapy
ist_cohort <- matrix( '', nrow = 5, ncol = 4 ) %>% as.data.frame
names( ist_cohort ) <- cat_names
ist_cohort$Variable <- c( 'Immunosuppressive Therapy', 'Yes', 'No', 'Unknown', '' )
ist_cohort[2:4, 2:3] <- ist_props <- table( bsi_dat$CurrentImmunosuppressiveRx, bsi_dat$Exposure, useNA = 'ifany' )
ist_cohort[1, 4] <- round( prop.test( ist_props[-3, ] )$p.value, 2 )

# Vascular Access
va_cohort <- matrix( '', nrow = 5, ncol = 4 ) %>% as.data.frame
names( va_cohort ) <- cat_names
va_cohort$Variable <- c( 'Access Type', 'T-CVC', 'AVF', 'AVG', 'Unknown' )
va_dat <- split( select( bsi_dat, all_of( paste0('VascularAccess.', 2:3) ) ), bsi_dat$Exposure )
va_rowsums <- sapply( va_dat, rowSums )
va_cohort[2, 2:3] <- sapply( va_rowsums, function(x) { sum( x == 0, na.rm = TRUE ) } )
va_cohort[5, 2:3] <- sapply( va_rowsums, function(x) { sum( is.na(x) ) } )

for ( i in 3:4 ) {
  va_cohort[i, 2:3] <- sapply( va_dat, ind_sum, i = i-2 )
}

va_props <- lapply( va_cohort[2:4, 2:3], as.numeric ) %>% as.data.frame %>% as.matrix
va_cohort[1, 4] <- round( prop.test( va_props )$p.value, 2 )

categorical_table <- do.call( 'rbind', list(sex_cohort, ckd_cohort, pt_cohort, ist_cohort, va_cohort) )
write.csv( categorical_table, 'Categorical_table.csv' )

