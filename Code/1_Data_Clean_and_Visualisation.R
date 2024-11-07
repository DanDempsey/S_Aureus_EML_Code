##### Clean Clinical Data
##### Daniel Dempsey

### Load in required libraries and set colour scheme
library( readr ) # Reading and writing csv files
library( readxl ) # Reading excel files
library( dplyr ) # Useful data frame manipulation tools
library( caret ) # Creating dummy variables
library( ggforce ) # Sina plots
library( ggpubr ) # For placing multiple graphs on the same plot
library( naniar ) # Missing data plots
library( ggcorrplot ) # Correlation matrix plots
mycol <- c( '#56B4E9', '#E69F00' )

### Utility function: make directories if they don't already exist
mkdir <- function( x ) {
  if( !dir.exists( x ) ) {
    dir.create( x )
  }
}
mkdir( 'Output/Data_Visualisations' )

###### Data Cleaning
### Load in clinical data
clinic_dat_full <- read_xlsx( 'Data/HD-180 - FACS and excel matching -Feb 2023.xlsx', 
                              skip = 4, sheet = 'US & 3 antigens pooled' )

# Two columns list study numbers; check that they are consistent
all( clinic_dat_full$StudyNo...1 == clinic_dat_full$StudyNo...34 ) # All TRUE

# Delete one of the redundant StudyNo columns
studynum_inds <- which( substr(colnames(clinic_dat_full), 1, 7) == 'StudyNo' ) # Two columns 
clinic_dat <- select( clinic_dat_full, -studynum_inds[2] )
colnames( clinic_dat )[1] <- 'StudyNo'

# Check that the two Exposure variables are consistent
table( clinic_dat$HospitaldocumentedSAExposure, clinic_dat$Exposure ) # Consistent, remove the redundant column
clinic_dat$HospitaldocumentedSAExposure <- NULL

### Missing data check
png( 'Output/Data_Visualisations/Missing_Values_rawData.png' )
vis_miss( clinic_dat ) # Almost all of the raw data is observed
dev.off()
sort( sapply( clinic_dat_full, function(x) { sum(is.na(x)) } ), decreasing = TRUE )
# Myeloma has a very large number of missing values and will be removed.
# Other than that only 5 missing values: 2 in VascularAcess, and 1 each in MonthsSinceCurrentAccessEstablished, TransplantEver and CurrentImmunosuppressiveRx 

### Correct issues with the data
# CurrentImmunosuppressiveRx loaded as a character with inconsistent casing - convert it to binary
clinic_dat$CurrentImmunosuppressiveRx <- ifelse( tolower(clinic_dat$CurrentImmunosuppressiveRx) == 'y', 1, 0 )

# Diabetes loaded in as a character - should be numeric binary. Also one value is a 4. This is a misprint, should have been a 1. Correct these.
clinic_dat$Diabetes <- as.numeric( clinic_dat$Diabetes )
clinic_dat$Diabetes[which(clinic_dat$Diabetes == 4)] <- 1
clinic_dat$Diabetes <- as.factor( clinic_dat$Diabetes )

# Myeloma is binary with inconsistent codeding; y/n and 1/0. Convert it all to binary 
clinic_dat$Myeloma <- ifelse( clinic_dat$Myeloma == 'Y', '1', clinic_dat$Myeloma )
clinic_dat$Myeloma <- ifelse( clinic_dat$Myeloma %in% c('N', 'n'), '0', clinic_dat$Myeloma )
clinic_dat$Myeloma <- as.numeric( clinic_dat$Myeloma )
clinic_dat$Myeloma <- as.factor( clinic_dat$Myeloma )

# VascularAccess uses an ampersand to denote multiple groups. Code this in a more consistent way
clinic_dat$VascularAccess[which(clinic_dat$VascularAccess == '1&2')] <- '12'
clinic_dat$VascularAccess <- as.numeric( clinic_dat$VascularAccess )
clinic_dat$VascularAccess <- as.factor( clinic_dat$VascularAccess )

# AetiologyofCKD uses and ampersand to denote multiple groups. Code this in a more consistent way
clinic_dat$AetiologyofCKD[which(clinic_dat$AetiologyofCKD == '1&2')] <- '12'
clinic_dat$AetiologyofCKD[which(clinic_dat$AetiologyofCKD == '2&4')] <- '24'
clinic_dat$AetiologyofCKD[which(clinic_dat$AetiologyofCKD == '2&5')] <- '25'
clinic_dat$AetiologyofCKD[which(clinic_dat$AetiologyofCKD == '2 &7')] <- '27'
clinic_dat$AetiologyofCKD <- as.numeric( clinic_dat$AetiologyofCKD )
clinic_dat$AetiologyofCKD <- as.factor( clinic_dat$AetiologyofCKD )

# Convert remaining categorical data to a factor
clinic_dat$NasalScreen <- as.factor( clinic_dat$NasalScreen )

### Remove irrelevant variables
clinic_dat_trim <- clinic_dat
clinic_dat_trim$CollectionSite <- NULL # This is simply telling us where the patient was treated
clinic_dat_trim$Myeloma <- NULL # Too many missing values

### One-hot encoding of non-binary categorical variables
dv <- dummyVars( ~ VascularAccess + AetiologyofCKD + NasalScreen, clinic_dat_trim )
new_facs <- predict( dv, clinic_dat_trim ) %>% as.data.frame

# Assign codings to where there are multiple groups
va12_inds <- which( new_facs$VascularAccess.12 == 1 )
new_facs$VascularAccess.1[va12_inds] <- new_facs$VascularAccess.2[va12_inds] <- 1
new_facs$VascularAccess.12 <- NULL

ackd12_inds <- which( new_facs$AetiologyofCKD.12 == 1 )
new_facs$AetiologyofCKD.1[ackd12_inds] <- new_facs$AetiologyofCKD.2[ackd12_inds] <- 1
new_facs$AetiologyofCKD.12 <- NULL

ackd24_inds <- which( new_facs$AetiologyofCKD.24 == 1 )
new_facs$AetiologyofCKD.2[ackd24_inds] <- new_facs$AetiologyofCKD.4[ackd24_inds] <- 1
new_facs$AetiologyofCKD.24 <- NULL

ackd25_inds <- which( new_facs$AetiologyofCKD.25 == 1 )
new_facs$AetiologyofCKD.2[ackd25_inds] <- new_facs$AetiologyofCKD.5[ackd25_inds] <- 1
new_facs$AetiologyofCKD.25 <- NULL

ackd27_inds <- which( new_facs$AetiologyofCKD.27 == 1 )
new_facs$AetiologyofCKD.2[ackd27_inds] <- new_facs$AetiologyofCKD.7[ackd27_inds] <- 1
new_facs$AetiologyofCKD.27 <- NULL

### Combine one-hot-encoded variables and delete corresponding raw variables
clinic_clean <- cbind( clinic_dat_trim, new_facs )
clinic_clean$VascularAccess <- clinic_clean$AetiologyofCKD <-
  clinic_clean$NasalScreen <- clinic_clean$VascularAccess.1 <- 
  clinic_clean$AetiologyofCKD.1 <- clinic_clean$NasalScreen.0 <- NULL

### Remove non-corrected day 8 cells and other irrelevant data
drop_cols <- c( 'CD4Unstim', 'CD4HKPS80', 'CD4HKMRSA', 'CD4HKMSSA', 'MemCD4Unstim', 'MemCD4HKPS80', 
                'MemCD4HKMRSA', 'MemCD4HKMSSA', 'DivMem CD4Unstim', 'DivMemCD4HKPS80', 'DivMemCD4HKMRSA', 
                'DivMemCD4HKMSSA', 'IL10Unstim', 'IL10HKPS80', 'IL10HKMRSA', 'IL10HKMSSA', 'Th17Unstim', 
                'Th17HKPS80', 'Th17HKMRSA', 'Th17HKMSSA', 'TransTh17Unstim', 'TransTh17HKPS80', 'TransTh17HKMRSA', 
                'TransTh17HKMSSA', 'Th1Unstim', 'Th1HKPS80', 'Th1HKMRSA', 'Th1HKMSSA', 'exTh17Unstim', 'exTh17HKPS80', 
                'exTh17HKMRSA', 'exTh17HKMSSA', 'CD8', 'EffectorTcells', 'CD3', 'CD3_A',
                'MonthsSinceCurrentAccessEstablished', 'Diabetes', 'NasalScreen.1', 'NasalScreen.2' )

dat_reduced <- select( clinic_clean, -all_of(drop_cols) )

### Append a randomly generated column for model validation purposes
set.seed( 123 )
dat_reduced$random_num <- abs(10*rnorm(nrow(dat_reduced)))

### Correlation matrix of current data
cor_mat_raw <- cor( select( dat_reduced, -c(Exposure, StudyNo) ),  use = 'complete.obs' )
cor_mat_reversed <- cor_mat_raw[, nrow(cor_mat_raw):1 ]

pro_name_change <- gsub( 'Div', 'Pro', colnames(cor_mat_raw) )
colnames( cor_mat_reversed ) <- rev( pro_name_change )
rownames( cor_mat_reversed ) <- pro_name_change

png( 'Output/Data_Visualisations/Correlation_Matrix_Raw.png', height = 900, width = 900 )
cor_mat_raw_plot <- ggcorrplot( cor_mat_reversed, colors = c("#FFCC33", "#FFFFFF", "#6600CC") ) + 
  scale_y_discrete() +
  theme( plot.margin = unit(c(0, 0, 0, 0), 'centimeters') ) + 
  theme(plot.title = element_text(vjust = 3, size = 30))
print( cor_mat_raw_plot )
dev.off()

### Dimension reduction of different Day 8 strains using PCA
PCA_vars <- c( 'MemCD4corHKPS80', 'MemCD4corHKMRSA', 'MemCD4corHKMSSA',
               'DivMemCD4corHKPS80', 'DivMemCD4corHKMRSA', 'DivMemCD4corHKMSSA',
               'IL10corHKPS80', 'IL10corHKMRSA', 'IL10corHKMSSA', 
               'Th17corHKPS80', 'Th17corHKMRSA', 'Th17corHKMSSA',
               'TransTh17corHKPS80', 'TransTh17corHKMRSA', 'TransTh17corHKMSSA',
               'Th1corHKPS80', 'Th1corHKMRSA', 'Th1corHKMSSA',
               'exTh17corHKPS80', 'exTh17corHKMRSA', 'exTh17corHKMSSA' )

PCA_var_split <- split( PCA_vars, c(rep(1, 6), rep(2:6, each = 3)) )
names( PCA_var_split ) <- c( 'memCD4', 'IL10', 'Th17', 'transTh17', 'Th1', 'exTh17' )

PCA_fun <- function( x ) {
  prcomp( dat_reduced[, x], scale. = TRUE )
} 

PCA_list_raw <- lapply( PCA_var_split, PCA_fun )
PCA_list <- lapply( PCA_list_raw, function(x) { 
  if ( all(x$rotation[, 1] < 0) ) {
    x$x[, 1] <- -x$x[, 1]
  }
  x } ) # This is so the first principal component is easier to interpret

# Highlight principal components that explain > 80% of the variance
plot_PCA <- function( x, title, ... ) {
  
  var_ex <- x$sdev^2 / sum( x$sdev^2 ) 
  var_ex_cumsum <- cumsum( var_ex )
  last_pca <- which( var_ex_cumsum > 0.8 )[1]
  col_ind <- c( rep( 1, last_pca ), rep( 2, ncol(x$x) - last_pca ) )
  
  var_ex_dat <- data.frame( Var1 = paste0('PC', 1:length(col_ind)), Freq = var_ex )
  var_ex_cumsum_dat <- data.frame( Var1 = paste0('PC', 1:length(col_ind)), Freq = var_ex_cumsum )
  
  ggplot( var_ex_dat, aes( x = Var1, y = Freq ) ) +
    geom_bar(stat="identity", fill=c('cyan', 'grey')[col_ind], col = 'black') + ylim(0, 1) +
    ylab('Proportion of Variance Explained') + xlab('') + ggtitle(paste0(title, ' Principal Components')) +
    theme( text = element_text(size = 15) )
  
}

PCA_plots <- Map( plot_PCA, x = PCA_list, title = names( PCA_list ) )

png( 'Output/Data_Visualisations/PCA_Variance_Explained.png', height = 700, width = 1000 )
ggarrange( plotlist = PCA_plots, nrow = 2, ncol = 3 )
dev.off()

add_cols <- function( var_name, dat = dat_reduced ) {
  
  x <- PCA_list[[var_name]]
  var_ex <- x$sdev^2 / sum( x$sdev^2 ) 
  var_ex_cumsum <- cumsum( var_ex )
  var_ind <- 1:which( var_ex_cumsum > 0.8 )[1]
  
  var_pca <- x$x[, var_ind, drop = FALSE]
  
  colnames( var_pca ) <- paste0( paste0(var_name, '_PC'), var_ind )
  cbind( dat, var_pca )
  
}

for ( i in 1:length(PCA_list) ) {
  dat_reduced <- add_cols( names(PCA_list)[i] )
}

final_dat <- select( dat_reduced, -all_of(PCA_vars) )
dim( final_dat ) # 180 x 39 (including StudyNo and response variable)

### Final adjustment: change column name to remove hyphen as this causes glitches in imputation code
ncm_ind <- which( colnames(final_dat) == 'Non-Classical Monocytes' )
colnames( final_dat )[ncm_ind] <- 'NonClassicalMonocytes'

### Write to file
write_csv( final_dat, 'Data/bsi_dat.csv' )

###### Data Visualisation
setwd('Output/Data_Visualisations/')

### Plot response data
expose_dat <- as.data.frame( table( final_dat$Exposure ) )

exposure_plot <- ggplot( expose_dat, aes( x = Var1, y = Freq ) ) +
  geom_bar(stat="identity", fill=mycol)+
  ylab('') + xlab('Exposure Status') + ggtitle('Proportion of Exposure Status')

png( 'Exposure_Barplot.png' )
exposure_plot
dev.off()

### Plot categorical covariates
tab_fun <- function( x ) {
  tab <- as.data.frame( table(clinic_dat[[x]], clinic_dat$Exposure) )
  colnames(tab)[1:2] <- c( x, 'Exposure' )
  tab
}

cat_list <- list()

cat_list[[1]] <- tab_fun( 'Gender' )
g_nms <- c('Male', 'Female')
cat_list[[1]]$Gender <- rep( g_nms, 2 )
cat_list[[1]]$Gender <- factor( cat_list[[1]]$Gender, levels = g_nms )
cat_list[[2]] <- tab_fun( 'CurrentImmunosuppressiveRx' )
cat_list[[3]] <- tab_fun( 'Diabetes' )
cat_list[[4]] <- tab_fun( 'TransplantEver' )
cat_list[[5]] <- tab_fun( 'Myeloma' )
cat_list[[6]] <- tab_fun( 'VascularAccess' )
v_nms <- c('T-CVC', 'AVF', 'AVG', 'T-CVC & AVF')
cat_list[[6]]$VascularAccess <- rep( v_nms, 2 )
cat_list[[6]]$VascularAccess <- factor( cat_list[[6]]$VascularAccess, levels = v_nms )
cat_list[[7]] <- tab_fun( 'AetiologyofCKD' )
a_nms <- c('Chronic Glomerulonephritis', 'Ischaemic Nephrology', 
           'Polycystic Kidney Disease', 'Diabetes', 'Congenital',
           'Other', 'Unknown', 'Chronic Glomerulonephritis & Ischaemic Nephrology',
           'Ischaemic Nephrology & Diabetes', 'Ischaemic Nephrology & Congenital',
           'Ischaemic Nephrology & Unknown')
cat_list[[7]]$AetiologyofCKD <- rep( a_nms, 2 )
cat_list[[7]]$AetiologyofCKD <- factor( cat_list[[7]]$AetiologyofCKD, levels = a_nms )
cat_list[[8]] <- tab_fun( 'NasalScreen' )
n_nms <- c('Not Done', 'Negative', 'Positive')
cat_list[[8]]$NasalScreen <- rep( n_nms, 2 )
cat_list[[8]]$NasalScreen <- factor( cat_list[[8]]$NasalScreen, n_nms )

mkdir( 'Categorical_Barplots' )
for( i in 1:length(cat_list) ) {
  
  nm <- colnames( cat_list[[i]] )[1]
  png( paste0('Categorical_Barplots/', nm, '_Barplot.png') )
  cat_plot <- ggplot(data=cat_list[[i]], aes(x=cat_list[[i]][, 1], y=Freq, fill=Exposure)) +
    geom_bar(stat="identity", position=position_dodge()) + ylab('') + xlab( nm ) +
    scale_fill_manual(values=mycol) + ggtitle( paste0('Distribution of ', nm) ) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  print( cat_plot )
  dev.off()
  
}

### Plot numerical covariates
cat_vars <- c( 'StudyNo', 'Gender', 'CurrentImmunosuppressiveRx', 'TransplantEver', 
               'TransplantEver', 'VascularAccess.2', 'VascularAccess.3', 
               'AetiologyofCKD.2', 'AetiologyofCKD.3', 'AetiologyofCKD.4', 
               'AetiologyofCKD.5', 'AetiologyofCKD.6', 'AetiologyofCKD.7' )

num_dat <- select( final_dat, -all_of(cat_vars) ) 
num_dat$Exposure <- factor( num_dat$Exposure )

group_cols <- mycol
use_cols <- group_cols[ match( num_dat$Exposure, c( 1, 0 ) ) ]

mkdir( 'Numerical_Violins' )
for ( i in 1:ncol( num_dat )  ) {
  
  nm <- colnames( num_dat )[i]
  if (nm == 'Exposure') { next }
  png( paste0('Numerical_Violins/', nm, '_Violin.png') )
  vplot <- ggplot(num_dat, aes(x = Exposure, y = num_dat[[i]])) +
    geom_violin(alpha = 0.3, size = 0.1, fill = 'grey', col = 'grey') +
    geom_sina( aes(x = Exposure, y = num_dat[[i]], color = use_cols), method = "counts", maxwidth = 0.7, 
               alpha = 0.7 ) +
    theme(legend.position = "none") +
    ggtitle( paste0('Distribution of ', nm, ' by Exposure Status') ) + ylab( nm ) + xlab( 'Exposure Status' )
  print( vplot )
  dev.off()
  
}

### Plot of correlation matrix
cor_mat_all <- cor( select( final_dat, -c(Exposure, StudyNo) ),  use = 'complete.obs' )
cor_mat_num <- cor( select( num_dat, -Exposure ), use = 'complete.obs' )

png( 'Correlation_Matrix_All.png', height = 750, width = 750 )
cor_mat_all_plot <- ggcorrplot( cor_mat_all, colors = c("#FFCC33", "#FFFFFF", "#6600CC") ) + 
  theme( plot.margin = unit(c(0, 0, 0, 0), 'centimeters') ) + 
  ggtitle("Processed Data") + theme(plot.title = element_text(vjust = 3, size = 30))
print( cor_mat_all_plot )
dev.off()

png( 'Correlation_Matrix_Numeric.png', height = 750, width = 750 )
cor_mat_num_plot <- ggcorrplot( cor_mat_num, colors = c("#FFCC33", "#FFFFFF", "#6600CC") ) + 
  theme( plot.margin = unit(c(0, 0, 0, 0), 'centimeters') ) + 
  ggtitle("Processed Data") + theme(plot.title = element_text(vjust = 3, size = 30))
print( cor_mat_num_plot )
dev.off()

### Missing values
png( 'Missing_Values_processedData.png' )
vis_miss( final_dat )
dev.off()

apply( final_dat, 1, function(x) { any(is.na(x))} ) %>% sum # 4 patients have missing values

### exTh17
dat_reduced$Exposure <- factor( dat_reduced$Exposure )
ggplot(dat_reduced, aes(x = factor(Exposure), y = exTh17corHKMRSA)) +
  geom_violin(alpha = 0.3, size = 0.1, fill = 'grey', col = 'grey') +
  geom_sina( aes(x = Exposure, y = exTh17corHKMRSA, color = use_cols), method = "counts", maxwidth = 0.7, 
             alpha = 0.7 ) +
  theme(legend.position = "none") +
  ggtitle( paste0('Distribution of exTh17corHKMRSA by Exposure Status') ) + ylab( nm ) + xlab( 'Exposure Status' )

ggplot(dat_reduced, aes(x = factor(Exposure), y = exTh17corHKMSSA)) +
  geom_violin(alpha = 0.3, size = 0.1, fill = 'grey', col = 'grey') +
  geom_sina( aes(x = Exposure, y = exTh17corHKMSSA, color = use_cols), method = "counts", maxwidth = 0.7, 
             alpha = 0.7 ) +
  theme(legend.position = "none") +
  ggtitle( paste0('Distribution of exTh17corHKMSSA by Exposure Status') ) + ylab( nm ) + xlab( 'Exposure Status' )

ggplot(dat_reduced, aes(x = factor(Exposure), y = exTh17corHKPS80)) +
  geom_violin(alpha = 0.3, size = 0.1, fill = 'grey', col = 'grey') +
  geom_sina( aes(x = Exposure, y = exTh17corHKPS80, color = use_cols), method = "counts", maxwidth = 0.7, 
             alpha = 0.7 ) +
  theme(legend.position = "none") +
  ggtitle( paste0('Distribution of exTh17corHKPS80 by Exposure Status') ) + ylab( nm ) + xlab( 'Exposure Status' )

