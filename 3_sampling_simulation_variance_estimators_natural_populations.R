#######################################################
# 
# Spatially aggregated populations sampling simulation
# 
# 3_sampling_simulation_variance_estimators_natural_populations.R
# 
# jan.perret@gmail.com
#######################################################

### 
### simulate sampling surveys on the three natural populations to assess the performance of variance estimators
### 

### clean workspace
rm(list = ls())

### load packages
source("./R/imports.R")
library(parallel)
library(doParallel)
library(foreach)
library(pscl)

### load functions
files.sources <- list.files("./R", full.names = TRUE)
invisible(sapply(files.sources, source))

### load the data from the 3 natural populations
matrix_SANMIN <- read.csv2(file = "./data/tidy/map_sanguisorba_minor.csv", header = FALSE)
matrix_LIMGIR <- read.csv2(file = "./data/tidy/map_limonium_girardianum.csv", header = FALSE)
matrix_BELSYL <- read.csv2(file = "./data/tidy/map_bellis_sylvestris.csv", header = FALSE)

# group cells to switch from 10*10 cm to 20*20 cm cells
group_cells_matrix <- function(old_matrix) {
  
  matrix_new <- data.frame(matrix(NA, nrow = 100, ncol = 100))
  row_vect <- seq(from = 1, to = 200, by = 2)
  col_vect <- seq(from = 1, to = 200, by = 2)
  
  for (i in row_vect) {
    for (j in col_vect) {
      
      matrix_new[ceiling(i/2), ceiling(j/2)] <- old_matrix[i, j] + old_matrix[i+1, j] + old_matrix[i, j+1] + old_matrix[i+1, j+1]
      
    }
  }
  
  return(matrix_new)
  
}

matrix_SANMIN_new <- group_cells_matrix(old_matrix = matrix_SANMIN)
matrix_LIMGIR_new <- group_cells_matrix(old_matrix = matrix_LIMGIR)
matrix_BELSYL_new <- group_cells_matrix(old_matrix = matrix_BELSYL)

list_natural_pop <- list(matrix_SANMIN_new, matrix_LIMGIR_new, matrix_BELSYL_new)

# table with the names of the 3 populations
coef_tab_total <- data.frame(pattern = c("SANMIN", "LIMGIR", "BELSYL"),
                             intensity = c(1.4696, 1.6568, 0.0789),
                             pop_size = c(14696, 16568, 789),
                             dispersion_index = c(7.631599, 8.421541, 2.868155))


### general settings for the simulations
nsim = 1000
my_window <- owin(c(0, 100), c(0, 100))
window_size <- area(my_window)
sample_size_vect <- c(9, 16, 25, 36, 49, 64, 81, 100, 121, 144, 169, 196, 225,
                      256, 289, 324, 361, 400, 441, 484, 529, 576, 625)


# set the folder to save the outputs
storage_folder <- "RUN_2022-07-26_variance_estimators_natural_populations_with_ZIP"


# store and print start time
start.time_total <- Sys.time()
start.time_total

# set the number of cores to use for parallel computation
registerDoParallel(2)


### loop over the populations
# for (num_pattern in 1:length(list_natural_pop)) {
foreach (num_pattern = 1:length(list_natural_pop), .packages = c("spatstat", "sp", "spdep")) %dopar% {
    
  
  # get the current population
  mymatrix <- list_natural_pop[[num_pattern]]
  
  # get the population names
  coef_tab <- coef_tab_total[num_pattern, ]
  
  # store current intensity
  current_intensity <- coef_tab$intensity[1]
  
  
  ### loop over sample sizes
  for (sample_size in sample_size_vect) {

    
    # create the tables to store the results
    result_random <- data.frame(matrix(NA, nrow = nsim, ncol = 1))
    result_syst <- data.frame(matrix(NA, nrow = nsim, ncol = 1))
    result_BAS <- data.frame(matrix(NA, nrow = nsim, ncol = 1))
    prop_zero_random <- data.frame(matrix(NA, nrow = nsim, ncol = 1))
    prop_zero_syst <- data.frame(matrix(NA, nrow = nsim, ncol = 1))
    prop_zero_BAS <- data.frame(matrix(NA, nrow = nsim, ncol = 1))
    estim_vtextbook_random <- data.frame(matrix(NA, nrow = nsim, ncol = 1))
    estim_vtextbook_syst <- data.frame(matrix(NA, nrow = nsim, ncol = 1))
    estim_vtextbook_BAS <- data.frame(matrix(NA, nrow = nsim, ncol = 1))
    estim_v8_syst <- data.frame(matrix(NA, nrow = nsim, ncol = 1))
    estim_vw_syst <- data.frame(matrix(NA, nrow = nsim, ncol = 1))
    estim_vstr2_syst <- data.frame(matrix(NA, nrow = nsim, ncol = 1))
    estim_vnbh_BAS <- data.frame(matrix(NA, nrow = nsim, ncol = 1))
    estim_vZIP_syst <- data.frame(matrix(NA, nrow = nsim, ncol = 1))
    estim_vZIP_BAS <- data.frame(matrix(NA, nrow = nsim, ncol = 1))
    
    # rename the column
    colnames(result_random) <- coef_tab$pattern[1]
    colnames(result_syst) <- coef_tab$pattern[1]
    colnames(result_BAS) <- coef_tab$pattern[1]
    colnames(prop_zero_random) <- coef_tab$pattern[1]
    colnames(prop_zero_syst) <- coef_tab$pattern[1]
    colnames(prop_zero_BAS) <- coef_tab$pattern[1]
    colnames(estim_vtextbook_random) <- coef_tab$pattern[1]
    colnames(estim_vtextbook_syst) <- coef_tab$pattern[1]
    colnames(estim_vtextbook_BAS) <- coef_tab$pattern[1]
    colnames(estim_v8_syst) <- coef_tab$pattern[1]
    colnames(estim_vw_syst) <- coef_tab$pattern[1]
    colnames(estim_vstr2_syst) <- coef_tab$pattern[1]
    colnames(estim_vnbh_BAS) <- coef_tab$pattern[1]
    colnames(estim_vZIP_syst) <- coef_tab$pattern[1]
    colnames(estim_vZIP_BAS) <- coef_tab$pattern[1]
    
    
    # make nsim iterations
    for (i in 1:nsim) {
      
      # draw the samples
      mysample_random <- sample_random(pop = mymatrix, nx = my_window$xrange[2], ny = my_window$yrange[2], sample_size = sample_size)
      mysample_syst <- sample_systematic_fixed_size(pop = mymatrix, nx = my_window$xrange[2], ny = my_window$yrange[2], sample_size = sample_size)
      mysample_BAS <- sample_BAS(pop = mymatrix, nx = my_window$xrange[2], ny = my_window$yrange[2], sample_size = sample_size, sequence = "halton")
      
      # store the sample mean
      result_random[i, 1] <- mean(mysample_random)
      result_syst[i, 1] <- mean(mysample_syst)
      result_BAS[i, 1] <- mean(mysample_BAS)
      
      # store the proportion of zeros in the sample
      prop_zero_random[i, 1] <- sum(mysample_random == 0) / sample_size
      prop_zero_syst[i, 1] <- sum(mysample_syst == 0) / sample_size
      prop_zero_BAS[i, 1] <- sum(mysample_BAS == 0) / sample_size
      
      # calculate variance estimators
      # 1. v_textbook estimator for the 3 sampling methods
      estim_vtextbook_random[i, 1] <- var_estim_vtextbook(y = mysample_random, n = sample_size, N = window_size)
      estim_vtextbook_syst[i, 1] <- var_estim_vtextbook(y = mysample_syst, n = sample_size, N = window_size)
      estim_vtextbook_BAS[i, 1] <- var_estim_vtextbook(y = mysample_BAS, n = sample_size, N = window_size)
      
      # 2. other estimators for syst and BAS
      estim_v8_syst[i, 1] <- var_estim_v8(y = mysample_syst, n = sample_size, N = window_size)
      estim_vw_syst[i, 1] <- var_estim_vw(y = mysample_syst, n = sample_size, N = window_size)
      estim_vstr2_syst[i, 1] <- var_estim_vstr2(y = mysample_syst, n = sample_size, N = window_size, nx = my_window$xrange[2], ny = my_window$yrange[2])
      estim_vnbh_BAS[i, 1] <- var_estim_vnbh(y = mysample_BAS, n = sample_size, N = window_size)
      estim_vZIP_syst[i, 1] <- var_estim_ZIP(y = mysample_syst, n = sample_size, N = window_size)
      estim_vZIP_BAS[i, 1] <- var_estim_ZIP(y = mysample_BAS, n = sample_size, N = window_size)
      
      
    } # i
    
    
    # calculate the synthetic indicators and fill the result table
    result_tab <- cbind(coef_tab,
                        num_pattern = paste0("pattern_", num_pattern),
                        nsim = nsim,
                        sample_size = sample_size)
    
    result_tab <- cbind(result_tab,
                        mean_random = mean(result_random[, 1]),
                        median_random = median(result_random[, 1]),
                        var_random = var(result_random[, 1]),
                        mean_syst = mean(result_syst[, 1]),
                        median_syst = median(result_syst[, 1]),
                        var_syst = var(result_syst[, 1]),
                        mean_BAS = mean(result_BAS[, 1]),
                        median_BAS = median(result_BAS[, 1]),
                        var_BAS = var(result_BAS[, 1]),
                        prop_zero_random = mean(prop_zero_random[, 1]),
                        prop_zero_syst = mean(prop_zero_syst[, 1]),
                        prop_zero_BAS = mean(prop_zero_BAS[, 1]),
                        mean_estim_vtextbook_random = mean(estim_vtextbook_random[, 1]),
                        mean_estim_vtextbook_syst = mean(estim_vtextbook_syst[, 1]),
                        mean_estim_vtextbook_BAS = mean(estim_vtextbook_BAS[, 1]),
                        mean_estim_v8_syst = mean(estim_v8_syst[, 1]),
                        mean_estim_vw_syst = mean(estim_vw_syst[, 1]),
                        mean_estim_vstr2_syst = mean(estim_vstr2_syst[, 1]),
                        mean_estim_vnbh_BAS = mean(estim_vnbh_BAS[, 1]),
                        mean_estim_vZIP_syst = mean(estim_vZIP_syst[, 1]),
                        mean_estim_vZIP_BAS = mean(estim_vZIP_BAS[, 1]))
    
    result_tab <- cbind(result_tab,
                        biais_random = result_tab$mean_random - result_tab$intensity,
                        biais_syst = result_tab$mean_syst - result_tab$intensity,
                        biais_BAS = result_tab$mean_BAS - result_tab$intensity,
                        MSE_random = sapply(((result_random - current_intensity)^2), mean),
                        MSE_syst = sapply(((result_syst - current_intensity)^2), mean),
                        MSE_BAS = sapply(((result_BAS - current_intensity)^2), mean),
                        RMSE_random = sqrt(sapply(((result_random - current_intensity)^2), mean)),
                        RMSE_syst = sqrt(sapply(((result_syst - current_intensity)^2), mean)),
                        RMSE_BAS = sqrt(sapply(((result_BAS - current_intensity)^2), mean)),
                        var_ratio_syst = result_tab$var_syst / result_tab$var_random,
                        var_ratio_BAS = result_tab$var_BAS / result_tab$var_random)
    
    # compute the coverage obtained with each variance estimator
    # 1. textbook estimator
    my_Ybars <- c(rep(current_intensity, times = nsim)) # vector with nsim times the intensity of the current population
    coverage_vtextbook_random <- compute_coverage(Ybars = my_Ybars, ybars = result_random[, 1], estim_vars = estim_vtextbook_random[, 1], n = sample_size, N = window_size)
    coverage_vtextbook_syst <- compute_coverage(Ybars = my_Ybars, ybars = result_syst[, 1], estim_vars = estim_vtextbook_syst[, 1], n = sample_size, N = window_size)
    coverage_vtextbook_BAS <- compute_coverage(Ybars = my_Ybars, ybars = result_BAS[, 1], estim_vars = estim_vtextbook_BAS[, 1], n = sample_size, N = window_size)
    
    # 2. other estimators for syst and BAS
    coverage_v8_syst <- compute_coverage(Ybars = my_Ybars, ybars = result_syst[, 1], estim_vars = estim_v8_syst[, 1], n = sample_size, N = window_size)
    coverage_vw_syst <- compute_coverage(Ybars = my_Ybars, ybars = result_syst[, 1], estim_vars = estim_vw_syst[, 1], n = sample_size, N = window_size)
    coverage_vstr2_syst <- compute_coverage(Ybars = my_Ybars, ybars = result_syst[, 1], estim_vars = estim_vstr2_syst[, 1], n = sample_size, N = window_size)
    coverage_vnbh_BAS <- compute_coverage(Ybars = my_Ybars, ybars = result_BAS[, 1], estim_vars = estim_vnbh_BAS[, 1], n = sample_size, N = window_size)
    coverage_vZIP_syst <- compute_coverage(Ybars = my_Ybars, ybars = result_syst[, 1], estim_vars = estim_vZIP_syst[, 1], n = sample_size, N = window_size)
    coverage_vZIP_BAS <- compute_coverage(Ybars = my_Ybars, ybars = result_BAS[, 1], estim_vars = estim_vZIP_BAS[, 1], n = sample_size, N = window_size)
    
    # join the coverage to the result table
    result_tab <- cbind(result_tab,
                        coverage_vtextbook_random = coverage_vtextbook_random,
                        coverage_vtextbook_syst = coverage_vtextbook_syst,
                        coverage_vtextbook_BAS = coverage_vtextbook_BAS,
                        coverage_v8_syst = coverage_v8_syst,
                        coverage_vw_syst = coverage_vw_syst,
                        coverage_vstr2_syst = coverage_vstr2_syst,
                        coverage_vnbh_BAS = coverage_vnbh_BAS,
                        coverage_vZIP_syst = coverage_vZIP_syst,
                        coverage_vZIP_BAS = coverage_vZIP_BAS)
    
    ### save the results
    write.csv2(result_tab, file = paste0("./output/text/", storage_folder, "/simul_var_estimator_natural_populations_", coef_tab$pattern[1],
                                         "_sample-size_", sample_size, ".csv"), row.names = FALSE)
    
  } # end of loop over sample_size
  
} # end of loop over num_pattern

# close the cluster
stopImplicitCluster()

# print total execution time
end.time_total <- Sys.time()
time.taken_total <- difftime(end.time_total, start.time_total, units = "mins")
end.time_total
time.taken_total



#### Load and merge the simulation result tables ####

# load and merge the simulation output csv files

folder_path <- paste0("./output/text/", storage_folder, "/")
filenames <- paste0(folder_path, list.files(path = folder_path, pattern = "*.csv"))
files <- lapply(filenames, read.csv2)
merged_data <- data.table::rbindlist(files, use.names = FALSE)

# sort the table
merged_data <- merged_data[with(merged_data, order(pattern, sample_size)), ]

# save the dataframe
write.csv2(merged_data, file = "./data/tidy/RUN_2022-07-26_variance_estimators_natural_populations_with_ZIP_results_merged.csv", row.names = FALSE)

