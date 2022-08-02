#######################################################
# 
# Spatially aggregated populations sampling simulation
# 
# 3_sampling_simulation_variance_estimators_virtual_populations.R
# 
# jan.perret@gmail.com
#######################################################

### 
### simulate sampling surveys on the virtual populations to assess the performance of variance estimators, including a new estimator based on random sampling of a Zero-inflated Poisson (ZIP) distribution
### 

### clean workspace
rm(list = ls())

### load packages
source("./R/imports.R")
library(parallel)
library(doParallel)
library(foreach)

### load functions
files.sources <- list.files("./R", full.names = TRUE)
invisible(sapply(files.sources, source))



#### run sampling simulations ####

### general settings for the simulations
vect_num_pattern = c(1:20)
nsim = 600
my_window <- owin(c(0, 100), c(0, 100))
window_size <- area(my_window)
sample_size_vect = c(9, 25, 49, 100, 196)

# load parameter values for the 3 point processes
coef_tab_Matern <- read.csv2("./data/tidy/Matern_parameters_coef_tab_total_2021-06-13.csv") # modality 1
coef_tab_Matern <- cbind(coef_tab_Matern[, 1:2], process_type = "Matern_cluster", coef_tab_Matern[, 3:5], my_r = NA, expect_disp_ind = coef_tab_Matern[, 6])
coef_tab_Poisson <- read.csv2("./data/tidy/Poisson_parameters_coef_tab_total.csv")
coef_tab_Poisson <- cbind(coef_tab_Poisson[, 1:2], process_type = "poisson", coef_tab_Poisson[, 3:5], my_r = NA, expect_disp_ind = coef_tab_Poisson[, 6])
coef_tab_SSI <- read.csv2("./data/tidy/SSI_parameters_coef_tab_total_2021-10-06.csv")
coef_tab_SSI <- cbind(coef_tab_SSI[, 1:2], process_type = "SSI", my_kappa = NA, my_mu = NA, my_scale = NA, coef_tab_SSI[, 3:4])
coef_tab_total <- rbind(coef_tab_SSI, coef_tab_Poisson, coef_tab_Matern)

# add a column with the theoretical intensity value written as character
my_intensity <- paste0("intensity = ", coef_tab_total$theoretical_intensity)
coef_tab_total <- cbind(my_intensity = my_intensity, coef_tab_total)

# set the folder to save the outputs
storage_folder <- "RUN_2022-07-27_variance_estimator_ZIP_results_int1-10-20_nsim600_20_patterns"

# create a vector with the names of all files in the directory containing the previously generated count matrices
base_dir <- "./output/text/RUN_2022-06-08-simul_aggregation_intensity_patterns_modality_1/"
vect_all_files <- paste0(base_dir, list.files(base_dir, recursive = TRUE))
vect_all_files <- vect_all_files[seq(from = 1, to = length(vect_all_files), by = 2)] # keep every 2nd file to save computation time



# store and print start time
start.time_total <- Sys.time()
start.time_total

# set the number of cores to use for parallel computation
registerDoParallel(2)


### loop over the files containing the count matrices in vect_all_files
# for (my_file in vect_all_files) {
foreach (my_file = vect_all_files, .packages = c("spatstat", "sp", "spdep", "pscl")) %dopar% {
  
  # get the intensity and aggregation level of the current pattern
  vect_file_name <- as.vector(strsplit(my_file, split = "[-_.]+"))[[1]]
  my_intensity <- as.numeric(vect_file_name[10])
  
  if (vect_file_name[13] %in% c("04", "06", "08")){
    if (vect_file_name[13] == "04") {my_disp_ind <- 0.4}
    if (vect_file_name[13] == "06") {my_disp_ind <- 0.6}
    if (vect_file_name[13] == "08") {my_disp_ind <- 0.8}
  } else {
    my_disp_ind <- as.numeric(vect_file_name[13]) 
  }
  
  # get point process parameter values for the current intensity and aggregation levels
  my_row_index <- which(coef_tab_total$theoretical_intensity == my_intensity & round(coef_tab_total$expect_disp_ind, digits = 1) == my_disp_ind)
  coef_tab <- coef_tab_total[my_row_index, ]
  
  # load the count matrices
  vect_matrices <- readRDS(file = my_file)
  
  ### loop over the populations
  for (num_pattern in vect_num_pattern) {

    # get the current population
    mymatrix <- vect_matrices[[num_pattern]]

    # get intensity
    current_intensity <- mean(mymatrix)
    
    # get dispersion index
    current_disp_index <- DispersionIndex(pop = mymatrix, nx = my_window$xrange[2], ny = my_window$yrange[2])
    
    
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
                          biais_random = result_tab$mean_random - current_intensity,
                          biais_syst = result_tab$mean_syst - current_intensity,
                          biais_BAS = result_tab$mean_BAS - current_intensity,
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
      write.csv2(result_tab, file = paste0("./output/text/", storage_folder, "/simul_result_tab_intensity_", my_intensity, 
                                           "_disp_", round(coef_tab$expect_disp_ind, digits = 1), "_sample-size_", sample_size,
                                           "_num_pattern_", num_pattern, ".csv"), row.names = FALSE)
      
    } # end of loop over sample_size
    
  } # end of loop over num_pattern 
  
} # end of loop over my_file


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
merged_data <- merged_data[with(merged_data, order(theoretical_intensity, expect_disp_ind, num_pattern, sample_size)), ]

# join the cluster diameter
merged_data <- cbind(merged_data, cluster_diam = merged_data$my_scale * 2)

# save the dataframe
write.csv2(merged_data, file = "./data/tidy/RUN_2022-07-27_variance_estimator_ZIP_results_int1-10-20_nsim600_20_patterns_merged.csv", row.names = FALSE)

# make a version of the dataframe averaged over the 10 patterns
merged_data <- merged_data %>% select(-num_pattern)

merged_data_averaged <- merged_data %>% 
  group_by(pattern, theoretical_intensity, sample_size) %>% 
  summarise(across(my_kappa:cluster_diam, ~ mean(.x, na.rm = TRUE)))

# save the dataframe
write.csv2(merged_data_averaged, file = "./data/tidy/RUN_2022-07-27_variance_estimator_ZIP_results_int1-10-20_nsim600_20_patterns_merged_averaged.csv", row.names = FALSE)


