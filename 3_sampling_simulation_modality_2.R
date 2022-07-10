#######################################################
# 
# Sampling simulation of spatially aggregated populations
# 
# 3_sampling_simulation_modality_2.R
# 
# jan.perret@gmail.com
#######################################################

###
### script objective : simulate sampling surveys on previously generated patterns to estimate the true variance of various sampling designs - for modality 2
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

### general settings
vect_num_pattern = c(1:60)
nsim = 1000
my_window <- owin(c(0, 100), c(0, 100))
window_size <- area(my_window)
sample_size_vect = c(9, 15, 25, 49, 100, 150, 196, 300, 400)

# load parameter values for the 3 point processes
coef_tab_Matern <- read.csv2("./data/tidy/Matern_parameters_coef_tab_total_2021-06-24_kappa.csv") # modality 2
coef_tab_Matern <- cbind(coef_tab_Matern[, 1:2], process_type = "Matern_cluster", coef_tab_Matern[, 3:5], my_r = NA, expect_disp_ind = coef_tab_Matern[, 6])
coef_tab_total <- coef_tab_Matern

# add a column with the theoretical intensity value written as character
my_intensity <- paste0("intensity = ", coef_tab_total$theoretical_intensity)
coef_tab_total <- cbind(my_intensity = my_intensity, coef_tab_total)


# set the folder to save the outputs
storage_folder <- "RUN_2022-06-08_simul_aggregation_intensity_results_modality_2"

# create a vector with the names of all files in the directory containing the previously generated count matrices
base_dir <- "./output/text/RUN_2022-06-08-simul_aggregation_intensity_patterns_modality_2/"
vect_all_files <- paste0(base_dir, list.files(base_dir, recursive = TRUE))

# store and print start time
start.time_total <- Sys.time()
start.time_total

# set the number of cores to use
registerDoParallel(2)

# set the random seed for result reproducibility
# set.seed(20220608)


### loop over the files containing the count matrices in vect_all_files
# for (my_file in vect_all_files) {
foreach (my_file = vect_all_files, .packages = "spatstat") %dopar% {
  
  # # print progress
  # cat(paste(my_file, "; started at", Sys.time(), "\n"))
  
  # get the intensity and aggregation level of the current pattern
  vect_file_name <- as.vector(strsplit(my_file, split = "[-_.]+"))[[1]]
  my_intensity <- as.numeric(vect_file_name[13])
  my_disp_ind <- as.numeric(vect_file_name[16])
  
  # get point process parameter values for the current intensity and aggregation levels
  my_row_index <- which(coef_tab_total$theoretical_intensity == my_intensity & round(coef_tab_total$expect_disp_ind, digits = 1) == my_disp_ind)
  coef_tab <- coef_tab_total[my_row_index, ]
  
  # load the count matrices
  vect_matrices <- readRDS(file = my_file)
  
  
  ### loop over the previously generated patterns
  for (num_pattern in vect_num_pattern) {
    
    # get the current count matrix
    mymatrix <- vect_matrices[[num_pattern]]
    
    # get intensity
    current_intensity <- mean(mymatrix)
    
    # get dispersion index
    current_disp_index <- DispersionIndex(pop = mymatrix, nx = my_window$xrange[2], ny = my_window$yrange[2])
    
    
    ### loop over sample sizes
    for (sample_size in sample_size_vect) {
      
      # # print progress
      # cat(paste("sample size =", sample_size, "; num pattern =", num_pattern, "; started at", Sys.time(), "\n"))
      
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
      
      
      # make nsim iterations
      for (i in 1:nsim) {
        
        # if mypattern is empty we assign the output values and go to the next iteration
        if (current_intensity == 0) {
          result_random[i, 1] <- 0
          result_syst[i, 1] <- 0
          result_BAS[i, 1] <- 0
          prop_zero_random[i, 1] <- 1
          prop_zero_syst[i, 1] <- 1
          prop_zero_BAS[i, 1] <- 1
          next
        }
        
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
        
        # compute the textbook variance estimator for the 3 sampling methods
        estim_vtextbook_random[i, 1] <- var_estim_vtextbook(y = mysample_random, n = sample_size, N = window_size)
        estim_vtextbook_syst[i, 1] <- var_estim_vtextbook(y = mysample_syst, n = sample_size, N = window_size)
        estim_vtextbook_BAS[i, 1] <- var_estim_vtextbook(y = mysample_BAS, n = sample_size, N = window_size)
        
        
      } # i
      
      
      # calculate the synthetic indicators and fill the result table
      result_tab <- cbind(coef_tab,
                          num_pattern = paste0("pattern_", num_pattern),
                          nsim = nsim,
                          sample_size = sample_size,
                          intensity = current_intensity,
                          pop_size = (current_intensity * area(my_window)),
                          dispersion_index = current_disp_index)
      
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
                          mean_estim_vtextbook_BAS = mean(estim_vtextbook_BAS[, 1]))
      
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
      
      # compute the coverage obtained with the textbook variance estimator
      coverage_vtextbook_random <- compute_coverage(Ybars = c(rep(current_intensity, times = nsim)), ybars = result_random[, 1], estim_vars = estim_vtextbook_random[, 1], n = sample_size, N = window_size)
      coverage_vtextbook_syst <- compute_coverage(Ybars = c(rep(current_intensity, times = nsim)), ybars = result_syst[, 1], estim_vars = estim_vtextbook_syst[, 1], n = sample_size, N = window_size)
      coverage_vtextbook_BAS <- compute_coverage(Ybars = c(rep(current_intensity, times = nsim)), ybars = result_BAS[, 1], estim_vars = estim_vtextbook_BAS[, 1], n = sample_size, N = window_size)
      
      # join the coverage to the result table
      result_tab <- cbind(result_tab,
                          coverage_vtextbook_random = coverage_vtextbook_random,
                          coverage_vtextbook_syst = coverage_vtextbook_syst,
                          coverage_vtextbook_BAS = coverage_vtextbook_BAS)
      
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

