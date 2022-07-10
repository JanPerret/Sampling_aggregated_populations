#########################################################
# 
# Sampling simulation of spatially aggregated populations
# 
# 3_sampling_simulation_natural_populations.R
# 
# jan.perret@gmail.com
#########################################################

###
### script objective : sampling survey simulation for the 3 natural populations that we have entirely mapped
###

# clean workspace
rm(list = ls())

# load packages
library(ggplot2)
library(ggpubr)
library(spatstat)
library(dplyr)
library(tidyr)

# load functions
files.sources <- list.files("./R", full.names = TRUE)
invisible(sapply(files.sources, source))


# load the data from the 3 natural populations
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

# pivot to long table format
matrix_SANMIN_long <- as.data.frame.table(as.matrix(matrix_SANMIN_new), responseName = "count") %>%
  mutate_if(is.factor, as.integer)
colnames(matrix_SANMIN_long) <- c("y", "x", "count")
matrix_LIMGIR_long <- as.data.frame.table(as.matrix(matrix_LIMGIR_new), responseName = "count") %>%
  mutate_if(is.factor, as.integer)
colnames(matrix_LIMGIR_long) <- c("y", "x", "count")
matrix_BELSYL_long <- as.data.frame.table(as.matrix(matrix_BELSYL_new), responseName = "count") %>%
  mutate_if(is.factor, as.integer)
colnames(matrix_BELSYL_long) <- c("y", "x", "count")

# save data for the heatmaps
write.csv2(matrix_SANMIN_long, file = "./data/tidy/matrix_SANMIN_long.csv")
write.csv2(matrix_LIMGIR_long, file = "./data/tidy/matrix_LIMGIR_long.csv")
write.csv2(matrix_BELSYL_long, file = "./data/tidy/matrix_BELSYL_long.csv")



# function to compute variance without the n-1 correction
var_N <- function(x){var(x)*(length(x)-1)/length(x)}

# sample sizes to test (we selected only sample sizes with entire square roots)
sample_sizes_vect <- c(9, 16, 25, 36, 49, 64, 81, 100, 121, 144, 169, 196, 225, 256, 289, 324, 361, 400,
                       441, 484, 529, 576, 625, 676, 729, 784, 841, 900, 961, 1024, 1089, 1156, 1225,
                       1296, 1369, 1444, 1521, 1600, 1681, 1764, 1849, 1936, 2025, 2116, 2209, 2304,
                       2401, 2500, 2601, 2704, 2809, 2916, 3025, 3136, 3249, 3364, 3481, 3600)

# number of iterations for the simulations
nsim = 500


### get SRS sampling variance for every sample size 
SRS_variance_vect_SANMIN <- c(rep(NA, times = length(sample_sizes_vect)))
SRS_variance_vect_LIMGIR <- c(rep(NA, times = length(sample_sizes_vect)))
SRS_variance_vect_BELSYL <- c(rep(NA, times = length(sample_sizes_vect)))

SRS_list_means_SANMIN <- vector(mode = "list", length = length(sample_sizes_vect))
SRS_list_means_LIMGIR <- vector(mode = "list", length = length(sample_sizes_vect))
SRS_list_means_BELSYL <- vector(mode = "list", length = length(sample_sizes_vect))

for (i in 1:length(sample_sizes_vect)) {
  
  SRS_samples_SANMIN <- list()
  SRS_samples_LIMGIR <- list()
  SRS_samples_BELSYL <- list()
  
  for (j in 1:nsim) {
    
    SRS_samples_SANMIN[[j]] <- sample_random(pop = matrix_SANMIN_new, nx = 100, ny = 100, sample_size = sample_sizes_vect[i])
    SRS_samples_LIMGIR[[j]] <- sample_random(pop = matrix_LIMGIR_new, nx = 100, ny = 100, sample_size = sample_sizes_vect[i])
    SRS_samples_BELSYL[[j]] <- sample_random(pop = matrix_BELSYL_new, nx = 100, ny = 100, sample_size = sample_sizes_vect[i])
    
  }
  
  # get the means of every sample
  means_SANMIN <- sapply(SRS_samples_SANMIN, mean)
  means_LIMGIR <- sapply(SRS_samples_LIMGIR, mean)
  means_BELSYL <- sapply(SRS_samples_BELSYL, mean)
  
  # store the means
  SRS_list_means_SANMIN[[i]] <- means_SANMIN
  SRS_list_means_LIMGIR[[i]] <- means_LIMGIR
  SRS_list_means_BELSYL[[i]] <- means_BELSYL
  
  # get sampling variance
  SRS_variance_vect_SANMIN[i] <- var_N(means_SANMIN)
  SRS_variance_vect_LIMGIR[i] <- var_N(means_LIMGIR)
  SRS_variance_vect_BELSYL[i] <- var_N(means_BELSYL)
  
}

# name vector elements with the corresponding sample size
names(SRS_variance_vect_SANMIN) <- as.character(sample_sizes_vect)
names(SRS_variance_vect_LIMGIR) <- as.character(sample_sizes_vect)
names(SRS_variance_vect_BELSYL) <- as.character(sample_sizes_vect)



### get SBS sampling variance for every sample size 
SBS_variance_vect_SANMIN <- c(rep(NA, times = length(sample_sizes_vect)))
SBS_variance_vect_LIMGIR <- c(rep(NA, times = length(sample_sizes_vect)))
SBS_variance_vect_BELSYL <- c(rep(NA, times = length(sample_sizes_vect)))

SBS_list_means_SANMIN <- vector(mode = "list", length = length(sample_sizes_vect))
SBS_list_means_LIMGIR <- vector(mode = "list", length = length(sample_sizes_vect))
SBS_list_means_BELSYL <- vector(mode = "list", length = length(sample_sizes_vect))

for (i in 1:length(sample_sizes_vect)) {
  
  SBS_samples_SANMIN <- list()
  SBS_samples_LIMGIR <- list()
  SBS_samples_BELSYL <- list()
  
  for (j in 1:nsim) {
    
    SBS_samples_SANMIN[[j]] <- sample_BAS(pop = matrix_SANMIN_new, nx = 100, ny = 100, sample_size = sample_sizes_vect[i])
    SBS_samples_LIMGIR[[j]] <- sample_BAS(pop = matrix_LIMGIR_new, nx = 100, ny = 100, sample_size = sample_sizes_vect[i])
    SBS_samples_BELSYL[[j]] <- sample_BAS(pop = matrix_BELSYL_new, nx = 100, ny = 100, sample_size = sample_sizes_vect[i])
    
  }
  
  # get the means of every sample
  means_SANMIN <- sapply(SBS_samples_SANMIN, mean)
  means_LIMGIR <- sapply(SBS_samples_LIMGIR, mean)
  means_BELSYL <- sapply(SBS_samples_BELSYL, mean)
  
  # store the means
  SBS_list_means_SANMIN[[i]] <- means_SANMIN
  SBS_list_means_LIMGIR[[i]] <- means_LIMGIR
  SBS_list_means_BELSYL[[i]] <- means_BELSYL
  
  # get sampling variance
  SBS_variance_vect_SANMIN[i] <- var_N(means_SANMIN)
  SBS_variance_vect_LIMGIR[i] <- var_N(means_LIMGIR)
  SBS_variance_vect_BELSYL[i] <- var_N(means_BELSYL)
  
}

# name vector elements with the corresponding sample size
names(SBS_variance_vect_SANMIN) <- as.character(sample_sizes_vect)
names(SBS_variance_vect_LIMGIR) <- as.character(sample_sizes_vect)
names(SBS_variance_vect_BELSYL) <- as.character(sample_sizes_vect)


### get SYS sampling variance for every sample size 
SYS_variance_vect_SANMIN <- c(rep(NA, times = length(sample_sizes_vect)))
SYS_variance_vect_LIMGIR <- c(rep(NA, times = length(sample_sizes_vect)))
SYS_variance_vect_BELSYL <- c(rep(NA, times = length(sample_sizes_vect)))

SYS_list_means_SANMIN <- vector(mode = "list", length = length(sample_sizes_vect))
SYS_list_means_LIMGIR <- vector(mode = "list", length = length(sample_sizes_vect))
SYS_list_means_BELSYL <- vector(mode = "list", length = length(sample_sizes_vect))

for (i in 1:length(sample_sizes_vect)) {
  
  SYS_samples_SANMIN <- list()
  SYS_samples_LIMGIR <- list()
  SYS_samples_BELSYL <- list()
  
  for (j in 1:nsim) {
    
    SYS_samples_SANMIN[[j]] <- sample_systematic_fixed_size(pop = matrix_SANMIN_new, nx = 100, ny = 100, sample_size = sample_sizes_vect[i])
    SYS_samples_LIMGIR[[j]] <- sample_systematic_fixed_size(pop = matrix_LIMGIR_new, nx = 100, ny = 100, sample_size = sample_sizes_vect[i])
    SYS_samples_BELSYL[[j]] <- sample_systematic_fixed_size(pop = matrix_BELSYL_new, nx = 100, ny = 100, sample_size = sample_sizes_vect[i])
    
  }
  
  # get the means of every sample
  means_SANMIN <- sapply(SYS_samples_SANMIN, mean)
  means_LIMGIR <- sapply(SYS_samples_LIMGIR, mean)
  means_BELSYL <- sapply(SYS_samples_BELSYL, mean)
  
  # store the means
  SYS_list_means_SANMIN[[i]] <- means_SANMIN
  SYS_list_means_LIMGIR[[i]] <- means_LIMGIR
  SYS_list_means_BELSYL[[i]] <- means_BELSYL
  
  # get sampling variance
  SYS_variance_vect_SANMIN[i] <- var_N(means_SANMIN)
  SYS_variance_vect_LIMGIR[i] <- var_N(means_LIMGIR)
  SYS_variance_vect_BELSYL[i] <- var_N(means_BELSYL)
  
}

# name vector elements with the corresponding sample size
names(SYS_variance_vect_SANMIN) <- as.character(sample_sizes_vect)
names(SYS_variance_vect_LIMGIR) <- as.character(sample_sizes_vect)
names(SYS_variance_vect_BELSYL) <- as.character(sample_sizes_vect)



### compute the variance ratios
# SYS
vect_SYS_var_ratio_SANMIN <- round(SYS_variance_vect_SANMIN / SRS_variance_vect_SANMIN, digits = 3)
vect_SYS_var_ratio_LIMGIR <- round(SYS_variance_vect_LIMGIR / SRS_variance_vect_LIMGIR, digits = 3)
vect_SYS_var_ratio_BELSYL <- round(SYS_variance_vect_BELSYL / SRS_variance_vect_BELSYL, digits = 3)

# SBS
vect_SBS_var_ratio_SANMIN <- round(SBS_variance_vect_SANMIN / SRS_variance_vect_SANMIN, digits = 3)
vect_SBS_var_ratio_LIMGIR <- round(SBS_variance_vect_LIMGIR / SRS_variance_vect_LIMGIR, digits = 3)
vect_SBS_var_ratio_BELSYL <- round(SBS_variance_vect_BELSYL / SRS_variance_vect_BELSYL, digits = 3)

df_natural_pop_var_ratio <- data.frame(rbind(vect_SYS_var_ratio_SANMIN, vect_SYS_var_ratio_LIMGIR, vect_SYS_var_ratio_BELSYL,
                                             vect_SBS_var_ratio_SANMIN, vect_SBS_var_ratio_LIMGIR, vect_SBS_var_ratio_BELSYL))
df_natural_pop_var_ratio <- cbind(species = c("SANMIN", "LIMGIR", "BELSYL", "SANMIN", "LIMGIR", "BELSYL"),
                                  sampling_method = c("SYS", "SYS", "SYS", "SBS", "SBS", "SBS"),
                                  df_natural_pop_var_ratio)
tibble_natural_pop_var_ratio <- tibble(df_natural_pop_var_ratio)

colnames(tibble_natural_pop_var_ratio) <- c("species", "sampling_method", "9", "16", "25", "36",
                                            "49", "64", "81", "100", "121", "144", "169", "196",
                                            "225", "256", "289", "324", "361", "400", "441", "484",
                                            "529", "576", "625", "676", "729", "784", "841", "900",
                                            "961", "1024", "1089", "1156", "1225", "1296", "1369",
                                            "1444", "1521", "1600", "1681", "1764", "1849", "1936",
                                            "2025", "2116", "2209", "2304", "2401", "2500", "2601",
                                            "2704", "2809", "2916", "3025", "3136", "3249", "3364",
                                            "3481", "3600")

tibble_natural_pop_var_ratio_long <- tibble_natural_pop_var_ratio %>%
  pivot_longer(cols = !c(species, sampling_method),
               names_to = "sample_size", values_to = "var_ratio")
tibble_natural_pop_var_ratio_long$sample_size <- as.integer(tibble_natural_pop_var_ratio_long$sample_size)


# save data
write.csv2(tibble_natural_pop_var_ratio_long, file = "./data/tidy/variance_ratio_table_natural_populations_60_sample_sizes.csv", row.names = FALSE)


