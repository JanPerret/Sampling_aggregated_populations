#######################################################
# 
# Sampling simulation of spatially aggregated populations
# 
# 2_pattern_creation_Matern_modalities_2_3.R
# 
# jan.perret@gmail.com
#######################################################

###
### script objective : generate the Matern and Poisson patterns on which the sampling surveys will be simulated - for modalities 2 and 3
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
nsim = 100
my_window <- owin(c(0, 100), c(0, 100))
intensity_vect <- c(1, 5, 10, 20)
disp_index_vect <- seq(from = 2, to = 74, by = 2)


### load parameter values for the Matern process
coef_tab_modality_2 <- read.csv2("./data/tidy/Matern_parameters_coef_tab_total_2021-06-24_kappa.csv")
coef_tab_modality_2 <- subset(coef_tab_modality_2, coef_tab_modality_2$theoretical_intensity %in% intensity_vect & 
                           round(coef_tab_modality_2$expect_disp_ind, digits = 1) %in% disp_index_vect)

coef_tab_modality_3 <- read.csv2("./data/tidy/Matern_parameters_coef_tab_total_2021-06-24_mu.csv")
coef_tab_modality_3 <- subset(coef_tab_modality_3, coef_tab_modality_3$theoretical_intensity %in% intensity_vect & 
                                round(coef_tab_modality_3$expect_disp_ind, digits = 1) %in% disp_index_vect)

# set to the number of cores to use
registerDoParallel(2)

# store and print start time
start.time_total <- Sys.time()
start.time_total



### MODALITY 2
# folder where the patterns will be saved
storage_folder <- "RUN_2022-06-08-simul_aggregation_intensity_patterns_modality_2"

# set the random seed for result reproducibility
set.seed(20220608)

### generate the patterns with the Matern process
# loop over all parameter values available in coef_tab for the current intensity level
# for (j in 1:nrow(coef_tab_modality_2)) {
foreach (j = 1:nrow(coef_tab_modality_2), .packages = "spatstat") %dopar% {
  
  # # print progress
  # cat(paste("Intensity =", coef_tab_modality_2$theoretical_intensity[j],
  #           "; expected dispersion index =", round(coef_tab_modality_2$expect_disp_ind[j], digits = 1),
  #           "started at", Sys.time(), "\n"))
  
  # create list for matrices
  matrix_list <- vector(mode = "list", length = nsim)
  
  
  # make nsim iterations
  for (i in 1:nsim) {
    
    # simulate the pattern
    mypattern <- rMatClust(kappa = coef_tab_modality_2$my_kappa[j], scale = coef_tab_modality_2$my_scale[j],
                           mu = coef_tab_modality_2$my_mu[j], win = my_window)
    
    # /!\ Sometimes rMatClust() creates points falling outside of the simulation window. 
    # Either the x or y coordinate is slightly outside the window's limits (e.g. the 9th decimal place).
    # To see the difference : sprintf("%.10f", mypattern_bug$y[which(mypattern_bug$y > 100)])
    # When such a point pattern occurs, the function quadratcount() crashes.
    # So I added this step below to remove the points falling outside the window :
    if (any(mypattern$x < 0)) {mypattern <- subset(mypattern, !mypattern$x < 0)}
    if (any(mypattern$x > my_window$xrange[2])) {mypattern <- subset(mypattern, !mypattern$x > my_window$xrange[2])}
    if (any(mypattern$y < 0)) {mypattern <- subset(mypattern, !mypattern$y < 0)}
    if (any(mypattern$y > my_window$yrange[2])) {mypattern <- subset(mypattern, !mypattern$y > my_window$yrange[2])}
    
    # create the matrix with the counts
    mymatrix <- quadratcount(mypattern, nx = my_window$xrange[2], ny = my_window$yrange[2])
    
    # add current count matrix to list
    matrix_list[[i]] <- mymatrix
    
  } # i
  
  # save the count matrices
  saveRDS(matrix_list, file = paste0("./output/text/", storage_folder,
                                     "/matrices_intensity-", coef_tab_modality_2$theoretical_intensity[j],
                                     "_disp_index-", round(coef_tab_modality_2$expect_disp_ind[j], digits = 1),
                                     "_modality_2.rds"))
  
} # j




### MODALITY 3
# folder where the patterns will be saved
storage_folder <- "RUN_2022-06-08-simul_aggregation_intensity_patterns_modality_3"

# set the random seed for result reproducibility
set.seed(20220608)

### generate the patterns with the Matern process
# loop over all parameter values available in coef_tab for the current intensity level
# for (j in 1:nrow(coef_tab_modality_3)) {
foreach (j = 1:nrow(coef_tab_modality_3), .packages = "spatstat") %dopar% {
  
  # # print progress
  # cat(paste("Intensity =", coef_tab_modality_3$theoretical_intensity[j],
  #           "; expected dispersion index =", round(coef_tab_modality_3$expect_disp_ind[j], digits = 1),
  #           "started at", Sys.time(), "\n"))
  
  # create list for matrices
  matrix_list <- vector(mode = "list", length = nsim)
  
  
  # make nsim iterations
  for (i in 1:nsim) {
    
    # simulate the pattern
    mypattern <- rMatClust(kappa = coef_tab_modality_3$my_kappa[j], scale = coef_tab_modality_3$my_scale[j],
                           mu = coef_tab_modality_3$my_mu[j], win = my_window)
    
    # /!\ Sometimes rMatClust() creates points falling outside of the simulation window. 
    # Either the x or y coordinate is slightly outside the window's limits (e.g. the 9th decimal place).
    # To see the difference : sprintf("%.10f", mypattern_bug$y[which(mypattern_bug$y > 100)])
    # When such a point pattern occurs, the function quadratcount() crashes.
    # So I added this step below to remove the points falling outside the window :
    if (any(mypattern$x < 0)) {mypattern <- subset(mypattern, !mypattern$x < 0)}
    if (any(mypattern$x > my_window$xrange[2])) {mypattern <- subset(mypattern, !mypattern$x > my_window$xrange[2])}
    if (any(mypattern$y < 0)) {mypattern <- subset(mypattern, !mypattern$y < 0)}
    if (any(mypattern$y > my_window$yrange[2])) {mypattern <- subset(mypattern, !mypattern$y > my_window$yrange[2])}
    
    # create the matrix with the counts
    mymatrix <- quadratcount(mypattern, nx = my_window$xrange[2], ny = my_window$yrange[2])
    
    # add current count matrix to list
    matrix_list[[i]] <- mymatrix
    
  } # i
  
  # save the count matrices
  saveRDS(matrix_list, file = paste0("./output/text/", storage_folder,
                                     "/matrices_intensity-", coef_tab_modality_3$theoretical_intensity[j],
                                     "_disp_index-", round(coef_tab_modality_3$expect_disp_ind[j], digits = 1),
                                     "_modality_3.rds"))
  
} # j





# close the cluster
stopImplicitCluster()


# print total execution time
end.time_total <- Sys.time()
time.taken_total <- difftime(end.time_total, start.time_total, units = "mins")
end.time_total
time.taken_total

