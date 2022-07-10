#######################################################
# 
# Sampling simulation of spatially aggregated populations
# 
# 2_pattern_creation_Matern_modality_1.R
# 
# jan.perret@gmail.com
#######################################################

###
### script objective : generate the Matern and Poisson patterns on which the sampling surveys will be simulated - for modality 1
###

### clean workspace
rm(list = ls())

### load packages
source("./R/imports.R")

### load functions
files.sources <- list.files("./R", full.names = TRUE)
invisible(sapply(files.sources, source))

### general settings
nsim = 200
my_window <- owin(c(0, 100), c(0, 100))
intensity_vect <- c(1, 10, 20)
disp_index_vect <- seq(from = 2, to = 75, by = 1)


### load parameter values for the Matern process
coef_tab_total <- read.csv2("./data/tidy/Matern_parameters_coef_tab_total_2021-06-13.csv")
coef_tab_total <- subset(coef_tab_total, coef_tab_total$theoretical_intensity %in% intensity_vect & 
                           round(coef_tab_total$expect_disp_ind, digits = 1) %in% disp_index_vect)

# folder where the patterns will be saved
storage_folder <- "RUN_2022-06-08-simul_aggregation_intensity_patterns_modality_1"


# store and print start time
start.time_total <- Sys.time()
start.time_total

# set the random seed for result reproducibility
set.seed(20220608)

### generate the patterns with the Matern process
# loop over all parameter values available in coef_tab for the current intensity level
for (j in 1:nrow(coef_tab_total)) {
  
  # print progress
  cat(paste("Intensity =", coef_tab_total$theoretical_intensity[j],
            "; expected dispersion index =", round(coef_tab_total$expect_disp_ind[j], digits = 1),
            "started at", Sys.time(), "\n"))
  
  # create list for patterns
  pattern_list <- vector(mode = "list", length = nsim)
  
  # create list for matrices
  matrix_list <- vector(mode = "list", length = nsim)
  
  
  # make nsim iterations
  for (i in 1:nsim) {
    
    # simulate the pattern
    mypattern <- rMatClust(kappa = coef_tab_total$my_kappa[j], scale = coef_tab_total$my_scale[j],
                           mu = coef_tab_total$my_mu[j], win = my_window)
    
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
    
    # # add current pattern to list
    # pattern_list[[i]] <- mypattern
    
    # add current count matrix to list
    matrix_list[[i]] <- mymatrix
    
  } # i
  
  # # save the point patterns
  # saveRDS(pattern_list, file = paste0("./output/text/", storage_folder,
  #                                     "/patterns_intensity-", coef_tab_total$theoretical_intensity[j],
  #                                     "_disp_index-", round(coef_tab_total$expect_disp_ind[j], digits = 1),
  #                                     ".rds"))
  
  # save the count matrices
  saveRDS(matrix_list, file = paste0("./output/text/", storage_folder,
                                     "/matrices_intensity-", coef_tab_total$theoretical_intensity[j],
                                     "_disp_index-", round(coef_tab_total$expect_disp_ind[j], digits = 1),
                                     ".rds"))
  
} # j



### generate the patterns with the Poisson process
# set the random seed for result reproducibility
set.seed(20220528)

# loop over all intensities in intensity_vect
for (k in 1:length(intensity_vect)) {
  
  # print progress
  cat(paste("Intensity =", intensity_vect[k],
            "; expected dispersion index = 1 (poisson process)",
            "started at", Sys.time(), "\n"))
  
  # # create list for patterns
  # pattern_list <- vector(mode = "list", length = nsim)
  
  # create list for matrices
  matrix_list <- vector(mode = "list", length = nsim)
  
  
  # make nsim iterations
  for (i in 1:nsim) {
    
    # simulate the pattern
    mypattern <- rpoispp(lambda = intensity_vect[k], win = my_window)
    
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
    
    # # add current pattern to list
    # pattern_list[[i]] <- mypattern
    
    # add current count patrix to list
    matrix_list[[i]] <- mymatrix
    
  } # i
  
  # # save the point patterns
  # saveRDS(pattern_list, file = paste0("./output/text/", storage_folder,
  #                                     "/patterns_intensity-", intensity_vect[k],
  #                                     "_disp_index-1",
  #                                     ".rds"))
  
  # save the count matrices
  saveRDS(matrix_list, file = paste0("./output/text/", storage_folder,
                                     "/matrices_intensity-", intensity_vect[k],
                                     "_disp_index-1",
                                     ".rds"))
  
} # k


# print total execution time
end.time_total <- Sys.time()
time.taken_total <- difftime(end.time_total, start.time_total, units = "mins")
end.time_total
time.taken_total

