#######################################################
# 
# Sampling simulation of spatially aggregated populations
# 
# 2_pattern_creation_SSI.R
# 
# jan.perret@gmail.com
#######################################################

###
### script objective : generate the SSI patterns on which the sampling surveys will be simulated
###

### clean workspace
rm(list = ls())

### load packages
source("./R/imports.R")

### load functions
files.sources <- list.files("./R", full.names = TRUE)
invisible(sapply(files.sources, source)) # invisible is to avoid sapply() to print anything for the cluster

### general settings
nsim = 200
my_window <- owin(c(0, 100), c(0, 100)) # window in which we are going to simulate the populations


### load parameter values for the SSI process
coef_tab_total <- read.csv2("./data/tidy/SSI_parameters_coef_tab_total_2021-10-06.csv")
coef_tab_total <- subset(coef_tab_total, round(coef_tab_total$expect_disp_ind, digits = 1) %in% c(0.4, 0.6, 0.8))
coef_tab_total <- subset(coef_tab_total, coef_tab_total$theoretical_intensity %in% c(1, 10, 20))


# folder where the patterns will be saved
storage_folder <- "RUN_2022-06-08-simul_aggregation_intensity_patterns_modality_1"


# store and print start time
start.time_total <- Sys.time()
start.time_total


# loop over all parameter values available in coef_tab for the current intensity level
for (j in 1:nrow(coef_tab_total)) {
  
  # print progress
  cat(paste("Intensity =", coef_tab_total$theoretical_intensity[j],
            "; expected dispersion index =", round(coef_tab_total$expect_disp_ind[j], digits = 1),
            "started at", Sys.time(), "\n"))
  
  # # create list for patterns
  # pattern_list <- vector(mode = "list", length = nsim)
  
  # create list for matrices
  matrix_list <- vector(mode = "list", length = nsim)
  
  
  # make nsim iterations
  for (i in 1:nsim) {
    
    mypattern <- rSSI(r = coef_tab_total$my_r[j],
                      n = coef_tab_total$theoretical_intensity[j] * my_window$xrange[2] * my_window$yrange[2], # wanted intensity * area of the simulation window
                      win = my_window, giveup = 2000)
    
    # /!\ Sometimes rMatClust() creates points falling outside of the simulation window. 
    # Either the x or y coordinate is slightly outside the window's limits (e.g. the 9th decimal place).
    # To see the difference : sprintf("%.10f", mypattern_bug$y[which(mypattern_bug$y > 100)])
    # When such a point pattern occurs, the function quadratcount() crashes.
    # So I added this step below to remove the points falling outside the window :
    if (any(mypattern$x < 0)) {mypattern <- subset(mypattern, !mypattern$x < 0)}
    if (any(mypattern$x > 100)) {mypattern <- subset(mypattern, !mypattern$x > 100)}
    if (any(mypattern$y < 0)) {mypattern <- subset(mypattern, !mypattern$y < 0)}
    if (any(mypattern$y > 100)) {mypattern <- subset(mypattern, !mypattern$y > 100)}
    
    # create the matrix with the counts
    mymatrix <- quadratcount(mypattern, nx = 100, ny = 100)
    
    # # add current pattern to list
    # pattern_list[[i]] <- mypattern
    
    # add current count matrix to list
    matrix_list[[i]] <- mymatrix
    
  } # i
  
  # # save the point patterns
  # saveRDS(pattern_list, file = paste0("./output/text/", storage_folder,
  #                                     "/patterns_intensity-", coef_tab_total$theoretical_intensity[j],
  #                                     "_disp_index-", round(coef_tab_total$expect_disp_ind[j], digits = 1),
  #                                     "_script1.rds"))
  
  # save the count matrices
  saveRDS(matrix_list, file = paste0("./output/text/", storage_folder,
                                     "/matrices_intensity-", coef_tab_total$theoretical_intensity[j],
                                     "_disp_index-", round(coef_tab_total$expect_disp_ind[j], digits = 1),
                                     "_script1.rds"))
  
} # j


# print total execution time
end.time_total <- Sys.time()
time.taken_total <- difftime(end.time_total, start.time_total, units = "mins")
end.time_total
time.taken_total


