#######################################################
# 
# Sampling simulation of spatially aggregated populations
# 
# 1_estimate_point_process_parameters_Matern_modality_2.R
# 
# jan.perret@gmail.com
#######################################################

###
### script objective : estimate the Matern process parameters to get the wanted levels of dispersion index - for modality 2 (intensity increased by increasing only kappa)
###

### clean workspace
rm(list = ls())

### load packages
source("./R/imports.R")

### load functions
files.sources <- list.files("./R", full.names = TRUE)
sapply(files.sources, source)


### get Matern cluster process parameters values to get the following phase space :
# intensity from 1 to 30
# dispersion index from 1 to 75
# step : every unit

# Reminder of the parameters of the Matern process :
# kappa = intensity of the Poisson process of cluster centers
# scale = radius of the clusters
#    mu = mean number of points per cluster
# Expected intensity = kappa * mu

# general settings :
wanted_intensities = c(1:30) # vector with the wanted levels of intensity
wanted_disp_index = c(1:75) # vector with the wanted levels of dispersion index
my_window = owin(c(0, 100), c(0, 100)) # window in which the populations will be simulated
nsim = 500 # the number of iterations to estimate parameters values


# get the values of kappa for each level of intensity given that the value of mu is a constant = 100
mu_vect = rep(100, times = 30)
kappa_vect = wanted_intensities / mu_vect
intensity_vect = kappa_vect * mu_vect

# get the vector of initial values for scale parameter
scale_vect = c(seq(from = 0.1, to = 30, by = 0.1))
# with scale = 0.5 each cluster approximately fits in one 1x1 cell
# with scale = 50 each cluster covers all the simulation window,
# but I choose scale = 30 as limit value because that's enough
# for the gradient of dispersion index of interest (there is no level
# of intensity for which the dispersion index exceeds 2 for scale > 30)

# create simulation input table
input_table <- data.frame(cbind(
  intensity = rep(intensity_vect, each = length(scale_vect)),
  kappa = rep(kappa_vect, each = length(scale_vect)),
  mu = rep(mu_vect, each = length(scale_vect)),
  scale = rep(scale_vect, times = length(intensity_vect)),
  dispersion_index = NA))


# start execution timer
start.time <- Sys.time()

# loop over values of kappa and mu
for (j in 1:nrow(input_table)) {
  
  cat(paste("input_table line", j, "/", nrow(input_table)))
  cat("\n")
  
  # create an empty vector to store the values of dispersion index of the nsim simulations
  disp_index_vect <- c(rep(NA, times = nsim))
  
  # make the nsim simulations
  for (i in 1:nsim) {
    
    # create the pattern
    mypattern <- rMatClust(kappa = input_table$kappa[j], scale = input_table$scale[j],
                           mu = input_table$mu[j], win = my_window)
    
    # /!\ Sometimes rMatClust() creates points falling outside of the simulation window.
    # Either the x or y coordinate is slightly outside the window's limits (e.g. the 9th decimal place).
    # To see the difference : sprintf("%.10f", mypattern_bug$y[which(mypattern_bug$y > 100)])
    # When such a point pattern occurs, the function quadratcount() crashes.
    # So I added this step below to remove the points falling outside the window :
    if (any(mypattern$x < 0)) {mypattern <- subset(mypattern, !mypattern$x < 0)}
    if (any(mypattern$x > 100)) {mypattern <- subset(mypattern, !mypattern$x > 100)}
    if (any(mypattern$y < 0)) {mypattern <- subset(mypattern, !mypattern$y < 0)}
    if (any(mypattern$y > 100)) {mypattern <- subset(mypattern, !mypattern$y > 100)}
    
    # if mypattern is empty we assign the output values and go to the next iteration
    if (intensity(mypattern) == 0) {
      disp_index_vect[i] <- NaN
      next
    }
    
    # create the matrix with the counts
    mymatrix <- quadratcount(mypattern, nx = 100, ny = 100)
    
    # calculate the dispersion index
    disp_index_vect[i] <- DispersionIndex(pop = mymatrix, nx = 100, ny = 100)
    
    
  } # i
  
  # store the mean dispersion index
  input_table$dispersion_index[j] <- mean(disp_index_vect, na.rm = TRUE)
  
} # j

# store execution time
end.time <- Sys.time()
time.taken <- difftime(end.time, start.time, units = "mins")
time.taken

write.csv2(input_table, file = "./data/tidy/Matern_parameters_2021-06-24_kappa.csv", row.names = FALSE)



### refine the estimates
# load the table
Matern_parameters_tab <- read.csv2("./data/tidy/Matern_parameters_2021-06-24_kappa.csv")

# add the column row_status (1rst_step / 2nd_step / final_value)
Matern_parameters_tab <- cbind(Matern_parameters_tab, row_status = "1rst_step", stringsAsFactors = FALSE)

# remove rows with dispersion index below 1.6 because we won't need them
Matern_parameters_tab <- subset(Matern_parameters_tab, Matern_parameters_tab$dispersion_index > 1.6)

# check if the column dispersion_index is monotonic and if not, remove the corresponding rows
Matern_parameters_tab <- Matern_parameters_tab[with(Matern_parameters_tab, order(intensity, scale)), ]

# loop until all values remaining in the table are monotonic
while (!all(Matern_parameters_tab[with(Matern_parameters_tab, order(intensity, scale)), ] ==
            Matern_parameters_tab[with(Matern_parameters_tab, order(intensity, -dispersion_index)), ])) {
  
  # vector of rows to delete
  rows_to_delete <- c()
  
  for (i in 2:nrow(Matern_parameters_tab)) {
    
    if (Matern_parameters_tab$intensity[i-1] != Matern_parameters_tab$intensity[i]) {
      next
    }
    
    if (Matern_parameters_tab$dispersion_index[i-1] < Matern_parameters_tab$dispersion_index[i]) {
      cat(paste("Non-monotonic value at line", i))
      cat("\n")
      
      # store row indexes to delete
      rows_to_delete <- c(rows_to_delete, i)
      
    }
    
  } # i
  
  # remove the non monotonic values
  Matern_parameters_tab <- Matern_parameters_tab[-rows_to_delete, ]
  
} # end of the while loop


### settings for the next step
parameter_tab = Matern_parameters_tab
parameter_tab$dispersion_index <- as.numeric(parameter_tab$dispersion_index)
# wanted_intensities = c(1:30)
# wanted_disp_index = c(2:75)
my_window = owin(c(0, 100), c(0, 100))
nsim = 500

### check if all wanted values of dispersion index are between extreme values for every intensity level
# loop over levels of intensity
for (my_intensity in wanted_intensities) {
  
  # table with parameter values for current level of intensity
  current_intensity_tab <- subset(Matern_parameters_tab, Matern_parameters_tab$intensity == my_intensity)
  
  if(min(wanted_disp_index) < min(current_intensity_tab$dispersion_index)) {
    cat(paste("For intensity =", my_intensity, ", the lowest wanted value of dispersion index is outside the limits"))
    cat("\n")
  }
  
  if(max(wanted_disp_index) > max(current_intensity_tab$dispersion_index)) {
    cat(paste("For intensity =", my_intensity, ", the highest wanted value of dispersion index is outside the limits"))
    cat("\n")
  }
  
} # my_intensity



### Start the refining of scale estimates

# start execution timer
start.time <- Sys.time()

# start iteration counter
iteration_count <- 0

# table to store the values that couldn't be estimated
convergence_fails <- data.frame(intensity = NA, dispersion_index = NA)

# loop over levels of intensity
for (my_intensity in wanted_intensities) {
  
  # loop over levels of dispersion index
  for (my_disp_index in wanted_disp_index) {
    
    # print number of iterations needed
    cat(paste("Needed", iteration_count, "iterations"))
    cat("\n")
    iteration_count <- 0
    
    # print progress
    cat(paste("Intensity =", my_intensity, "; Dispersion index =", my_disp_index))
    cat("\n")
    
    # get row index for current intensity level
    int_rows <- which(parameter_tab$intensity == my_intensity)
    
    # loop while there is no value in parameter_tab equal to the wanted value of dispersion index
    while (!any(my_disp_index == round(parameter_tab[int_rows, ]$dispersion_index, digits = 1) )) {
      
      # count number of iterations needed to estimate the parameter
      iteration_count <- iteration_count + 1
      
      # break if the number of iterations becomes too high
      if (iteration_count > 100) {
        cat(paste("Intensity =", my_intensity, "; Dispersion index =", my_disp_index, "-----------> couldn't converge"))
        cat("\n")
        convergence_fails <- rbind(convergence_fails, c(my_intensity, my_disp_index))
        break
      }
      
      # get index of rows with dispersion index directly below and above the wanted value
      # index_below <- which.min(abs(parameter_tab[int_rows, ]$dispersion_index - my_disp_index)) + min(int_rows) - 1
      index_below <- Position(function(x) x < my_disp_index, parameter_tab[int_rows, ]$dispersion_index) + min(int_rows) - 1
      index_above <- index_below - 1
      
      ### refine the value of scale by iterative dichotomy
      below_scale <- parameter_tab$scale[index_below]
      above_scale <- parameter_tab$scale[index_above]
      
      # get the middle value between the 2 scale values
      scale_to_test <- (below_scale + above_scale) / 2
      
      # create an empty vector to store the values of dispersion index of the nsim simulations
      disp_index_vect <- c(rep(NA, times = nsim))
      
      # make the nsim simulations for this value
      for (i in 1:nsim) {
        
        # create the pattern
        mypattern <- rMatClust(kappa = parameter_tab$kappa[index_below], scale = scale_to_test,
                               mu = parameter_tab$mu[index_below], win = my_window)
        
        # remove the points falling outside the window
        if (any(mypattern$x < 0)) {mypattern <- subset(mypattern, !mypattern$x < 0)}
        if (any(mypattern$x > 100)) {mypattern <- subset(mypattern, !mypattern$x > 100)}
        if (any(mypattern$y < 0)) {mypattern <- subset(mypattern, !mypattern$y < 0)}
        if (any(mypattern$y > 100)) {mypattern <- subset(mypattern, !mypattern$y > 100)}
        
        # if mypattern is empty we assign the output values and go to the next iteration
        if (intensity(mypattern) == 0) {
          disp_index_vect[i] <- NaN
          next
        }
        
        # create the matrix with the counts
        mymatrix <- quadratcount(mypattern, nx = 100, ny = 100)
        
        # calculate the dispersion index
        disp_index_vect[i] <- DispersionIndex(pop = mymatrix, nx = 100, ny = 100)
        
        
      } # i
      
      # store the mean dispersion index
      if (my_disp_index == round(mean(disp_index_vect, na.rm = TRUE), digits = 1)) {
        parameter_tab <- rbind(parameter_tab,
                               c(my_intensity,
                                 parameter_tab$kappa[index_below],
                                 parameter_tab$mu[index_below],
                                 scale_to_test,
                                 as.numeric(mean(disp_index_vect, na.rm = TRUE)),
                                 NA))
        parameter_tab$row_status[nrow(parameter_tab)] <- "final_estimate"
        
      } else {
        
        parameter_tab <- rbind(parameter_tab,
                               c(my_intensity,
                                 parameter_tab$kappa[index_below],
                                 parameter_tab$mu[index_below],
                                 scale_to_test,
                                 as.numeric(mean(disp_index_vect, na.rm = TRUE)),
                                 NA))
        parameter_tab$row_status[nrow(parameter_tab)] <- "2nd_step"
      }
      
      # sort parameter_tab and update int_rows
      parameter_tab <- parameter_tab[with(parameter_tab, order(intensity, scale)), ]
      int_rows <- which(parameter_tab$intensity == my_intensity)
      parameter_tab$dispersion_index <- as.numeric(parameter_tab$dispersion_index)
      
    }
    
    
  } # my_disp_index
  
} # my_intensity

# store execution time
end.time <- Sys.time()
time.taken <- difftime(end.time, start.time, units = "mins")
time.taken

# save results
write.csv2(parameter_tab, file = "./data/tidy/Matern_parameters_2021-06-24_kappa_refined.csv", row.names = FALSE)
write.csv2(convergence_fails, file = "./data/tidy/Matern_parameters_2021-06-24_kappa_convergence_fails.csv", row.names = FALSE)


### create the dataframe with all the parameter values coef_tab_total

wanted_intensities = c(1:30)
wanted_disp_index = c(2:75)

# load the table
Matern_parameters_tab <- read.csv2("./data/tidy/Matern_parameters_2021-06-24_kappa_refined.csv")

# rename columns and remove column row_status
colnames(Matern_parameters_tab) <- c("theoretical_intensity", "my_kappa", "my_mu", "my_scale", "expect_disp_ind", "row_status")
Matern_parameters_tab <- subset(Matern_parameters_tab, select = -row_status)

# keep only wanted levels
# make an empty dataframe with the same columns than Matern_parameters_tab
coef_tab_total <- Matern_parameters_tab
coef_tab_total <- coef_tab_total[0, ]

# loop over levels of intensity
for (my_intensity in wanted_intensities) {
  
  # loop over levels of dispersion index
  for (my_disp_index in wanted_disp_index) {
    
    inter_df <- subset(Matern_parameters_tab, Matern_parameters_tab$theoretical_intensity == my_intensity &
                         round(Matern_parameters_tab$expect_disp_ind, digits = 1) == my_disp_index)
    
    # if there is more than one row in inter_df, keep only the row with the closest disp_index
    if (nrow(inter_df) > 1) {
      row_to_keep <- which.min(abs(inter_df$expect_disp_ind - my_disp_index))
      inter_df <- inter_df[row_to_keep, ]
    }
    
    # store kept rows
    coef_tab_total <- rbind(coef_tab_total, inter_df)
    
  } # my_disp_index
} # my_intensity

# sort the dataframe
coef_tab_total <- coef_tab_total[with(coef_tab_total, order(theoretical_intensity, -my_scale)), ]

# add column "pattern" with the name of each pattern
pattern <- c(paste0("int", coef_tab_total$theoretical_intensity, "_disp", round(coef_tab_total$expect_disp_ind, digits = 1)))
coef_tab_total <- cbind(pattern = pattern, coef_tab_total)

# save data
write.csv2(coef_tab_total, file = "./data/tidy/Matern_parameters_coef_tab_total_2021-06-24_kappa.csv", row.names = FALSE)

# # plot scale values
# myplot <- ggplot(coef_tab_total,
#                  aes(x = my_scale, y = expect_disp_ind, color = as.factor(theoretical_intensity))) +
#   geom_line(size = 1) +
#   xlab('Scale parameter') +
#   ylab("Dispersion index") +
#   ggtitle('Dispersion index as a function of the "scale" parameter of the Matern process (mu fixed at 100)') +
#   labs(color = "Intensity") +
#   theme(text = element_text(size = 20))
# myplot


