#######################################################
# 
# Sampling simulation of spatially aggregated populations
# 
# 1_estimate_point_process_parameters_SSI.R
# 
# jan.perret@gmail.com
#######################################################


### clean workspace
rm(list = ls())

### load packages
source("./R/imports.R")

### load functions
files.sources <- list.files("./R", full.names = TRUE)
invisible(sapply(files.sources, source))

### general settings
my_window <- owin(c(0, 50), c(0, 50)) # simulation window
wanted_intensities = c(1, 5, 10, 20)
wanted_disp_index = c(0.9, 0.8, 0.7, 0.6, 0.5, 0.4)
nsim = 200


# as all radius values are not possible for every intensity, we set adapted limits per intensity level
radius_vect_int1 = seq(from = 0.05, to = 0.75, by = 0.02)
radius_vect_int5 = seq(from = 0.05, to = 0.27, by = 0.01)
radius_vect_int10 = seq(from = 0.01, to = 0.19, by = 0.01)
radius_vect_int20 = seq(from = 0.01, to = 0.12, by = 0.01)

# create simulation input table
input_table <- data.frame(cbind(
  intensity = c(rep(wanted_intensities[1], times = length(radius_vect_int1)),
                rep(wanted_intensities[2], times = length(radius_vect_int5)),
                rep(wanted_intensities[3], times = length(radius_vect_int10)),
                rep(wanted_intensities[4], times = length(radius_vect_int20))),
  my_r = c(radius_vect_int1, radius_vect_int5, radius_vect_int10, radius_vect_int20),
  dispersion_index = NA,
  proportion_NA = NA
))



### first function to estimate the dispersion index for the values of "r" in input_table
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
    mypattern <- rSSI(r = input_table$my_r[j],
                      n = input_table$intensity[j] * my_window$xrange[2] * my_window$yrange[2], # wanted intensity * area of the simulation window
                      win = my_window, giveup = 2000)

    # just in case there is the same bug than in rMatClust() with sometimes points falling outside
    # of the simulation window, which creates an error in the function quadratcount()
    if (any(mypattern$x < 0)) {mypattern <- subset(mypattern, !mypattern$x < 0)}
    if (any(mypattern$x > 50)) {mypattern <- subset(mypattern, !mypattern$x > 50)}
    if (any(mypattern$y < 0)) {mypattern <- subset(mypattern, !mypattern$y < 0)}
    if (any(mypattern$y > 50)) {mypattern <- subset(mypattern, !mypattern$y > 50)}

    # if mypattern is empty we assign the output values and go to the next iteration
    if (intensity(mypattern) == 0) {
      disp_index_vect[i] <- NaN
      next
    }

    # if the requested number of points is not in the pattern, don't use it
    if (mypattern$n == input_table$intensity[j] * my_window$xrange[2] * my_window$yrange[2]) {

      # create the matrix with the counts
      mymatrix <- quadratcount(mypattern, nx = 50, ny = 50)

      # calculate the dispersion index
      disp_index_vect[i] <- DispersionIndex(pop = mymatrix, nx = 50, ny = 50)

    } else {

      disp_index_vect[i] <- NA

    }

  } # i

  # store the mean dispersion index
  input_table$dispersion_index[j] <- mean(disp_index_vect, na.rm = TRUE)
  input_table$proportion_NA[j] <- sum(is.na(disp_index_vect)) / nsim

} # j

# store execution time
end.time <- Sys.time()
time.taken <- difftime(end.time, start.time, units = "mins")
time.taken

write.csv2(input_table, file = "./data/tidy/SSI_parameters_2021-07-26.csv", row.names = FALSE)



### second function to refine the estimates
# load the table 
SSI_parameters_tab <- read.csv2("./data/tidy/SSI_parameters_2021-07-26.csv")

# add the column row_status (1rst_step / 2nd_step / final_value)
SSI_parameters_tab <- cbind(SSI_parameters_tab, row_status = "1rst_step", stringsAsFactors = FALSE)

### check if the column dispersion_index is monotonic and if not, remove the corresponding rows 
SSI_parameters_tab <- SSI_parameters_tab[with(SSI_parameters_tab, order(intensity, my_r)), ]

# loop until all values remaining in the table are monotonic
while (!all(SSI_parameters_tab[with(SSI_parameters_tab, order(intensity, my_r)), ] == 
            SSI_parameters_tab[with(SSI_parameters_tab, order(intensity, -dispersion_index)), ])) {
  
  # vector of rows to delete
  rows_to_delete <- c()
  
  for (i in 2:nrow(SSI_parameters_tab)) {
    
    if (SSI_parameters_tab$intensity[i-1] != SSI_parameters_tab$intensity[i]) {
      next
    }
    
    if (SSI_parameters_tab$dispersion_index[i-1] < SSI_parameters_tab$dispersion_index[i]) {
      cat(paste("Non-monotonic value at line", i))
      cat("\n")
      
      # store row indexes to delete
      rows_to_delete <- c(rows_to_delete, i)
      
    }
    
  } # i
  
  # remove the non monotonic values
  SSI_parameters_tab <- SSI_parameters_tab[-rows_to_delete, ]
  
} # end of the while loop



### check if all wanted values of dispersion index are between extreme values for every intensity level

parameter_tab = SSI_parameters_tab
parameter_tab$dispersion_index <- as.numeric(parameter_tab$dispersion_index)

# loop over levels of intensity
for (my_intensity in wanted_intensities) {
  
  # table with parameter values for current level of intensity
  current_intensity_tab <- subset(SSI_parameters_tab, SSI_parameters_tab$intensity == my_intensity)
  
  if(min(wanted_disp_index) < min(current_intensity_tab$dispersion_index)) {
    cat(paste("For intensity =", my_intensity, ", the lowest wanted value of dispersion index is outside the limits"))
    cat("\n")
  }
  
  if(max(wanted_disp_index) > max(current_intensity_tab$dispersion_index)) {
    cat(paste("For intensity =", my_intensity, ", the highest wanted value of dispersion index is outside the limits"))
    cat("\n")
  }
  
} # my_intensity



### Start the refining of r estimates

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
      index_below <- Position(function(x) x < my_disp_index, parameter_tab[int_rows, ]$dispersion_index) + min(int_rows) - 1
      index_above <- index_below - 1
      
      ### refine the value of my_r by iterative dichotomy
      below_r <- parameter_tab$my_r[index_below]
      above_r <- parameter_tab$my_r[index_above]
      
      # get the middle value between the 2 r values
      r_to_test <- (below_r + above_r) / 2
      
      # create an empty vector to store the values of dispersion index of the nsim simulations
      disp_index_vect <- c(rep(NA, times = nsim))
      
      # make the nsim simulations for this value
      for (i in 1:nsim) {
        
        # create the pattern
        mypattern <- rSSI(r = r_to_test,
                          n = parameter_tab$intensity[index_below] * my_window$xrange[2] * my_window$yrange[2], # wanted intensity * area of the simulation window
                          win = my_window, giveup = 2000)
        
        # remove the points falling outside the window
        if (any(mypattern$x < 0)) {mypattern <- subset(mypattern, !mypattern$x < 0)}
        if (any(mypattern$x > 50)) {mypattern <- subset(mypattern, !mypattern$x > 50)}
        if (any(mypattern$y < 0)) {mypattern <- subset(mypattern, !mypattern$y < 0)}
        if (any(mypattern$y > 50)) {mypattern <- subset(mypattern, !mypattern$y > 50)}
        
        # if mypattern is empty we assign the output values and go to the next iteration
        if (intensity(mypattern) == 0) {
          disp_index_vect[i] <- NaN
          next
        }
        
        # if the requested number of points is not in the pattern, don't use it
        if (mypattern$n == input_table$intensity[j] * my_window$xrange[2] * my_window$yrange[2]) {
          
          # create the matrix with the counts
          mymatrix <- quadratcount(mypattern, nx = 50, ny = 50)
          
          # calculate the dispersion index
          disp_index_vect[i] <- DispersionIndex(pop = mymatrix, nx = 50, ny = 50)
          
        } else {
          
          disp_index_vect[i] <- NA
          
        }
        
      } # i
      
      # store the mean dispersion index
      if (my_disp_index == round(mean(disp_index_vect, na.rm = TRUE), digits = 1)) {
        parameter_tab <- rbind(parameter_tab,
                               c(my_intensity,
                                 r_to_test,
                                 as.numeric(mean(disp_index_vect, na.rm = TRUE)),
                                 sum(is.na(disp_index_vect)) / nsim,
                                 NA))
        parameter_tab$row_status[nrow(parameter_tab)] <- "final_estimate"
        
      } else {
        
        parameter_tab <- rbind(parameter_tab,
                               c(my_intensity,
                                 parameter_tab$kappa[index_below],
                                 parameter_tab$mu[index_below],
                                 r_to_test,
                                 as.numeric(mean(disp_index_vect, na.rm = TRUE)),
                                 NA))
        parameter_tab$row_status[nrow(parameter_tab)] <- "2nd_step"
      }
      
      # sort parameter_tab and update int_rows
      parameter_tab <- parameter_tab[with(parameter_tab, order(intensity, my_r)), ]
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
write.csv2(parameter_tab, file = "./data/tidy/SSI_parameters_2021-07-26_refined.csv", row.names = FALSE)
write.csv2(convergence_fails, file = "./data/tidy/SSI_parameters_2021-07-26_convergence_fails.csv", row.names = FALSE)

# # plot scale values
# myplot <- ggplot(parameter_tab,
#                  aes(x = my_r, y = dispersion_index, color = as.factor(intensity))) +
#   geom_line(size = 1) +
#   xlab('r parameter') +
#   ylab("Dispersion index") +
#   ggtitle('Dispersion index as a function of the "r" parameter of the SSI process') +
#   labs(color = "Intensity") +
#   theme(text = element_text(size = 20))
# myplot



### create the dataframe with all the parameter values coef_tab_total

wanted_intensities = c(1, 5, 10, 20)
wanted_disp_index = c(0.9, 0.8, 0.7, 0.6, 0.5, 0.4)

# load the table 
SSI_parameters_tab <- read.csv2("./data/tidy/SSI_parameters_2021-07-26_refined.csv")

# rename columns and remove column row_status
colnames(SSI_parameters_tab) <- c("theoretical_intensity", "my_r", "expect_disp_ind", "proportion_NA", "row_status")
SSI_parameters_tab <- subset(SSI_parameters_tab, select = -c(proportion_NA, row_status))

# keep only wanted levels
# make an empty dataframe with the same columns than SSI_parameters_tab
coef_tab_total <- SSI_parameters_tab
coef_tab_total <- coef_tab_total[0, ]

# loop over levels of intensity
for (my_intensity in wanted_intensities) {
  
  # loop over levels of dispersion index
  for (my_disp_index in wanted_disp_index) {
    
    inter_df <- subset(SSI_parameters_tab, SSI_parameters_tab$theoretical_intensity == my_intensity &
                         round(SSI_parameters_tab$expect_disp_ind, digits = 1) == my_disp_index)
    
    # if there is more than one row in inter_df, keep only the row with the closest disp_index
    if (nrow(inter_df) > 1) {
      row_to_keep <- which.min(abs(inter_df$expect_disp_ind - my_disp_index))
      inter_df <- inter_df[row_to_keep, ]
    }
    
    # store keept rows
    coef_tab_total <- rbind(coef_tab_total, inter_df)
    
  } # my_disp_index
} # my_intensity

# sort the dataframe
coef_tab_total <- coef_tab_total[with(coef_tab_total, order(theoretical_intensity, -my_r)), ]

# add column "pattern" with the name of each pattern
pattern <- c(paste0("int", coef_tab_total$theoretical_intensity, "_disp", round(coef_tab_total$expect_disp_ind, digits = 1)))
coef_tab_total <- cbind(pattern = pattern, coef_tab_total)

# save data
write.csv2(coef_tab_total, file = "./data/tidy/SSI_parameters_coef_tab_total_2021-07-26.csv", row.names = FALSE)

# # plot scale values
# myplot <- ggplot(coef_tab_total,
#                  aes(x = my_r, y = expect_disp_ind, color = as.factor(theoretical_intensity))) +
#   geom_line(size = 1) +
#   xlab('r parameter') +
#   ylab("Dispersion index") +
#   ggtitle('Dispersion index as a function of the "r" parameter of the SSI process') +
#   labs(color = "Intensity") +
#   theme(text = element_text(size = 15))
# myplot




