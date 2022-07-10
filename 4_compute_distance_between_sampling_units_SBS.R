#######################################################
# 
# Sampling simulation of spatially aggregated populations
# 
# compute_distance_between_sampling_units_SBS.R
# 
# jan.perret@gmail.com
#######################################################

###
### script objective : compute the mean distance between the sampling units for Spatially Balanced Sampling (SBS, also sometimes called BAS in the scripts)
###

### clean workspace
rm(list = ls())

### load packages
source("./R/imports.R")

### load functions
files.sources <- list.files("./R", full.names = TRUE)
invisible(sapply(files.sources, source))


### get mean distance to nearest neighbours for Halton sequence ####
nsim = 5000
nx = 100
ny = 100
sample_size_vect <- c(9, 16, 25, 36, 49, 64, 81, 100, 121, 144, 169, 196, 225,
                      256, 289, 324, 361, 400, 441, 484, 529, 576, 625)

# make table to save results
dist_result_tab <- data.frame(sample_size = sample_size_vect,
                              mean_dist_1rst_neighbour = rep(NA, times = length(sample_size_vect)),
                              mean_dist_2nd_neighbour = rep(NA, times = length(sample_size_vect)),
                              mean_dist_3rd_neighbour = rep(NA, times = length(sample_size_vect)),
                              mean_dist_4th_neighbour = rep(NA, times = length(sample_size_vect)),
                              mean_dist_5th_neighbour = rep(NA, times = length(sample_size_vect)),
                              mean_dist_6th_neighbour = rep(NA, times = length(sample_size_vect)),
                              mean_dist_7th_neighbour = rep(NA, times = length(sample_size_vect)),
                              mean_dist_8th_neighbour = rep(NA, times = length(sample_size_vect)),
                              mean_dist_9th_neighbour = rep(NA, times = length(sample_size_vect)),
                              mean_dist_10th_neighbour = rep(NA, times = length(sample_size_vect)),
                              mean_dist_11th_neighbour = rep(NA, times = length(sample_size_vect)),
                              mean_dist_12th_neighbour = rep(NA, times = length(sample_size_vect)),
                              mean_dist_13th_neighbour = rep(NA, times = length(sample_size_vect)),
                              mean_dist_14th_neighbour = rep(NA, times = length(sample_size_vect)))

# loop over sample sizes
for (sample_size in sample_size_vect) {
  
  # print progress
  cat(paste("n =", sample_size, "started at", Sys.time(), "\n"))
  
  # initialize empty vector to store mean distance to nearest units for each simulation
  dist_1rst_neighbour_vect <- c(rep(NA, times = nsim))
  dist_2nd_neighbour_vect <- c(rep(NA, times = nsim))
  dist_3rd_neighbour_vect <- c(rep(NA, times = nsim))
  dist_4th_neighbour_vect <- c(rep(NA, times = nsim))
  dist_5th_neighbour_vect <- c(rep(NA, times = nsim))
  dist_6th_neighbour_vect <- c(rep(NA, times = nsim))
  dist_7th_neighbour_vect <- c(rep(NA, times = nsim))
  dist_8th_neighbour_vect <- c(rep(NA, times = nsim))
  dist_9th_neighbour_vect <- c(rep(NA, times = nsim))
  dist_10th_neighbour_vect <- c(rep(NA, times = nsim))
  dist_11th_neighbour_vect <- c(rep(NA, times = nsim))
  dist_12th_neighbour_vect <- c(rep(NA, times = nsim))
  dist_13th_neighbour_vect <- c(rep(NA, times = nsim))
  dist_14th_neighbour_vect <- c(rep(NA, times = nsim))
  
  
  for (i in 1:nsim) {
    
    # draw an halton sequence
    halt.samp <- halton(n = sample_size*10, dim = 2, start = floor(runif(n = 2, min = 0, max = 100000)))
    colnames(halt.samp) <- c("x", "y")
    
    # add a column with rownum
    halt.samp <- cbind(halt.samp, rownum = seq(1:(sample_size*10)))
    
    # convert from the [0, 1] square to a square box covering [bb]
    halt.samp[,1] <- halt.samp[,1] * nx
    halt.samp[,2] <- halt.samp[,2] * ny
    
    # round the grid coordinates to transform them to row and column indexes
    halt.samp_rounded <- ceiling(halt.samp)
    
    # replace all remaining zeros by 1 (zeros are an "errors" due to R's floating point arithmetic in function ceiling)
    halt.samp_rounded[halt.samp_rounded == 0] <- 1
    
    # delete duplicates
    halt.samp_rounded <- subset(halt.samp_rounded, !duplicated(halt.samp_rounded[, c(1, 2)]))
    
    # keep only the needed number of cells
    halt.samp_rounded <- halt.samp_rounded[c(1:sample_size), ]
    
    # filter halt.samp in the same manner than halt.samp_rounded
    halt.samp <- subset(halt.samp, halt.samp[,3] %in% halt.samp_rounded[,3])
    
    # convert halt.samp to a ppp object
    halt.samp_ppp <- ppp(x = halt.samp[,1], y = halt.samp[,2], window = owin(c(0, 100), c(0, 100)))
    
    # get distance to neirest neighbour
    dist_1rst_neighbour_vect[i] <- mean(nndist(halt.samp_ppp, k = 1))
    dist_2nd_neighbour_vect[i] <- mean(nndist(halt.samp_ppp, k = 2))
    dist_3rd_neighbour_vect[i] <- mean(nndist(halt.samp_ppp, k = 3))
    dist_4th_neighbour_vect[i] <- mean(nndist(halt.samp_ppp, k = 4))
    dist_5th_neighbour_vect[i] <- mean(nndist(halt.samp_ppp, k = 5))
    dist_6th_neighbour_vect[i] <- mean(nndist(halt.samp_ppp, k = 6))
    dist_7th_neighbour_vect[i] <- mean(nndist(halt.samp_ppp, k = 7))
    dist_8th_neighbour_vect[i] <- mean(nndist(halt.samp_ppp, k = 8))
    dist_9th_neighbour_vect[i] <- mean(nndist(halt.samp_ppp, k = 9))
    dist_10th_neighbour_vect[i] <- mean(nndist(halt.samp_ppp, k = 10))
    dist_11th_neighbour_vect[i] <- mean(nndist(halt.samp_ppp, k = 11))
    dist_12th_neighbour_vect[i] <- mean(nndist(halt.samp_ppp, k = 12))
    dist_13th_neighbour_vect[i] <- mean(nndist(halt.samp_ppp, k = 13))
    dist_14th_neighbour_vect[i] <- mean(nndist(halt.samp_ppp, k = 14))
    
    
  } # i loop
  
  # store result for current sample_size
  dist_result_tab$mean_dist_1rst_neighbour[which(dist_result_tab$sample_size == sample_size)] <- mean(dist_1rst_neighbour_vect)
  dist_result_tab$mean_dist_2nd_neighbour[which(dist_result_tab$sample_size == sample_size)] <- mean(dist_2nd_neighbour_vect)
  dist_result_tab$mean_dist_3rd_neighbour[which(dist_result_tab$sample_size == sample_size)] <- mean(dist_3rd_neighbour_vect)
  dist_result_tab$mean_dist_4th_neighbour[which(dist_result_tab$sample_size == sample_size)] <- mean(dist_4th_neighbour_vect)
  dist_result_tab$mean_dist_5th_neighbour[which(dist_result_tab$sample_size == sample_size)] <- mean(dist_5th_neighbour_vect)
  dist_result_tab$mean_dist_6th_neighbour[which(dist_result_tab$sample_size == sample_size)] <- mean(dist_6th_neighbour_vect)
  dist_result_tab$mean_dist_7th_neighbour[which(dist_result_tab$sample_size == sample_size)] <- mean(dist_7th_neighbour_vect)
  dist_result_tab$mean_dist_8th_neighbour[which(dist_result_tab$sample_size == sample_size)] <- mean(dist_8th_neighbour_vect)
  dist_result_tab$mean_dist_9th_neighbour[which(dist_result_tab$sample_size == sample_size)] <- mean(dist_9th_neighbour_vect)
  dist_result_tab$mean_dist_10th_neighbour[which(dist_result_tab$sample_size == sample_size)] <- mean(dist_10th_neighbour_vect)
  dist_result_tab$mean_dist_11th_neighbour[which(dist_result_tab$sample_size == sample_size)] <- mean(dist_11th_neighbour_vect)
  dist_result_tab$mean_dist_12th_neighbour[which(dist_result_tab$sample_size == sample_size)] <- mean(dist_12th_neighbour_vect)
  dist_result_tab$mean_dist_13th_neighbour[which(dist_result_tab$sample_size == sample_size)] <- mean(dist_13th_neighbour_vect)
  dist_result_tab$mean_dist_14th_neighbour[which(dist_result_tab$sample_size == sample_size)] <- mean(dist_14th_neighbour_vect)
  
  
} # sample_size loop

# save results
write.csv2(dist_result_tab, file = "./data/tidy/BAS_mean_distance_nearest_units_14_neighbours.csv", row.names = FALSE)


