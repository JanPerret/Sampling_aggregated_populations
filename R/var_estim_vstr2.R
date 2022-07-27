
# function to compute the v_str2 variance estimator. Ref : D'Orazio (2003) eq. (4.1) and (4.2) p. 285
# uses functions from the package {spdep}

var_estim_vstr2 <- function(y, n, N, nx, ny) {
  
  # Check that y is a vector and has n elements:
  if( !is.vector(y) ) cat(" !(is.vector(y)), ERROR: y is not an R vector. ")
  if( length(y)!=n ) cat(" length(y)!=n, ERROR: y is not of length n. ")
  
  f = n/N # sampling fraction
  ybar = mean(y) # sample mean
  
  # sum over j to compute s2
  s2sum = 0.
  for (j in 1:n) {
    s2sum = s2sum + (y[j] - ybar)^2
  }
  s2 = s2sum /(n-1)
  
  
  ### compute Gery's C
  # get the dimensions of the grid formed by the sampling units
  grid_dim <- get_factor_pair_closest_to_square(length(y))
  # grid_dim[1] # nrow
  # grid_dim[2] # ncol
  if(grid_dim[1] != grid_dim[2]) cat("ERROR: vstr2 estimator doesn't work if the sampling units don't form a square grid")
  
  # get the distance between sampling units on the x and y axis
  d <- c(nx/grid_dim[1], ny/grid_dim[2])
  if(d[1] != d[2]) cat("ERROR: vstr2 estimator doesn't work if the simulation window is not a square")
  
  # compute the max distance at which we consider units as neighbours (i.e. max distance is the unit in diagonal)
  max_neighbour_dist <- sqrt(d[1]^2 + d[1]^2)
  
  # get sampling units coordinates
  coord_list <- strsplit(names(y), split = ",")
  coord_y <- as.numeric(unlist(lapply(coord_list, `[[`, 1)))
  coord_x <- as.numeric(unlist(lapply(coord_list, `[[`, 2)))
  df_xy <- data.frame(x = coord_x, y = coord_y)
  
  # compute the list of neighbours (object of type 'nb' from package spdep)
  xy.nb <- spdep::dnearneigh(as.matrix(df_xy), d1 = 0, d2 = max_neighbour_dist, longlat = FALSE)
  
  # add spatial distance based weights to the neighbours 
  my_points <- sp::SpatialPoints(df_xy)
  my_listw <- spdep::nb2listwdist(neighbours = xy.nb, x = my_points, type = "idw", alpha = 1, longlat = TRUE)
  
  # get the unique values in my_listw$weights
  vect_weight <- unique(unlist(my_listw$weights))
  weight_direct_neighbour <- max(vect_weight)
  weight_diagonal_neighbour <- min(vect_weight)
  
  # replace the weight of direct neighbours by 1 and the weight of diagonal neighbours by 1/sqrt(2) ad indicated in Magnussen et al. (2020)
  my_listw$weights <- rapply(my_listw$weights, function(x) ifelse(x == weight_direct_neighbour, 1, x), how = "replace")
  my_listw$weights <- rapply(my_listw$weights, function(x) ifelse(x == weight_diagonal_neighbour, 1/sqrt(2), x), how = "replace")
  
  # compute geary's C
  geary_C <- spdep::geary(x = y, listw = my_listw, n = n, n1 = n-1, S0 = sum(unlist(my_listw$weights)))
  
  # for the case if there are only zeros in the sample
  if(is.nan(geary_C$C)){
    geary_C$C <- 0
  }
  
  ### compute D'Orazio's v_str2 variance estimator
  v_str2 = (1-f) * (s2/n) * geary_C$C
  
  return(v_str2)
}

