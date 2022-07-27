
# function to compute the v_w variance estimator. Ref : D'Orazio (2003) table 2 & eq. (4.3) p. 284 & 285
# code adapted from the supplementary material of McGarvey et al. (2016)

var_estim_vw <- function(y, n, N) {
  
  # Check that y is a vector and has n elements:
  if( !is.vector(y) ) cat("!(is.vector(y)), ERROR: y is not an R vector.")
  if( length(y) != n ) cat("length(y) != n, ERROR: y is not of length n.")
  
  f = n/N # sampling fraction
  ybar = mean(y) # sample mean
  
  # sum over j to compute s2
  s2sum = 0.
  for (j in 1:n) {  
    s2sum = s2sum + (y[j] - ybar)^2  
  } # end sum over j to compute s2
  s2 = s2sum /(n-1)
  
  if (s2 == 0) { # if the sample contains only zeros the estimated variance is set to zero
    
    v_w = (1-f) * (s2/n)
    
  } else {
    
    ### compute Moran's I
    # first we compute the weight matrix
    # get the dimensions of the grid formed by the sampling units
    grid_dim <- get_factor_pair_closest_to_square(length(y))
    # grid_dim[1] # nrow 
    # grid_dim[2] # ncol
    
    # create matrix of dimension n*n filled with zeros
    weight_matrix <- matrix(data = 0, nrow = length(y), ncol = length(y), byrow = TRUE)
    
    
    # set the cells of the weight matrix to 1 for the neighbouring sampling units of every sampling unit 
    for (i in 1:length(y)) {
      
      # column neighbours
      if (i %% grid_dim[1] != 1) {
        weight_matrix[i, i-1] <- 1
        weight_matrix[i-1, i] <- 1
      }
      
      if (i %% grid_dim[1] != 0) {
        weight_matrix[i, i+1] <- 1
        weight_matrix[i+1, i] <- 1
      }
      
      # row neighbours
      if (i > grid_dim[1]) {
        weight_matrix[i, i-grid_dim[1]] <- 1
        weight_matrix[i-grid_dim[1], i] <- 1
      }
      
      if (i <= (grid_dim[2]-1)*grid_dim[1]) {
        weight_matrix[i, i+grid_dim[1]] <- 1
        weight_matrix[i+grid_dim[1], i] <- 1
      }
      
    } # i
    
    # compute Moran's I using function Moran.I() from package {ape}
    IMoran <- ape::Moran.I(y, weight_matrix)
    IMoran <- IMoran$observed
    
    ### compute D'Orazio v_w variance estimator
    if (IMoran > 0) {
      v_w = (1-f) * (s2/n) * (1 + (2/log(IMoran)) + 2/((1/IMoran) - 1))
    } else {
      v_w = (1-f) * (s2/n)
    }
    
  }
  
  return(v_w)
}

