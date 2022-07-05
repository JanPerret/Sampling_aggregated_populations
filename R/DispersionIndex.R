

DispersionIndex <- function(pop, nx, ny){

  if(class(pop)[1] == "ppp") {
    
    # check if the study area is a square
    win_type <- pop$window$type
    win_dim <- c(pop$window$xrange, pop$window$yrange)
    dim_boolean <- any(win_dim[c(1,2)] != win_dim[c(3,4)])
    if( win_type != "rectangle" | dim_boolean ){
      stop("This function only works for square study areas")
    }
    
    # get the matrix with the counts
    pop_counts <- quadratcount(pop, nx = nx, ny = ny)
    
  } else if(class(pop)[1] == "quadratcount" | class(pop)[1] == "data.frame" | class(pop)[1] == "matrix") {
    
    # check if the matrix is square
    if(nrow(pop) != ncol(pop)) {
      stop("This function only works for square study areas")
    }
    
    # check if the matrix is of dimensions nx*ny
    if(ncol(pop) != nx | nrow(pop) != ny){
      stop("If 'pop' is a of type 'quadratcount' it has to be of same dimensions as nx and ny.")
    }
    
    # assign the matrix into pop_counts
    pop_counts <- pop
    
  } else {
    
    stop("Object given to argument 'pop' has to be of type 'ppp' or 'quadratcount'.")
  }
  
  
  # compute dispersion index (see Baddeley et al. 2016 page 200-201)
  disp_ind <- var(c(pop_counts)) / mean(pop_counts)
  
  return(disp_ind)
}

