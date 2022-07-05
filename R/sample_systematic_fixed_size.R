
### function to draw a systematic sample of the exact specified sample size
# NB : if sample_size is a prime, the function will return an error.
# To change this you just have to replace the second version of get_pair_factors()
# by the first, so for primes it will return 1 and the specified prime instead of
# only 1.
# NB2 : I changed the default to the first version of get_pair_factors(), so if 
# sample_size if a prime, the function will arrange the plots as a line, and if
# sample_size is superior to the number of cells, it will return less than sample_size values.

sample_systematic_fixed_size <- function(pop, nx, ny, sample_size) {
  
  # check if sample_size is superior to the number of possible quadrats
  if( sample_size > nx*ny ) {
    sample_size <- nx*ny
    warning("Sample size is bigger than the number of possible quadrats, that's not sampling any more !
            Sample size has been reset to the number of quadrats (nx*ny).")
  }
  
  
  if(class(pop)[1] == "ppp") {
    
    # check if the study area is a square
    win_type <- pop$window$type
    win_dim <- c(pop$window$xrange, pop$window$yrange)
    dim_boolean <- any(win_dim[c(1, 2)] != win_dim[c(3, 4)])
    if( win_type != "rectangle" | dim_boolean ){
      stop("This function only works for square study areas")
    }
    
    # get the matrix with the counts
    pop_counts <- quadratcount(pop, nx = nx, ny = ny)
    
  } else if(class(pop)[1] == "quadratcount" | class(pop)[1] == "data.frame") {
    
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
  
  
  # get the dimensions of the grid closest to a square
  grid_dim <- get_factor_pair_closest_to_square(sample_size)
  
  # compute spacing on the x and y axis
  delta <- c(nx/grid_dim[1], ny/grid_dim[2])
  
  # draw the random start point inside a rectangle of dimension delta[1]*delta[2]
  # with the (0,0) point as corner
  m.x <- runif(n = 1, min = 0, max = delta[1] )
  m.y <- runif(n = 1, min = 0, max = delta[2] )
  
  # get the grid extent
  d <- max(nx, ny)
  
  # get the coordinates of the "rows" and the "columns" of the grid
  seq.x <- seq( 0, d, by = delta[1] ) + m.x
  seq.y <- seq( 0, d, by = delta[2] ) + m.y
  
  # remove points that fall outside the study area
  seq.x <- seq.x[seq.x < nx]
  seq.y <- seq.y[seq.y < ny]
  
  # make the grid
  grd <- expand.grid( x = seq.x, y = seq.y )
  
  # round the grid coordinates to transform them to row and column indexes
  grd <- ceiling(grd)
  
  # check for duplicates
  if(any(duplicated(grd))) {
    stop("Some quadrats where drawn more than once in the sample, probably because
    the chosen sample_size leeds to a spatial arrangement of the quadrats that
         exeeds the available length or width of the study area.")
  }
  
  
  # get the values from the selected cells
  mysample <- pop_counts[cbind(grd$x, grd$y)]
  
  # give names of the cells indexes to the count values
  grd_vect <- paste(grd$x, grd$y, sep = ',')
  names(mysample) <- grd_vect
  
  return(mysample)
}
