
sample_BAS <- function(pop, nx, ny, sample_size, sequence = "halton") {
  
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
  
  
  # check if nx and ny are integer
  if( !nx == as.integer(nx) ) {
    stop("nx has to be an integer")
  }
  if( !ny == as.integer(ny) ) {
    stop("ny has to be an integer")
  }
  
  # check if sample size is less than one
  if( sample_size < 1 ) {
    sample_size <- 1
    warning("Sample size less than one has been reset to 1")
  }
  
  # check sample_size is superior to the number of possible quadrats
  if( sample_size > nx*ny ) {
    sample_size <- nx*ny
    warning("Sample size is bigger than the number of possible quadrats, that's not sampling any more !\nSample size has been reset to the number of quadrats (nx*ny)")
  }
  
  
  # draw a random-start 2-D Halton sequence of length sample_size
  if (sequence == "halton") {
    halt.samp <- halton(n = sample_size*10, dim = 2, start = floor(runif(n = 2, min = 0, max = 100000)))
  } else if (sequence == "R2") {
    halt.samp <- R2seq(n = sample_size*10, start = floor(runif(n = 1, min = 0, max = 100000)))
  } else {
    stop('The type of low discrepancy sequence is not recognized. Only "halton" and "R2" are implemented yet.')
  }
  
  colnames(halt.samp) <- c("x", "y")
  
  # convert from the [0,1] square to a square box covering [bb]
  halt.samp[,1] <- halt.samp[,1] * nx
  halt.samp[,2] <- halt.samp[,2] * ny
  
  # round the grid coordinates to transform them to row and colum indexes
  halt.samp <- ceiling(halt.samp)
  
  # replace all remaining zeros by 1 (zeros are an "errors" due to R's floating point arithmetic in function ceiling)
  halt.samp[halt.samp == 0] <- 1
  
  # delete duplicates
  halt.samp_2 <- unique(halt.samp)
  
  # keep only the needed number of cells
  if (nrow(halt.samp_2) >= sample_size) {
    # select the first sample_size cells
    halt.samp3 <- halt.samp_2[c(1:sample_size), ]
  } else {
    # this condition is here to avoid any risk of error : if even 10*sample_size cells
    # is not enough, we simply keep the duplicates, so there is no error. This should
    # never happen with the values of sample_size we choose for the simulations.
    halt.samp3 <- halt.samp[c(1:sample_size), ]
  }
  
  # get the values from the selected cells
  mysample <- pop_counts[cbind(halt.samp3[,1], halt.samp3[,2])]
  
  # give names of the cells indexes to the count values
  halt.samp_vect <- paste(halt.samp3[,1], halt.samp3[,2], sep = ',')
  names(mysample) <- halt.samp_vect
  
  return(mysample)

}

