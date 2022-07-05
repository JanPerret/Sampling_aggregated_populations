
sample_random <- function(pop, nx, ny, sample_size) {
  
  if(class(pop)[1] == "ppp") {
    
    pop_counts <- quadratcount(pop, nx = nx, ny = ny)
    
  } else if(class(pop)[1] == "quadratcount" | class(pop)[1] == "data.frame") {
    
    # check if the matrix is of dimensions nx*ny
    if(ncol(pop) != nx | nrow(pop) != ny){
      stop("If 'pop' is a of type 'quadratcount' it has to be of same dimensions as nx and ny.")
    }
    
    pop_counts <- pop
    
  } else {
    
    stop("Object given to argument 'pop' has to be of type 'ppp' or 'quadratcount'.")
    
  }
  
  coord <- expand.grid(x = seq_len(nrow(pop_counts)), y = seq_len(ncol(pop_counts)))
  
  mysample_coord <- coord[sample(nrow(coord), replace = FALSE, size = sample_size), ]
  
  mysample <- pop_counts[cbind(mysample_coord$x, mysample_coord$y)]
  
  # give names of the cells indexes to the count values
  index_vect <- paste(mysample_coord$x, mysample_coord$y, sep = ',')
  names(mysample) <- index_vect
  
  return(mysample)
  
}
