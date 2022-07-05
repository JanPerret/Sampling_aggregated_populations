
### function to get the factor pair closest to a square
# NB : if x is a prime, the function returns c(NA, NA)

get_factor_pair_closest_to_square <- function(x) {
  
  # if the squareroot of x is integer its simple
  if(sqrt(x) == as.integer(sqrt(x))) {
    
    closest_pair <- c(sqrt(x), sqrt(x))
    
  } else {
    
    # get all pair factors
    factors_vect <- get_pair_factors(x)
    
    # function for the position argument
    p <- function(f, b) function(a) f(a, b)
    
    # find position of the first value in the vector superior to sqrt(x)
    pos <- Position(p(`>`, sqrt(x)), factors_vect)
    
    # get the factor pair closest to a square
    closest_pair <- c(factors_vect[pos-1], factors_vect[pos])
  }
  
  return(closest_pair)
}
