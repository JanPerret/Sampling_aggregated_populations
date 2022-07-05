
### function to get all pair factors of a given number

get_pair_factors <- function(x) {
  x <- as.integer(x)
  div <- seq_len(abs(x))
  factors <- div[x %% div == 0L]
  return(factors)
}
