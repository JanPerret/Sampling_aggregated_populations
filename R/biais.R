biais <- function(observed, real) {
  
  mean(observed) - real
  
}

# bias = average difference between sample and population
# relative bias = bias divided by the real value
# percent bias = relative bias * 100