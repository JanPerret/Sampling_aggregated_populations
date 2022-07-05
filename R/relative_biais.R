relative_biais <- function(observed, real) {
  
  (mean(observed) - real)/real
  
}

# bias = average difference between sample and population
# relative bias = bias divided by the true value
# percent bias = relative bias * 100