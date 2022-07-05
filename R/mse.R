mse <- function(observed, real) {
  
  mean((observed - real)^2)
  
}

# mean square error = average squared difference between the estimated values and the actual value