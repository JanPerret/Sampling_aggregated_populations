rmse <- function(observed, real) {
  
  sqrt(mean((observed - real)^2))
  
}

# square root of the MSE