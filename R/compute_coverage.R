
# function to compute the effective coverage ratio for a given variance estimator

compute_coverage <- function(Ybars, ybars, estim_vars, n, N) {
  
  f <- n/N # sampling fraction
  
  CI_lower <- ybars - 1.96*sqrt(estim_vars)
  CI_upper <- ybars + 1.96*sqrt(estim_vars)
  
  my_coverage <- sum((Ybars > CI_lower) & (Ybars < CI_upper)) / length(Ybars)
  
  return(my_coverage)
}