
#### function to compute our new estimator based on the hypothesis of random sampling from a ZIP distribution ####
var_estim_ZIP <- function(y, n, N) {
  
  if(n != length(y)) print("ERROR: n != length(y)")
  if(any(is.na(y))) print("ERROR: NA's in y")
  
  f = n/N # sampling fraction
  
  
  # if y are only zeros, variance is zero
  if (all(y == 0)) {
    
    v_ZIP = 0
    
  # if there are no zeros in y we consider y comes from a poisson distribution,
  # so the sample variance is an unbiased estimator of the population variance
  } else if (!any(y == 0)) {
    
    v_ZIP = (var(y) / n) * (1 - f)
    
  } else {
    
    # convert the vector with the counts to a dataframe
    df_sample_counts <- as.data.frame(y)
    
    # fit ZIP regression
    fm_zip <- pscl::zeroinfl(y ~ 1 , data = df_sample_counts)
    
    # backtransform estimates
    x <- unname(coef(fm_zip)[2])
    pi <- exp(x)/(1 + exp(x)) # inverse logit
    
    lambda <- unname(coef(fm_zip)[1])
    lambda <- exp(lambda)
    
    # compute population mean and variance
    estim_pop_mean <- (1 - pi) * lambda
    estim_pop_var <- lambda * (1 - pi) * (1 + pi * lambda)
    
    # compute sampling variance estimate
    v_ZIP = (estim_pop_var / n) * (1 - f)
    
  }
  return(v_ZIP)
}
