
# function to compute the v_8 variance estimator. Ref : Wolter (2007) eq. (8.2.9) p. 302
# code adapted from the supplementary material of McGarvey et al. (2016)

var_estim_v8 <- function(y, n, N) {
  
  if(n != length(y)) print("ERROR: n != length(y)")
  if(any(is.na(y))) print("ERROR: NA's in y")
  
  f = n/N # sampling fraction
  ybar = mean(y) # sample mean
  
  # sum over j to compute s2
  s2sum = 0.
  for (j in 1:n) {
    s2sum = s2sum + (y[j] - ybar)^2
  } # end sum over j to compute s2
  s2 = s2sum /(n-1)
  
  if (s2 == 0) { # if the sample contains only zeros the estimated variance is set to zero
    
    v_8 = (1-f) * (s2/n)
    
  } else {
    
    # sum over j to compute rhohatp
    rhohatpsum = 0.
    for (j in 2:n) {
      rhohatpsum = rhohatpsum + (y[j] - ybar) * (y[j-1] - ybar)
    } # end sum over j to compute rhohatp
    rhohatp = rhohatpsum / ((n-1)*s2)
    
    if (rhohatp > 0) {v_8 = (1-f) * (s2/n) * (1 + (2/log(rhohatp)) + 2/((1/rhohatp) - 1))
    #  if (rhohatp <= 0) {
    #    #print("PROBLEM:  rhohatp <= 0 ")
    #    v_8 = (1-f) * (s2/n)
    #  }
    } else {v_8 = (1-f) * (s2/n)}
    
  }
  
  
  return(v_8)
}

