
# function to compute the basic SRS variance estimator. Ref : Cochran (1977) page 26 equation 2.20

var_estim_vtextbook <- function(y, n, N) {
  
  if(n != length(y)) print("ERROR: n != length(y)")
  if(any(is.na(y))) print("ERROR: NA's in y")
  
  f = n/N # sampling fraction
  ybar = mean(y) # sample mean
  sqdifs = (y-ybar)^2
  sumsqdifs = sum(sqdifs)
  s2 = sumsqdifs/(n-1)
  v_textbook = (s2/n)*(1-f)
  
  return(v_textbook)
}