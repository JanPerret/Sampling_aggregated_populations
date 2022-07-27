
# function to compute the v_nbh variance estimator. Ref : Steven and Olsen (2003) last equation on page 601
# code adapted from the {spsurvey} package : https://cran.r-project.org/web/packages/spsurvey/index.html

var_estim_vnbh <- function(y, n, N) {
  
  if(n != length(y)) print("ERROR: n != length(y)")
  if(any(is.na(y))) print("ERROR: NA's in y")
  
  f = n/N # sampling fraction
  ybar = mean(y) # sample mean
  
  # approximate inclusion probability
  inclusion_prob <- rep(n/N, times = n)
  
  # other stuff we need
  wgt <- 1/inclusion_prob
  
  # calculate the weighted residuals values
  rv_mean <- wgt * (y - ybar)
  
  # get sampling units coordinates
  coord_list <- strsplit(names(y), split = ",")
  coord_y <- as.numeric(unlist(lapply(coord_list, `[[`, 1)))
  coord_x <- as.numeric(unlist(lapply(coord_list, `[[`, 2)))
  
  # compute the index values of neighboring points and associated weights required by the local mean variance estimator.
  weight_1st <- spsurvey::localmean_weight(x = coord_x, y = coord_y, prb = inclusion_prob)
  
  # compute the local mean estimator for population mean
  v_nbh <- spsurvey::localmean_var(z = rv_mean, weight_1st = weight_1st) / (N^2)
  # the firts term is the variance estimate for the population total, so we divide by N^2 to get the variance estimate for the population mean
  
  
  return(v_nbh)
}
