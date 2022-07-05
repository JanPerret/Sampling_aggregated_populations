
### Compute a new low discrepancy sequence (additive recurrence R-sequence)
### SOURCE : http://extremelearning.com.au/unreasonable-effectiveness-of-quasirandom-sequences/#GeneralizingGoldenRatio

R2seq <- function (n, start = 0) {
  
  g = 1.32471795724474602596 # approximate value of the plastic constant with 20 digits
  # g = ((9+sqrt(69))/18)^(1/3) + ((9-sqrt(69))/18)^(1/3) # exact value of the plastic constant
  
  a1 = 1.0/g
  a2 = 1.0/(g*g)
  x = c()
  y = c()
  count = 1
  
  for (i in start:(start + n - 1)) {
    x[count] = (0 + a1*i) %% 1
    y[count] = (0 + a2*i) %% 1
    count = count + 1
  }
  
  R2_sequence <- matrix(data = c(x, y), ncol = 2)
  
  
  return(R2_sequence)
}
