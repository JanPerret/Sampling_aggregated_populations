
# a duplicate of the function SDraw::halton

halton <- function (n, dim = 1, start = 0, bases = NULL) {
  
  first.primes <- c(2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 
                    37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 
                    101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 
                    157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 
                    223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 
                    277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 
                    349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 
                    419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 
                    479, 487, 491, 499, 503, 509, 521, 523, 541)
  
  if (dim > length(first.primes)) {
    first.primes <- primes(dim)
  }
  if (length(start) == 1) {
    start <- rep(start, dim)
  }
  if (length(start) != dim) {
    stop("The start vector must either have length 1 or length equal to the number of dimensions")
  }
  if (is.null(bases)) {
    bases <- first.primes[1:dim]
  }
  pos <- t(sapply(start, FUN = function(x, k) {
    x:(x + k - 1)
  }, k = n))
  if ((n == 1) & all(start == 0)) {
    return(matrix(0, 1, dim))
  }
  n.sum.terms <- max(floor(log(start + n - 1)/log(bases)) + 
                       1)
  ans <- matrix(0, nrow(pos), ncol(pos))
  for (j in 0:n.sum.terms) {
    ans <- ans + ((pos%/%bases^j)%%bases)/bases^(j + 1)
  }
  if (n == 1) {
    return(matrix(ans, 1))
  }
  else {
    return(t(ans))
  }
}
