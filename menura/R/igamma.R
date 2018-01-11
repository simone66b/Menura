dinvgamma <- function(x, shape, scale = 1, log=TRUE) {
 # error checking
 if (shape <= 0 | scale <= 0) {
    stop("Shape or scale parameter negative in dinvgamma().\n")
 }
 alpha <- shape
 beta <- scale
 # done on log scale to allow for large alphas and betas
 log_density <- alpha * log(beta) - lgamma (alpha) -
       (alpha + 1) * log(x) - (beta / x)
 if (log) {
   return(log_density)
 } else {
   return(exp(log_density))
 }
}

rinvgamma <- function(n, shape, scale = 1) {
    return(1 / rgamma(n = n, shape = shape, scale = 1 / scale))
}
