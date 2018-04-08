# dc_fn Sets the method up
# milstein & euler only, stop message if not 

##Approximates the conditional density of a diffusion process (stochastic process)
# x is a vector of quantiles
# t is lag or time
#x0 the value of the process at t0
#t0 initial time
#theta parameter of the process
#log, if TRUE, probabilities are given as log(p)
#d is the drift coefficient as a fxn
#s is the diffusion coefficient as a fxn
dc_fn <- function(x, t, x0, t0, theta, model, log=TRUE, method) {

  #s, d, sx are all functions of t, x, and theta
  if (method == "milstein") {
    ret <- sde::dcElerian(x, t, x0, t0, theta, d = model$d,
                      s = model$s, sx = model$s_x, log=log)
    return(ret)
  } else if (method == "euler") {
    ret <- sde::dcEuler(x, t, x0, t0, theta, d = model$d,
                      s = model$s, log=log)
    # if ((is.na(ret))){
    #  print("dcfn")
    #   print(ret)
    #   print("x")
    #   print(x)
    #   print("t")
    #   print(t)
    #   print("x0")
    #   print(x0)
    #   print("t0")
    #   print(t0)
    # }
   
    return(ret)
  }

  stop("Presently, method has to be either \"euler\" or \"milstein\"")
}

# Chooses the appropriate log likelihood method 
logl_fn <- function (X, theta, model, log = TRUE, method) {

  if (length(X) == 1) {
    return (0)
  }
  if (method == "milstein") {
    l0 <- 0
    tvec <- time(X)
    for (j in 2:length(X))
      l0 <- l0 + dc_fn(x=X[j], t=tvec[j], x0=X[j-1], t0=tvec[j-1],
                  theta, model, log=log, method = method)

    return(l0)

  } else if (method == "euler") {
    ret <- sde::EULERloglik(X, theta, model$d, model$s, log=log)
    
    # if (is.na(ret)){
    #   print("logfn")
    #   print(ret)
    #   print("X")
    #   print(X)
    #   print("theta")
    #   print(theta)
    # }
    
    return(ret)
  }

  stop("Presently, method has to be either \"euler\" or \"milstein\"")

}

# if (method == "milstein") {
#
#   dc_fn <- function(x, t, x0, t0, theta, model, log=TRUE) {
#     sde::dcElerian(x, t, x0, t0, theta, d = model$d,
#                       s = model$s, sx = model$s_x, log=log)
#   }
#
#   logl_fn <- function (X, theta, model, log = TRUE) {
#     l0 <- 0
#     tvec <- time(X)
#     for (j in 2:length(X))
#       l0 <- l0 + dc_fn(x=X[j], t=tvec[j], x0=X[j-1], t0=tvec[j-1],
#                   theta, model, log=TRUE)
#     l0
#   }
#
# } else if (method == "euler") {
#
#   dc_fn <-  <- function(x, t, x0, t0, theta, model, log=TRUE) {
#     sde::dcEuler(x, t, x0, t0, theta, d = model$d, s = model$s, log=log)
#   }
#
#   logl_fn <- function (X, theta, model, log = TRUE) {
#     sde::EULERloglik(X, theta, model$d, model$s, log=log)
#   }
# }

