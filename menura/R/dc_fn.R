dc_fn <- function(x, t, x0, t0, theta, model, log=TRUE, method) {

  if (method == "milstein") {
    ret <- sde::dcElerian(x, t, x0, t0, theta, d = model$d,
                      s = model$s, sx = model$s_x, log=log)
    return(ret)
  } else if (method == "euler") {
    ret <- sde::dcEuler(x, t, x0, t0, theta, d = model$d,
                      s = model$s, log=log)
    return(ret)
  }

  stop("Presently, method has to be either \"euler\" or \"milstein\"")
}

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
