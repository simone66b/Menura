sde_model <- function(model, rt_value, tipdata, ...) {

  if (class(model) == "character") {
    if (model == "OU") {
      M <- list()
      M$d <- function(t, x, theta) {
                theta[1] * (theta[2] - x)
      }
      M$s <- function(t, x, theta) {
        theta[3]
      }
      M$s_x <- function(t, x, theta) {
        0
      }
      M$drift <- quote(alpha * (mu - x))
      M$diffusion <- quote(sigma)
      M$dx_diffusion <- quote(0)
      attr(M, "sde_class") <- "OU"

    } else if (model == "CIR") {

      M <- list()
      M$d <- function(t, x, theta) {
        ((theta[1] * theta[2] - 0.25 * theta[3]^2) / (2 * x)) - theta[1] * x / 2
      }
      M$s <- function(t, x, theta) {
        0.5 * theta[3]
      }
      M$s_x <- function(t, x, theta) {
        0
      }
      M$drift <- quote(((alpha * mu - 0.25 * sigma^2) / (2 * x)) - alpha * x / 2)
      M$diffusion <- quote(sigma / 2)
      M$dx_diffusion <- quote(0)
      tipdata  <- sqrt(tipdata)
      rt_value <- sqrt(rt_value)
      attr(M, "sde_class") <- "CIR"

    } else if (model == "Beta") {

      if (max(abs(tipdata)) > 1)
        stop("For the Beta model, tip values can take values in ",
              " the -1 to 1 range only.")

      M <- list()
      M$d <- function(t, x, theta) {
        ((2 * theta[1] * theta[2] - theta[2]) / sin(x)) -
          (theta[2] - theta[3] / 2) * (1 / tan(x))
      }
      M$s <- function(t, x, theta) {
        theta[3]
      }
      M$s_x <- function(t, x, theta) {
        0
      }
      M$drift <- quote((2 * alpha * mu - alpha) / sin(x) +
                          (alpha - sigma / 2) * (1 / tan(x)))
      M$diffusion <- quote(sigma)
      M$dx_diffusion <- quote(0)
      tipdata  <- 2 * asin(tipdata)
      rt_value <- 2 * asin(rt_value)
      attr(M, "sde_class") <- "Beta"
    } else {
      stop("If the model ", model, " is not any of OU, CIR, Beta,",
            "then the model parameters must be specified in a list.",
           " See documentation.", call. = FALSE)
    }
  } else {
    if (class(model) == "list") {
      M <- model
      Mparam <-  c("d", "s", "drift", "diffusion", "dx_diffusion")
      Mtest <- Mparam %in% names(M)
      if (!all(Mtest)) {
         stop(paste("The parameters ", Mparam[which(!Mtest)],
                    " are missing in the model"), call. = FALSE)
      }
      if( is.null(attr(M, "sde_class")) ) attr(M, "sde_class") <- "general"
    } else {
      stop("Incorrect specification of model parameter. See documentation.",
      call. = FALSE)
    }
  }

  return(list(M=M, tipdata=tipdata, rt_value = rt_value))
}
