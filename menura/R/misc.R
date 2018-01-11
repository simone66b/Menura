print.menura <- function(x, ...) {
  cat("Call:\n")
  print(x$call)

  cat("\nSummary: \n")
  print(summary(x$mcmctrace))

  cnam <- colnames(x$mcmctrace)
  for (ind in 1:length(cnam)) {
    plot(x$mcmctrace[, cnam[ind]], type = "l", ylab = cnam[ind],
         xlab = "MCMC trace")
    if (ind != length(cnam))
    readline(prompt = "Press [enter] to continue")
  }
}

summary.menura <- function(object, ...) {
  cat("Call:\n")
  print(object$call)

  cat("\nSummary: \n")
  print(summary(object$mcmctrace))
}

plot.menura <- function(x, ...) {
  cnam <- colnames(x$mcmctrace)
  for (ind in 1:length(cnam)) {
    plot(x$mcmctrace[, cnam[ind]], type = "l", ylab = cnam[ind],
         xlab = "MCMC trace")
    if (ind != length(cnam))
    readline(prompt = "Press [enter] to continue")
  }
}
