order_tree <- function(tr, alpha, mu, sigma, priors) {

  n_edges <- length(tr$edge.length)
  n_tips  <- length(tr$tip.label)
  old_tr <- NULL

  if (!any(attributes(tr) %in% "cladewise")) {
    old_tr <- tr
    tr <- reorder(old_tr, "cladewise")

    reorder_edge <- unlist(
      lapply(
        1:nrow(old_tr$edge),
        function(i)
          which(apply(old_tr$edge, 1, identical, tr$edge[i,])))
    )
  }

  para2est <- NULL
  if (!is.null(alpha)) {
    #alpha <- Alpha
    if (length(alpha) == 1)
      alpha <- rep(alpha, n_edges)
  } else {
    alpha <- rep(priors[["alpha"]]$rf(1), n_edges)
    para2est <- c(para2est, alpha = "alpha")
  }

  if (!is.null(mu)) {
    #mu <- Mu
    if (length(mu) == 1)
      mu <- rep(mu, n_edges)
  } else {
    mu <- rep(priors[["mu"]]$rf(1), n_edges)
    para2est <- c(para2est, mu = "mu")
  }

  if (!is.null(sigma)) {
    #sigma <- Sigma
    if (length(sigma) == 1)
      sigma <- rep(sigma, n_edges)
  } else {
    sigma <- rep(sqrt(priors[["sigma"]]$rf(1)), n_edges)
    para2est <- c(para2est, sigma = "sigma")
  }

  if (length(alpha) != n_edges)
    stop("If Alpha given, it must either be a numeric value or",
         " a vector of numeric of length equivalnet to number of",
         " edges in the tree.", call. = FALSE)

  if (length(mu) != n_edges)
    stop("If Mu given, it must either be a numeric value or",
         " a vector of numeric of length equivalnet to number of",
         " edges in the tree.", call. = FALSE)

  if (length(sigma) != n_edges)
    stop("If Sigma given, it must either be a numeric value or",
         " a vector of numeric of length equivalnet to number of",
         " edges in the tree.", call. = FALSE)

  theta <- cbind(alpha = alpha, mu = mu, sigma = sigma)

  return(list(tr = tr, theta = theta, para2est = para2est, old_tr = old_tr))
}
