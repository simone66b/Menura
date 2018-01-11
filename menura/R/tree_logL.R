tree_logL <- function(tr, tipdata, lst, alpha, mu, sigma, model,
              method, ...) {

tipdata <- as.numeric(tipdata)
n_tips  <- length(tr$tip.label)
rt_node <- n_tips + 1
logL <- NULL
rt_node_dist <- ape::dist.nodes(tr)[rt_node, ]

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

logL_edges <- function (node, tr, tipdata, lst, alpha, mu, sigma, model) {
  daughters <- tr$edge[which(tr$edge[, 1] == node), 2]
  for (ind_d in 1:2) {
    edge <- which((tr$edge[, 1] == node) & (tr$edge[, 2] == daughters[ind_d]))
    theta <- c(alpha[edge], mu[edge], sigma[edge])

    if (daughters[ind_d] > n_tips) {
      logL[edge] <<- logl_fn(X = lst[[edge]], theta = theta,
                            model = model, log = TRUE,
                            method = method)
    } else {
      logL[edge] <<- logl_fn(X = lst[[edge]], theta = theta,
                          model = model, log = TRUE, method = method) +
                     dc_fn(x = tipdata[daughters[ind_d]],
                          t = rt_node_dist[daughters[ind_d]],
                          x0 = lst[[edge]][length(lst[[edge]])],
                          t0 = tsp(lst[[edge]])[2],
                          theta = theta,
                          model = model,
                          log = TRUE,
                          method = method)
    }
    if (daughters[ind_d] > n_tips) {
      logL_edges(daughters[ind_d], tr, tipdata, lst, alpha, mu, sigma, model)
    }
  }
}

logL_edges(rt_node, tr, tipdata, lst, alpha, mu, sigma, model)
return(sum(logL))
}
