### Calculates the Log Likelihood of the Diffusion Process in a Tree of Life
## Method can only be 'euler' or 'milstein'
#see help file for argument descriptors
tree_logL <- function(tr, tipdata, lst, alpha, mu, sigma, model,
              method, ...) {
tipdata <- as.numeric(tipdata)
n_tips  <- length(tr$tip.label)
rt_node <- n_tips + 1
logL <- NULL

#list of the distance of each node from the root node
rt_node_dist <- ape::dist.nodes(tr)[rt_node, ]
#rt_node_dist


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

##Calculates the log likelihood at each edge given the specific parameters of the edge
##returns value logL which is the sum of the logLikelihood of all edges
logL_edges <- function (node, tr, tipdata, lst, alpha, mu, sigma, model) {
  #The daughter nodes are those that at first branch directily from rt_node, node is updated for subsequent splits
  daughters <- tr$edge[which(tr$edge[, 1] == node), 2]
  
  for (ind_d in 1:2) {
  
    edge <- which((tr$edge[, 1] == node) & (tr$edge[, 2] == daughters[ind_d]))
    theta <- c(alpha[edge], mu[edge], sigma[edge])

    #if the value of the node is greater than the total number of tips, the node is not a tip
    #this is logical due to the way trees are labelled (tips first)
    #if not a tip, calculate LogL and jump to recursive call
    if (daughters[ind_d] > n_tips) {
      logL[edge] <<- logl_fn(X = lst[[edge]], theta = theta,
                            model = model, log = TRUE,
                            method = method)
      #must be a tip
      #calls 2 functions within dc_fn
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
    #recursive call
    #reset the node to the "new root" 
    if (daughters[ind_d] > n_tips) {
      logL_edges(daughters[ind_d], tr, tipdata, lst, alpha, mu, sigma, model)
    }
  }
}

logL_edges(rt_node, tr, tipdata, lst, alpha, mu, sigma, model)
#logL allows us to sum the likelihoods rather than take the product
return(sum(logL))
}
logL
