update_tree <- function(lst, tr, tipdata, rt_value, N, method,
                         theta, model, mcmc_type, ...) {

  new_lst   <- lst
  n_tips <- length(tr$tip.label)
  rt_node <- n_tips + 1

  dotslist <- list(...)
  if ("pred.corr" %in% names(dotslist)) {
    pred.corr <- dotslist$pred.corr
  } else {
    pred.corr <- FALSE
  }

  if (method == "milstein") {
    pred.corr <- TRUE
  }

  if (pred.corr) {
    if (!exists("dx_diffusion", model)) {
      model$dx_diffusion <- D(model$diffusion, "x")
    }
  } else {
      model$dx_diffusion <- quote(NULL)
  }

  rt_node_dist <- ape::dist.nodes(tr)[rt_node, ]

  sde_edges <- function(tr, node, X0, t0) {

    daughters <- tr$edge[which(tr$edge[, 1] == node), 2]
    for (d_ind in 1:2) {
      edge <- which((tr$edge[, 1] == node) & (tr$edge[, 2] == daughters[d_ind])) 
      drift <- as.expression(force(eval(substitute(
        substitute(e, list(alpha = theta[edge, "alpha"], mu = theta[edge, "mu"],
                   sigma = theta[edge, "sigma"])), list(e = model$drift)))))

      diffusion <- as.expression(force(eval(substitute(
        substitute(e, list(alpha = theta[edge, "alpha"], mu = theta[edge, "mu"],
                   sigma = theta[edge, "sigma"])), list(e = model$diffusion)))))

      diffusion_x <- as.expression(force(eval(substitute(
        substitute(e, list(alpha = theta[edge, "alpha"], mu = theta[edge, "mu"],
                   sigma = theta[edge, "sigma"])),
                    list(e = model$dx_diffusion)))))

      tE <- tsp(lst[[edge]])[2]
      n_steps <- length(lst[[edge]]) - 1
      if (n_steps > 0) {
        new_lst[[edge]] <<- sde::sde.sim(X0 = X0, t0 = t0, T = tE, N = n_steps,
                                         method = method,
                                         drift = drift,
                                         sigma = diffusion,
                                         sigma.x = diffusion_x,
                                         pred.corr = pred.corr)
      } else {
        new_lst[[edge]][1] <<- X0
      }
      if (daughters[d_ind] > n_tips) {
        sde_edges(tr, daughters[d_ind], new_lst[[edge]][n_steps + 1], tE)
      }
    }
  }

  sde_edges(tr, rt_node, X0 = rt_value, t0 = 0)

  if (mcmc_type == "tanner-wong") {
    lst <- new_lst
    data_accept <- 1
  } else {

    # if (method == "milstein") {
    #   dc_fn <- sde::Elerian
    # } else if (method == "euler") {
    #   dc_fn <- sde::dcEuler
    # }

    data_accept <- 0
    logl_curr <- logl_new <- 0
    for (i_tip in 1:n_tips) {
      edgeoftip <- ape::which.edge(tr, tr$tip.label[i_tip])
      X_end <- tipdata[i_tip]
      T_end <- rt_node_dist[i_tip]
      X_start <- lst[[edgeoftip]][length(lst[[edgeoftip]])]
      T_start <- tsp(lst[[edgeoftip]])[2]
      logl_curr <- logl_curr +
              dc_fn(x = X_end, t = T_end, x0 = X_start, t0 = T_start,
                     theta = theta[edgeoftip, ],
                     model = model,
                     log = TRUE,
                     method = method)

      X_start <- new_lst[[edgeoftip]][length(new_lst[[edgeoftip]])]
      T_start <- tsp(new_lst[[edgeoftip]])[2]
      logl_new <- logl_new +
              dc_fn(x = X_end, t = T_end, x0 = X_start, t0 = T_start,
                     theta = theta[edgeoftip, ],
                     model = model,
                     log = TRUE,
                     method = method)
    }
    accept_prob <- exp(sum(logl_new - logl_curr))

    if (is.na(accept_prob))
      accept_prob <- 0

    if (runif(1) <= accept_prob) {
      lst <- new_lst
      data_accept <- 1
    }
  }
  return(list(lst = lst, data_accept = data_accept))
}
