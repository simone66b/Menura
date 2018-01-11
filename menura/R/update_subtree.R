update_subtree <- function(lst, tr, tipdata, rt_value, N, method = "euler",
                            theta, model, mcmc_type = "tanner-wong", ...) {
  # Changes
  # log L's  are calculated for subtr_tips instead of a clade
  # update from any node regardless edge length
  # ** if the selected subtree is has small edges, then no
  #    data update is done.

  # get the method from ...'s
  new_lst   <- lst
  subtr_tips <- NULL
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
      edge <- which((tr$edge[, 1] == node) & (tr$edge[ ,2] == daughters[d_ind]))
      drift <- as.expression(force(eval(substitute(
                  substitute(e, list(alpha = theta[edge, "alpha"],
                    mu = theta[edge, "mu"],
                    sigma = theta[edge, "sigma"])),
                  list(e = model$drift)))))

      diffusion <- as.expression(force(eval(substitute(
                    substitute(e, list(alpha = theta[edge, "alpha"],
                      mu = theta[edge, "mu"],
                      sigma = theta[edge, "sigma"])),
                    list(e = model$diffusion)))))

      diffusion_x <- as.expression(force(eval(substitute(
                      substitute(e, list(alpha = theta[edge, "alpha"],
                        mu = theta[edge, "mu"],
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
      } else {
        subtr_tips <<- c(subtr_tips, daughters[d_ind])
      }
    }
  }

  rnode <- sample((n_tips + 1):max(tr$edge), 1)

  # If seelcted subtree doesn't contain large enough edges
  sub_trs <- ape::subtrees(tr)
  if (max( sub_trs[[rnode - n_tips]]$edge.length ) * N < 2) {
    return(list(lst = lst, data_accept = NA))
  }

  # If the selected node is the root, we would update whole tree.
  if (rnode == rt_node) {
    new_lst <- update_tree(lst = lst, tr = tr, tipdata = tipdata,
                           rt_value = rt_value, N = N, method = method,
                           theta = theta, model = model, mcmc_type = mcmc_type,
                           ...)
    return(list(lst = new_lst$lst, data_accept = new_lst$data_accept))
  } else {

    # Update the sub-tree
    redge <- which(tr$edge[, 1] == rnode)[1]
    sde_edges(tr, rnode, X0 = lst[[redge]][1], t0 = tsp(lst[[redge]])[1])

    if (mcmc_type == "tanner-wong") {
      lst <- new_lst
      data_accept <- 1
    } else {

      # if (method == "milstein") {
      #   dc_fn <- sde::Elerian
      # } else if (method == "euler") {
      #   dc_fn <- sde::dcEuler
      # }

      # Calculate the acceptance probability
      data_accept <- 0
      curr_path_logLs <- new_path_logLs <- 0
      for (i_tip in subtr_tips) {
        i_tip_edge <- ape::which.edge(tr, tr$tip.label[i_tip])
        X_end <- tipdata[i_tip]
        T_end <- rt_node_dist[i_tip]
        X_start <- lst[[i_tip_edge]][length(lst[[i_tip_edge]])]
        T_start <- tsp(lst[[i_tip_edge]])[2]
        curr_path_logLs <- curr_path_logLs +
                  dc_fn(x = X_end, t = T_end, x0 = X_start, t0 = T_start,
                  theta = theta[i_tip_edge, ],
                  model = model,
                  log = TRUE,
                  method = method)

        X_start <- new_lst[[i_tip_edge]][length(new_lst[[i_tip_edge]])]
        T_start <- tsp(new_lst[[i_tip_edge]])[2]
        new_path_logLs <- new_path_logLs +
                  dc_fn(x = X_end, t = T_end, x0 = X_start, t0 = T_start,
                  theta = theta[i_tip_edge, ],
                  model = model,
                  log = TRUE,
                  method = method)
      }
      accept_prob <- exp(sum(new_path_logLs - curr_path_logLs))

      if (is.na(accept_prob))
        accept_prob <- 0

      if (runif(1) <= accept_prob) {
        lst <- new_lst
        data_accept <- 1
      }
    }
    return(list(lst = lst, data_accept = data_accept))
  }
}
