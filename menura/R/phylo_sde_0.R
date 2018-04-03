

phylo_sde_0 <- function(fossils, tr, rt_value, N, theta, model, method, ...) {

  lst <- list()
  tr <- ftr
  n_tips <- length(tr$tip.label)
  rt_node <- n_tips + 1

  #only relevant if specified
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

  sde_edges <- function(tr, node, X0, t0) {
    
    #node = rt_node
    # preceeding nodes
    daughters <- tr$edge[which(tr$edge[, 1] == node), 2]
  
    if (any(daughters %in% fossils)) {
      # do not use fossil edge (length = 0), use the sister node edge
      
      edge <- which((tr$edge[,1] == node) & !(tr$edge[, 2] %in% fossils))
      f_edge <- which((tr$edge[,1] == node) & (tr$edge[, 2] %in% fossils))
      root <- tr$edge[edge, 2]
      lst[[f_edge]] <<- 0
      
      
      drift <- as.expression(force(eval(substitute(substitute(e,
                              list(alpha = theta[edge, "alpha"],
                              mu = theta[edge, "mu"],
                              sigma = theta[edge, "sigma"])),
                              list(e = model$drift)))))
      print(drift)
      #diffusion = expression (1)
      # model$diffusion = sigma
      diffusion <- as.expression(force(eval(substitute(substitute(e,
                                  list(alpha = theta[edge, "alpha"],
                                  mu = theta[edge, "mu"],
                                  sigma = theta[edge, "sigma"])),
                                  list(e = model$diffusion)))))
      
      #diffusion_x = 0
      #model$dx_diffusion  
      diffusion_x <- as.expression(force(eval(substitute(substitute(e,
                                    list(alpha = theta[edge, "alpha"],
                                    mu = theta[edge, "mu"],
                                    sigma = theta[edge, "sigma"])),
                                    list(e = model$dx_diffusion)))))
      
      #number of steps is the length of the edge times the given frequency N (100)
      #####May cause an issue as the edge length of a fossil is 0, but maybe not considering you can have polytomys
      n_steps <- tr$edge.length[edge] * N 
      #time end is time start plus the length of the given edge
      ####This will also be 0 for fossils
      tE <- t0 + tr$edge.length[edge]
      
      #runs sde.sim
      lst[[edge]] <<- sde::sde.sim(X0 = X0, t0 = t0, T = tE, N = n_steps,
                                   method = method,
                                   drift = drift,
                                   sigma = diffusion,
                                   sigma.x = diffusion_x,
                                   pred.corr = pred.corr)
      tE <- tsp(lst[[edge]])[2]
      
        if (root > n_tips) {
          sde_edges(tr, root, lst[[edge]][n_steps + 1], tE)
       }
        
     }else{     
        
        for (d_ind in 1:2) {
        
        #what are node and d_ind
        #I think this is how you traverse the tree
        # where does (col 1 = rt_node && col 2 = daughters[d_ind])
        edge <- which((tr$edge[, 1] == node) & (tr$edge[, 2] == daughters[d_ind]))
        
        
        
        #mu -> long term mean level
        #alpha -> 
        #drift = expression (0.1 *(0-x))
        # model$drift = alpha * (mu - x)
        #sets the parameters as 2 vectors/lists
        drift <- as.expression(force(eval(substitute(substitute(e,
                                                                list(alpha = theta[edge, "alpha"],
                                                                     mu = theta[edge, "mu"],
                                                                     sigma = theta[edge, "sigma"])),
                                                     list(e = model$drift)))))
        
        #diffusion = expression (1)
        # model$diffusion = sigma
        diffusion <- as.expression(force(eval(substitute(substitute(e,
                                                                    list(alpha = theta[edge, "alpha"],
                                                                         mu = theta[edge, "mu"],
                                                                         sigma = theta[edge, "sigma"])),
                                                         list(e = model$diffusion)))))
        
        #diffusion_x = 0
        #model$dx_diffusion  
        diffusion_x <- as.expression(force(eval(substitute(substitute(e,
                                                                      list(alpha = theta[edge, "alpha"],
                                                                           mu = theta[edge, "mu"],
                                                                           sigma = theta[edge, "sigma"])),
                                                           list(e = model$dx_diffusion)))))
        
        #number of steps is the length of the edge times the given frequency N (100)
        #####May cause an issue as the edge length of a fossil is 0, but maybe not considering you can have polytomys
        n_steps <- tr$edge.length[edge] * N 
        #time end is time start plus the length of the given edge
        ####This will also be 0 for fossils
        tE <- t0 + tr$edge.length[edge]
        
        #runs sde.sim
        lst[[edge]] <<- sde::sde.sim(X0 = X0, t0 = t0, T = tE, N = n_steps,
                                     method = method,
                                     drift = drift,
                                     sigma = diffusion,
                                     sigma.x = diffusion_x,
                                     pred.corr = pred.corr)
        tE <- tsp(lst[[edge]])[2]
        
        if (daughters[d_ind] > n_tips) {
          sde_edges(tr, daughters[d_ind], lst[[edge]][n_steps + 1], tE)
        }
       
      }
        
     } 
    
  }
  sde_edges(tr, rt_node, X0 = rt_value, t0 = 0)

 
 #print(lst)  
 
   
  # Remove tip values (we have observed tip values)
  node_len <- ape::node.depth.edgelength(tr)
  print(node_len)
  for (nthtip in 1:n_tips) {
  
  
      nEdge <- which(tr$edge[, 2] == nthtip)
      ntsp <- tsp(lst[[nEdge]])
      print(ntsp)
      
      # If the end time for the node edges is less than T - 1/N,
      # the last sample is removed from the simulated data,
      # if there is more than one sample.
    if (nthtip %in% fossils){
     ntsp[2] <- 0
     lst[[nEdge]] <- 0
      }
      if (ntsp[2] > (node_len[nthtip] - 1 / N)) {
        if (length(lst[[nEdge]]) > 1) {
          if (length(lst[[nEdge]]) == 2) {
            lst[[nEdge]] <- ts(lst[[nEdge]][1], start = ntsp[1], end = ntsp[1],
                               frequency = N)
          } else {
            lst[[nEdge]] <- window(lst[[nEdge]], start = ntsp[1],
                                   end = ntsp[2] - 1 / N)
          }
        }
      }
  }
  # fossils in list are given value of 0, NULL causes problems if the fossil is right justified
  return(lst)
}

