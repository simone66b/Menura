#'phylo_sde
#'
#'Simulate a Cox-Ingersoll-Ross Diffusion Process in the Tree of Life
#'
#'@param tr An object of class \code{phylo} from the ape package.  In this version the CIR
#'process parameters alpha, mu, and sigma for each branch are within vectors in the same order
#'as the edge (branch) labelling.
#'@param rt_value Value at the root of \code{tr}.
#'@param N Data imputation frequency.
#'@param theta Matrix of parameter values for each edge of the tree.
#'@param model A list containing drift, diffusion, and the partial differentiation of diffusion as quoted
#'expressions using method quote. For the Euler scheme
#'the drift coefficient as \code{drift}, the diffusion coefficient as
#'\code{diffusion}, and the partial differentiation of
#'\code{diffusion} by \code{x} as \code{dx_diffusion} is required.
#'See the Examples.
#'@param method Specified as either "euler" or "milstein."
#'@param fossils A numeric vector containing the tip values for every fossil added to the tree.
#'@param ... Not used.
#'
#'@return A list of time series projects for each simulated path equal to the length of
#'the number of branches in the \code{tr} object.
#'
#'@export


phylo_sde <- function(tr, rt_value, N, theta, model, method,  traits, fossils=NULL, ...) {

  #stores list of points for each simulated edge
  lst <- list()
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

    sde_edges <- function(fossils, tr, node, X0, t0, traits) {
    
    #node is the root node
    # preceeding nodes
    daughters <- tr$edge[which(tr$edge[, 1] == node), 2]
  ## if (length(daughters) == 0) return()
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
      n_steps <- tr$edge.length[edge] * N 
      #time end is time start plus the length of the given edge
      tE <- t0 + tr$edge.length[edge]
      X0 <- traits[fossils[which(fossils %in% daughters)]]
                                        #runs sde.sim
       ##  print(t0)
      lst[[edge]] <<- sde::sde.sim(X0 = X0, t0 = t0, T = tE, N = n_steps,
                                  method = method,
                                   drift = drift,
                                   sigma = diffusion,
                                   sigma.x = diffusion_x,
                                   pred.corr = pred.corr
                                   )
      tE <- tsp(lst[[edge]])[2]
      
        if (root > n_tips) {
          sde_edges(fossils=fossils, tr=tr, node=root, X0=lst[[edge]][n_steps + 1], t0=tE, traits=traits)
       }
        
     }else{     
        
        for (d_ind in 1:2) {
       #identify the branch for simulation
        edge <- which((tr$edge[, 1] == node) & (tr$edge[, 2] == daughters[d_ind]))
        
        
        
        #drift = expression (0.1 *(0-x))
        # model$drift = alpha * (mu - x)
        #pulls specific edge parameters from vector theta
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
       
        n_steps <- tr$edge.length[edge] * N 
        #time end is time start plus the length of the given edge
      
        tE <- t0 + tr$edge.length[edge]
        
                                        #runs sde.sim
           if (t0 == tE) {
             lst[[edge]] <- 0
            } else { 
                
        lst[[edge]] <<- sde::sde.sim(X0 = X0, t0 = t0, T = tE, N = n_steps,
                                   method = method,
                                     drift = drift,
                                     sigma = diffusion,
                                     sigma.x = diffusion_x,
                                     pred.corr = pred.corr
                                     )
        }
            tE <- tsp(lst[[edge]])[2]
            
        #recursive call until a tip is met
        if (daughters[d_ind] > n_tips && lengths(lst[[edge]]) > 0) {
          sde_edges(fossils=fossils, tr=tr, node=daughters[d_ind], X0=lst[[edge]][n_steps + 1],t0= tE, traits=traits)
        } else {
            sde_edges(fossils=fossils, tr=tr, node=daughters[d_ind], X0=X0, t0=tE, traits=traits)
       
        }
        
        } 
    
     }
    }

  sde_edges(fossils=fossils, tr=tr, node=rt_node, X0 = rt_value, t0 = 0, traits=traits)

   
  # Remove tip values (we have observed tip values)
  # fossils in list are given value of 0, NULL causes problems if the fossil is right justified
  return(lst)
}





############################################################################
 sde_edges <- function(fossils, tr, node, X0, t0, traits) {
    
    #node is the root node
    # preceeding nodes
    daughters <- tr$edge[which(tr$edge[, 1] == node), 2]
  ## if (length(daughters) == 0) return()
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
      n_steps <- tr$edge.length[edge] * N 
      #time end is time start plus the length of the given edge
      tE <- t0 + tr$edge.length[edge]
      X0 <- traits[fossils[which(fossils %in% daughters)]]
                                        #runs sde.sim
      lst[[edge]] <<- sde::sde.sim(X0 = X0, t0 = t0, T = tE, N = n_steps,
                                   drift = drift,
                                   sigma = diffusion,
                                   sigma.x = diffusion_x, method="euler", pred.corr=F
                                   )
      tE <- tsp(lst[[edge]])[2]
      
        if (root > n_tips) {
          sde_edges(fossils=fossils, tr=tr, node=root, X0=lst[[edge]][n_steps + 1], t0=tE, traits=traits)
       }
        
     }else{     
        
        for (d_ind in 1:2) {
       #identify the branch for simulation
        edge <- which((tr$edge[, 1] == node) & (tr$edge[, 2] == daughters[d_ind]))
        
        
        
        #drift = expression (0.1 *(0-x))
        # model$drift = alpha * (mu - x)
        #pulls specific edge parameters from vector theta
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
       
        n_steps <- tr$edge.length[edge] * N 
        #time end is time start plus the length of the given edge
      
        tE <- t0 + tr$edge.length[edge]
        
                                        #runs sde.sim
           if (t0 == tE) {
             lst[[edge]] <<- 0
            } else { 
                
        lst[[edge]] <<- sde::sde.sim(X0 = X0, t0 = t0, T = tE, N = n_steps,
                                     drift = drift, method="euler",
                                     sigma = diffusion,
                                     sigma.x = diffusion_x, pred.corr=F
                                     )
                }
            tE <- tsp(lst[[edge]])[2]
            
        #recursive call until a tip is met
        if (daughters[d_ind] > n_tips) {
          sde_edges(fossils=fossils, tr=tr, node=daughters[d_ind], X0=lst[[edge]][n_steps + 1],t0= tE, traits=traits)
        } else {
            sde_edges(fossils=fossils, tr=tr, node=daughters[d_ind], X0=X0, t0=tE, traits=traits)
       
        }
        
        } 
    
     }
}
sde_edges(fossils=fossils, tr=tr, node=6, X0 = rt_value, t0 = 0, traits=traits)
