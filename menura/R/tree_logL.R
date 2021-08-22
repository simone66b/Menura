#'tree_logL
#'
#'Calculate the loglikelihood of a Cox-Ingersoll-Ross Diffusion Process in the Tree of Life
#'
#'  This function caluclates the log likelihood.
#'
#'@param tr An object of class \code{phylo} from the ape package.  In this version the CIR
#'process parameters alpha, mu, and sigma for each branch are within vectors in the same order
#'as the edge (branch) labelling.
#'@param tipdata A numeric vector containing tip values in the same order as the tip labels in \code{tr$tip.label}.
#'@param lst A list of time series projects for each simulated path equal to the length of
#'the number of branches in the \code{tr} object.
#'@param alpha Set to NULL if alpha is to be estimated,
#'otherwise set to a numeric value or a numeric vector specifying
#'the value of the parameter for all the branches/edges.
#'In the latter case, the values must be specifying in the same order
#'as the edges in the \code{tr} object.
#'@param mu As \code{alpha}.
#'@param sigma As \code{alpha}.
#'@param model A list containing drift, diffusion, and the partial differentiation of diffusion as quoted
#'expressions using method quote. For the Euler scheme
#'the drift coefficient as \code{drift}, the diffusion coefficient as
#'\code{diffusion}, and the partial differentiation of
#'\code{diffusion} by \code{x} as \code{dx_diffusion} is required.
#'@param method Specified as either "euler" or "milstein."
#'@param fossils A numeric vector containing the tip values for every fossil added to the tree.
#'@seealso \code{\link{fossil_id}}
#'@param ... Not used.
#'
#'@return logL, an integer.
#'
#'@export

tree_logL <- function(tr, tipdata, lst, alpha, mu, sigma, model,
              method, fossils=NULL, ...) {
  

tipdata <- as.numeric(tipdata)

n_tips  <- length(tr$tip.label)
rt_node <- n_tips + 1
logL <- NULL


#list of the distance of each node from the root node
rt_node_dist <- ape::dist.nodes(tr)[rt_node, ]



##Calculates the log likelihood at each edge given the specific parameters of the edge
##returns value logL which is the sum of the logLikelihood of all edges
logL_edges <- function (fossils, node, tr, tipdata, lst, alpha, mu, sigma, model) {
  #The daughter nodes are those that at first branch directily from rt_node, node is updated for subsequent splits
    daughters <- tr$edge[which(tr$edge[, 1] == node), 2]
    
  ###if a daughter node is in fossils, apply the tipdata of the fossil to the current node
  ###continue to estimate along the other side of the tree
  
    if (any(daughters %in% fossils)){
    #edge that does not belong to the fossil
       edge <- which((tr$edge[,1] == node) & !(tr$edge[, 2] %in% fossils))
       theta <- c(alpha[edge], mu[edge], sigma[edge])
       
    #fossil edge that will take 0   
       f_edge <- which((tr$edge[,1] == node) & (tr$edge[, 2] %in% fossils))
       
       logL[f_edge] <<- 0
    #tip associated with fossil
       f <- daughters[which(daughters %in% fossils)]
    #node associated with edge   
       n <- tr$edge[edge, 2]
       reroot <- n
    
    #internal branch with fossil treated as edge with tipdata
       logL[edge] <<- logl_fn(X = lst[[edge]], theta = theta,
                           model = model, log = TRUE, method = method) +
                      dc_fn(x = tipdata[f],
                           t = rt_node_dist[n],
                           #t = 1,
                           x0 = lst[[edge]][length(lst[[edge]])],
                           #t0 = tsp(lst[[edge]])[2],
                           t0 = tsp(lst[[edge]])[1],
                           theta = theta,
                           model = model,
                           log = TRUE,
                           method = method)
      
       #recursive call, if we are not at a tip, keep going to the end of the tree
       #uses the non-fossil node as the new root
       if (n > n_tips){
          logL_edges(fossils, reroot, tr, tipdata, lst, alpha, mu, sigma, model)
       }
       
    }else{
     
      #if there is no fossil, run both daughter node edges
      #assumes bifurcating tree
      for (ind_d in 1:2) {
         edge <- which((tr$edge[, 1] == node) & (tr$edge[, 2] == daughters[ind_d]))
         theta <- c(alpha[edge], mu[edge], sigma[edge])
         
         #if the edge is an internal branch
         if (daughters[ind_d] > n_tips) {
            logL[edge] <<- logl_fn(X = lst[[edge]], theta = theta,
                               model = model, log = TRUE,
                               method = method)
     
        #must be a tip
        #logL of a tip and its edge is the sum of the conditional density of the diffusion process and the logl 
         } else {

             if (is.null(lst[[edge]])) logL
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
      #reset the node to the "new root" and rerun logL_edges until you get to the tips
         if (daughters[ind_d] > n_tips) {
            logL_edges(fossils, daughters[ind_d], tr, tipdata, lst, alpha, mu, sigma, model)
         }
      }
    }
 }  
  


logL_edges(fossils, rt_node, tr, tipdata, lst, alpha, mu, sigma, model)
#logL allows us to sum the likelihoods rather than take the product

return(sum(logL))
}
