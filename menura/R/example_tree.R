#examples{
  #set seed for random num generator to ensure simulation is reproducible
  set.seed(1)
  # install these packages
  rpkgs <- c("sde", "ape", "msm")
  lapply(rpkgs, require, character.only = TRUE)

  # Random tree with n tips
  tr <-  compute.brlen(rtree(n=4))
  
  tr$edge.length
  tr$edge
  # plot the tree
  plot(tr)
  edgelabels()
  nodelabels()
  add.scale.bar()
  
  
  ##working on manually adding fossils to tree as new branches with length 0
  tip1 <- list(edge = matrix(c(2,1),1,2), 
              tip.label = "fossil",
              edge.length = 0.0,
              Nnode = 1)
  class(tip1)<- "phylo"
  f1t <- bind.tree(tr,tip1, where = 7, position = 0.1)
  plot(f1t)
  nodelabels()
  tip2 <- list(edge = matrix(c(2,1),1,2), 
              tip.label = "fossil",
              edge.length = 0.0,
              Nnode = 1)
  class(tip2)<- "phylo"
  ftree <- bind.tree(f1t,tip2, where = 7, position = 0.1)
  plot(ftree)
  nodelabels()
  
  ftree$edge.length

  #a way to add multiple fossils given node #, tip.label, edge_length
  #PROBLEM, when entering a node at a position on the branch, nodes get reset each time
  bind.tip <- function(tree, tip.label, edge.length = NULL, where = NULL, position = NULL){
    if(is.null(where)) where <- length(tree$tip.label)+1
    tip <- list(edge = matrix(c(2,1),1,2), 
                 tip.label = "fossil",
                 edge.length = 0.0,
                 Nnode = 1)
    class(tip)<- "phylo"
    obj <- bind.tree(tree, tip, where = where)
    return(obj)
  }
  
  ###attempt at rerooting from fossils
  
  
  
  
  #step through the tree looking for fossils
  #at each node check both sides edge length, if edge length == 0 and a tip edge, 
  #add tip to vector "fossils"
  br_class <- branchClasses(ftree)
  fossils <- vector(mode = "numeric")
  
  #indicates which edges have fossils attached
  br_zero <- which(ftree$edge.length[] == 0.0)
  #indicates which tips are fossils
  f_label <- which(ftree$tip.label[] == "fossil")
 
  n_edges <- Nedge(ftree)
  

  find_fossils <- function(tr, br_class, e) 
    if ((tr$edge.length[e] == 0.0) & (br_class[e] == dead)) 
       return e
  
  
  e = 1
  n_edges = Nedge(tr)
  for(e in n_edges){
    fossil <- find_fossils(ftree, br_class, e)
  }
  
  
  br_class$brlen.dead
  plot(ftree)
  ftree$edge.length
  edgelabels()
  nodelabels()                                                    
  ftree$edge.length
  
  
  
  
  # SDE parameters
  # set to size of the length of the edge.length vector
  Nedges <- length(tr$edge.length)  
  
  dclade <- max(which(tr$edge[,1] == tr$edge[1,1])) - 1  
  
  #these parameters are vectors with size = number of edges = number of simulations
  alpha <- mu <- sigma <- rep(0, Nedges)
  
  #numeric vector of selective constraint strength for ec. branch
  alpha[1:Nedges]  <- 0.1
  
  #mu == theta as per method description
  #represents numeric vector giving ec. branch optimum
  mu[1:Nedges] <- 0
  
  #numeric vector Std-dev of random component for ec. branch
  sigma[1:Nedges] <- 1
  
  rt_value <- 0
  
  # rTraitCont simulates the evolution of a continuous character along a phylogeny
  # OU specifies model type (Brownian, Orn-Uhl, fxn)
  # model type OU is sensitive to sigma, alpha, theta 
  tipdata <- rTraitCont(tr, "OU", sigma=sigma, alpha=alpha, theta=mu,
                        root.value=rt_value)
  model <- list()
  model$d <- function (t, x, theta) {
    theta[1] * (theta[2] - x)
  }
  
  model$s <- function(t, x, theta) {
    theta[3]
  }
  
  # quote() assigns an expression to a variable instead of the solution
  model$drift <- quote(alpha * (mu - x))
  model$diffusion <- quote(sigma)
  
  # the derivative of diffusion = sigma = 0
  model$dx_diffusion <- quote(0)
  
  #theta is a vector of the parameters
  theta <- cbind(alpha=alpha, mu=mu, sigma=sigma)
  N <- 100
  
  #calls the sde function
  lst <- phylo_sde_0 (tr=tr, rt_value=rt_value, theta=theta, model=model,
                    N=N, method="euler")
  
  #calls log likelihood using Euler, approximates for diffusion process in tree
  loglike <-  tree_logL (tr=tr, tipdata=tipdata, lst=lst,
                         alpha=theta[, "alpha"],
                         mu=theta[, "mu"],
                         sigma=theta[, "sigma"], model,
                         method = "euler")
  
  
}
  
  #lst stores the point path to each node
  lst
  
  #print lst to a file
sink("phylo_sde_test_output")  
print(lst)
sink()

#plot the point paths in list
#compare to plot(tr)
plot(NA, xlim=c(0,1), ylim=c(-2, 1.5), type="n")
lapply(lst, lines)

#The higher the better, relative value
loglike
tipdata

       