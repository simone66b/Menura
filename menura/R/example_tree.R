#examples{
  #set seed for random num generator to ensure simulation is reproducible
  set.seed(1)
  # install these packages
  rpkgs <- c("sde", "ape", "msm")
  lapply(rpkgs, require, character.only = TRUE)

  # Random tree with n tips
  tr <-  compute.brlen(rtree(n=4))

  # plot the tree
  plot(tr)
  edgelabels()
  nodelabels()
  add.scale.bar()
  ftr <- tr
  ##manually adding fossils to tree as new branches with length 0 ######################################
 
  tip <- list(edge = matrix(c(2,1),1,2), 
              tip.label = "fossil",
              edge.length = 0.0,
              Nnode = 1)
  class(tip)<- "phylo"
  
  ftr1 <- bind.tree(tr,tip, where = 4, position = 0.1)
  plot(ftr1)
  ftr <- ftr1
  tip <- list(edge = matrix(c(2,1),1,2), 
              tip.label = "fossil",
              edge.length = 0.0,
              Nnode = 1)
  class(tip)<- "phylo"
  
  ftr1 <- bind.tree(ftr1,tip, where = 7, position = 0.1)
  plot(ftr1)
  
  tip <- list(edge = matrix(c(2,1),1,2), 
              tip.label = "fossil",
              edge.length = 0.0,
              Nnode = 1)
  class(tip)<- "phylo"
  
  ftr <- bind.tree(ftr1,tip, where = 11, position = 0.1) 
  #number of tips in the new tree
  ftr_tips <- Ntip(ftr)
  #plot the tree
  plot(ftr)
  edgelabels()
  
 nodelabels()
 ftr <- ftr1
  plot(ftr)
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
  


  #takes in the tree with fossils and the original tree
  #returns a vector of the nodes numbers that are fossils
  fossils <- fossil_id(ftr)
  fossils 
 
  
  #edge length to test if == 0.0
  tr$edge.length
  tr$tip.label
  tr$Nnode
  
  
  
  # SDE parameters
  # set to size of the length of the edge.length vector
  Nedges <- length(ftr$edge.length)  
  
  dclade <- max(which(ftr$edge[,1] == ftr$edge[1,1])) - 1  
  
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
  tipdata <- rTraitCont(ftr, "OU", sigma=sigma, alpha=alpha, theta=mu,
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
  lst <- phylo_sde_0 (fossils = fossils, tr=ftr, rt_value=rt_value, theta=theta, model=model,
                    N=N, method="euler")
  
  lst
  #calls log likelihood using Euler, approximates for diffusion process in tree
  loglike <-  tree_logL (fossils = fossils, tr=ftr, tipdata=tipdata, lst=lst,
                         alpha=theta[, "alpha"],
                         mu=theta[, "mu"],
                         sigma=theta[, "sigma"], model=model,
                         method = "euler")
  
  
  lst
  loglike
  #lst stores the point path to each node
 
 
  #print lst to a file
#sink("phylo_sde_test_output")  
#print(lst)
#sink()

#plot the point paths in list
#compare to plot(tr)
plot(NA, xlim=c(0,1), ylim=c(-2, 1.5), type="n")
lapply(lst, lines)

#The higher the better, relative value

tipdata

       