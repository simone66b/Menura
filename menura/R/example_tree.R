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
  
  
  ##working on manually adding fossils to tree as new branches with length 0 ######################################
  tip <- list(edge = matrix(c(2,1),1,2), 
              tip.label = "fossil",
              edge.length = 0.0,
              Nnode = 1)
  class(tip)<- "phylo"
  
  ftr1 <- bind.tree(tr,tip, where = 6, position = 0.1)
  plot(ftr1)
  nodelabels()
 
  tip <- list(edge = matrix(c(2,1),1,2), 
              tip.label = "fake",
              edge.length = 0.2,
              Nnode = 1)
  class(tip)<- "phylo"
  
  ftr <- bind.tree(ftr1,tip, where = 7, position = 0.1) 
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
  
  
########################################################################################################################  
  
  
  
##new attempt at locating fossils
  tr = ftr
  plot(tr)
  edgelabels()
  nodelabels()
  
  #which branches have a length of 0
  br_zero = which(tr$edge.length[] == 0.0)
  br_zero
  
  #which branches have tips attached
  br_classes = branchClasses(tr)
  br_classes
  
  num_tips <- length(tr$tip.label)
  
  #assigns dead tips to tr, represents any edge that is connected to a tip that is not at present
  tr$brlen.dead = br_classes$brlen.dead
  #this should work
  
  ###IDENTIFIES THE NODE # OF A FOSSIL
  #easy to identify which branches have length 0
    br_zero <-  which(tr$edge.length[] == 0.0) 
    br_zero
  #from this set, need to identify which ones have tipdata
    #edges are labeled according to their row in the edge matrix
    is_tip <- tr$edge[br_zero,2]
    is_tip
    #how many other rows have the same node value in column 2
    #a tip will only show up once
    tip_check <- which(tr$edge[,2] == is_tip)
    tip_check
    length(tip_check)
    fossil = 0
    #if there are more than 1 edges with the same column 2 node value, the edge does not belong to a fossil
    #if there is only 1 edge with that node value, it must be a fossil
    if (length(tip_check) == 1)
      fossil = is_tip
    
    fossil
    
    
 
 
  
  #edge length to test if == 0.0
  tr$edge.length
  Nedges(tr)
  tr$tip.label
  tr$Nnode
  
  
  
  
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

       