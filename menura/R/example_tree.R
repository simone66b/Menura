#examples{
  #set seed for random num generator to ensure simulation is reproducible
  set.seed(1)
  # install these packages
  rpkgs <- c("sde", "ape", "msm")
  lapply(rpkgs, require, character.only = TRUE)
  # Number of tips
  # Random tree with n tips
  tr <-  compute.brlen(rtree(n=3))
  
  # plot the tree
  plot(tr)
  
  # SDE parameters
  # set to size of the length of the edge.length vector
  Nedges <- length(tr$edge.length)  
  
  dclade <- max(which(tr$edge[,1] == tr$edge[1,1])) - 1  
  
  #these parameters are vectors with size = number of edges = number of simulations
  alpha <- mu <- sigma <- rep(0, Nedges)
  
  #numeric vector of selective constraint strength for ec. branch
  alpha[1:Nedges]  <- 0.1
  
  #mu == theta 
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

loglike


       