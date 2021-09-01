
library(ape)
library(menura)

## Number of tips on the tree
ntips <- 128

# SDE parameters
true.alpha <- 3
true.mu <- 5
true.sigma <- 2

# Set a starting root value
t.root.value <- true.mu

# Iterations set to 50 to limit run time.  Most simulations require 200000 or more
iters <- 2000

# Create a random tree
tr <-  compute.brlen(stree(n=ntips, type="balanced"))
## tr <- compute.brlen(rtree(1000))
# Plot the tree to see the structure
plot(tr)
edgelabels()
nodelabels()


## Uncomment the tipdata generation method for your preferred model

# Generates tipdata for CIR method
#f_TrCir <- function(x, l){
#  x <- 1
#  l <- 0.1
#  rcCIR(n=1, Dt=l, x0=x, theta=c(true.alpha*true.mu, true.alpha, true.sigma))
#}

#t.tipdata <- rTraitCont(tr, f_TrCir, ancestor = FALSE, root.value = t.root.value)

# Generates tipdata for OU method
t.tipdata <- rTraitCont(tr, "OU", sigma=true.sigma, alpha=true.alpha, theta=true.mu,
                     root.value=t.root.value)
                     
                     ##<insert example 1 or 2>##

###Example 1: Run Menura

## Run the simulation and parameter estimation 
set.seed(12345678)
model.1 <- fit_model(tr=tr, tipdata=t.tipdata, rt_value=t.root.value, iters=iters,
                     model = "OU", alpha = 3,  mu = 5, sigma = NULL,
                     N=240, init_method = "sim", update_method = "subtree")


## Use coda package to analyse the mcmc chain
 require(coda)
 plot(model.1$mcmctrace)
summary(model.1$mcmctrace)

tr <- as.mcmc(model.1$mcmctrace)

codamenu() ## navigate the menus to get plots and statistical summaries

######################################################################

## Use phylo_sde or tree_logL separately
set.seed(1)

  # Random tree with n tips
  tr <-  compute.brlen(rtree(n=4))

  # plot the tree
  plot(tr)
  edgelabels()
  nodelabels()
  add.scale.bar()

  tr <- tree
  # SDE parameters
  Nedges <- length(tr$edge.length)  
  
  dclade <- max(which(tr$edge[,1] == tr$edge[1,1])) - 1  
  
  # These parameters are vectors with size = number of edges 
  alpha <- mu <- sigma <- rep(0, Nedges)
  
  # Numeric vector of selective constraint strength for ec. branch
  alpha[1:Nedges]  <- 0.1
  
  #represents numeric vector giving ec. branch optimum
  mu[1:Nedges] <- 0
  
  #numeric vector Std-dev of random component for ec. branch
  sigma[1:Nedges] <- 1
  
  rt_value <- 0
  
  #simulate tipdata
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
  model$dx_diffusion <- quote(0)

  #theta is a vector of the parameters
  theta <- cbind(alpha=alpha, mu=mu, sigma=sigma)
  N <- 100

#### USE THIS TO TEST
  lst2 <- phylo_sde(tr=tree, rt_value=rt_value, theta=theta, model=model,
                    N=1000, method="euler")

fossil.data <- 
  # lst stores the time series point path to each node
 
  # Calls log likelihood using Euler, approximates for diffusion process in tree
  loglike <-  tree_logL(tr=tree, tipdata=traits, lst=lst2,
                        alpha=theta[, "alpha"],
                        mu=theta[, "mu"],
                        sigma=theta[, "sigma"], model=model,
                        method = "euler", fossils=fossil.data)
  

 
lst3 <- lst2[]

#compare to plot(tr)
plot(tr)
                                        #plot the point paths in list
lst3 <- lst2[which(lengths(lst2) > 1)]
plot(NA, xlim=c(0,1), ylim=c(-2, 1.5), type="n")
lapply(lst3, lines)

X11loglike

