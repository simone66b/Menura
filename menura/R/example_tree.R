#examples{
  #set seed for random num generator to ensure simulation is reproducible
  set.seed(1)
  # install these packages
  rpkgs <- c("sde", "ape", "msm")
  lapply(rpkgs, require, character.only = TRUE)
  # Number of tips
  # Random tree with 64 tips
  tr <-  compute.brlen(rtree(n=20))
  
  plot(tr)
  # SDE parameters
  Nedges <- length(tr$edge.length)
  dclade <- max(which(tr$edge[,1] == tr$edge[1,1])) - 1
  alpha <- mu <- sigma <- rep(0, Nedges)
  alpha[1:Nedges]  <- 0.1
  mu[1:Nedges] <- 0
  sigma[1:Nedges] <- 1
  rt_value <- 0
  tipdata <- rTraitCont(tr, "OU", sigma=sigma, alpha=alpha, theta=mu,
                        root.value=rt_value)
  model <- list()
  model$d <- function (t, x, theta) {
    theta[1] * (theta[2] - x)
  }
  
  model$s <- function(t, x, theta) {
    theta[3]
  }
  model$drift <- quote(alpha * (mu - x))
  model$diffusion <- quote(sigma)
  model$dx_diffusion <- quote(0)
  theta <- cbind(alpha=alpha, mu=mu, sigma=sigma)
  N <- 100
  lst <- phylo_sde_0 (tr=tr, rt_value=rt_value, theta=theta, model=model,
                    N=N, method="euler")
}
  
  lst
sink("phylo_sde_test_output")  
print(lst)
sink()

plot(NA, xlim=c(-1,1), ylim=c(-2, 1), type="n")
lapply(lst, lines)
