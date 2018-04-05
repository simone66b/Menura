
set.seed(1)
rpkgs <- c("sde", "ape", "msm")
lapply(rpkgs, require, character.only = TRUE)
# Number of tips
ntips <- 128
# SDE parameters
true.alpha <- 10
true.mu <- 5
true.sigma <- 2
#why
t.root.value <- true.mu
iters <- 500
# Generate tip values
set.seed(1)
tr <-  compute.brlen(stree(n=ntips, type="balanced"))
tr <- ftr
plot(tr)
f_TrCir <- function(x, l){
  x <- 1
  l <- 0.1
  rcCIR(n=1, Dt=l, x0=x, theta=c(true.alpha*true.mu, true.alpha, true.sigma))
}

t.tipdata <- rTraitCont(tr, f_TrCir, ancestor = FALSE, root.value = t.root.value)
#t.tipdata <- rTraitCont(ftr, "OU", sigma=sigma, alpha=alpha, theta=mu,
                     # root.value=rt_value)


t.tipdata
#t.tipdata[6] <- 0.1836
fossils <- fossil_id(tr)
fossils
set.seed(1)
model.1 <- fit_model.default(fossils = fossils, tr=tr, tipdata=t.tipdata, rt_value=t.root.value, iters=iters,
                     model = "CIR", alpha = 10,  mu = 15, sigma = NULL,
                     N=240, init_method = "sim", update_method = "subtree")

# Look at the MCMC trace of the parameters
summary(model.1)
model.1

## Use coda package to analyse the mcmc chain
 require(coda)
 plot(model.1$mcmctrace)
 summary(model.1$mcmctrace)

## The same can be done in the following way
model <- list()
model$d <- function (t, x, theta) {
  ((theta[1]*theta[2] - 0.25* theta[3]^2) / (2 * x)) - theta[2] * x / 2
}
model$s <- function(t, x, theta) {
  0.5 * theta[3]
}
model$drift <- quote(((alpha * mu - 0.25 * sigma^2) / (2*x)) -
                       alpha * x / 2)
model$diffusion <- quote(sigma/2)
model$dx_diffusion <- quote(0)

tipdata  <- sqrt(t.tipdata)
rt_value <- sqrt(t.root.value)
set.seed(1)
model.2 <- fit_model.default(fossils = fossils, tr=tr, tipdata=tipdata, rt_value=rt_value, iters=iters,
                     model = model, alpha = 10,  mu = 15, sigma = NULL,
                     N=240, init_method = "sim", update_method = "subtree")

plot(model.2$mcmctrace)
}

