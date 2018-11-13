\name{tree_logL}
\alias{tree_logL}
\title{Calculates the Log Likelihood of the Diffusion Process in a Tree of Life}
\usage{
tree_logL(fossils, tr, tipdata, lst, alpha, mu, sigma, model, method, ...)
}
\arguments{
\item{fossils} {Name of nodes that are fossils.}
\item{tr}{object of the class "phylo" as in the \code{ape} R package which
contains the information of the tree structure.}

\item{tipdata}{a numeric vector containing tip values in the
same order of the tip labels in \code{tr$tip.label}.}

\item{lst}{list containing the diffusion paths of the \code{tr} object in the same
order as the edges in the \code{tr} object.}

\item{alpha}{vector containing the parameter value of theta_1 for each edge \code{tr} object,
and the parameters location corresponds to the same edge numbering in the \code{tr} object.}

\item{mu}{similar to definition of \code{alpha}.}

\item{sigma}{similar to definition of \code{alpha}.}

\item{model}{a list which contains functions \code{d}, \code{s} and, possibly, \code{s_x} which are drift component, diffusion component, and partial derivative of diffusion component of the diffusion process.
}

\item{method}{Numerical method apprximate the sde.}

\item{...}{not used.}
}
\value{
logL a number.
}
\description{
Euler approximated Log likelihood of the diffusion process in the tree.
}
\examples{
set.seed(1)
rpkgs <- c("sde", "ape", "msm")
lapply(rpkgs, require, character.only = TRUE)
# Number of tips
# Random tree with 64 tips
tr <-  compute.brlen(rtree(n=64))

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
lst <- phylo_sde (tr=tr, rt_value=rt_value, theta=theta, model=model,
                   N=N, method="euler")

loglike <-  tree_logL (tr=tr, tipdata=tipdata, lst=lst,
                              alpha=theta[, "alpha"],
                              mu=theta[, "mu"],
                              sigma=theta[, "sigma"], model,
                              method = "euler")

}
