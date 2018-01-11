\name{phylo_sde}
\alias{phylo_sde}
\title{Simulate a CIR Diffusion Process in the Tree of Life}
\usage{
phylo_sde(tr, rt_value, N, theta, model, method = "euler", ...)
}
\arguments{
\item{tr}{a modified object of class \code{phylo} as in the
\code{ape} R package. In this version, the CIR process
parameters alpha, mu and sigma for each of the branch
is included as vectors in the same order
as the edge (branch) labelling.}

\item{rt_value}{value at the root.}

\item{N}{data imputation frequency.}

\item{theta}{matrix of parameter values for each edge of the tree.}

\item{model}{a list containing drift, diffusion
and their partial differentiation as quotes. For the Euler scheme
the drift coefficient as \code{drift}, the diffusion coefficient as
\code{diffusion}, and the partial differentiation of
\code{diffusion} by \code{x} as \code{dx_diffusion} is required.
See the Examples.}

\item{method}{currently only the "euler" scheme is used.}

\item{...}{not used.}
}
\value{
lst a list of length equal to the number of
branches in the object \code{tr}. Elements of lst are time series
objects which are the simulated paths.
}
\description{
Starting from the root of the tree - which is assumed to start at
time 0, and the root value is known - a recursive scheme is used in
simulating a CIR process in the branches of the tree.
}
\details{
The number of samples imputed in a branch (edge) is proportional to the length
of the branch. First, the samples are imputed for the two root edges.
The end points of these are taken to be the starting points of
the successive branches. This process is done recursive for until
the tip nodes are reached.

The number of samples imputed on a branch is equal to
\code{round(N * branch_length)}.
As, the next branch starts from the end point of the previous,
the branch start and stop times would change,
which would depend on the data imputation frequency N.

In case, if length of an edge is small, no samples may be imputed for
such a branch. As such, the simulated output may contain branches
of zero length. This can be avoided by employing higher value for N.
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
}
