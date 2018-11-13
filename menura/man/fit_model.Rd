\name{fit_model}
\alias{fit_model}
\alias{fit_model.default}
\title{Bayesian Estimator of Parameters of Univariate Diffusion Models for
Continuous Trait Evolution}
\usage{
fit_model(tr, tipdata, rt_value, model, ...)

\method{fit_model}{default}(tr, tipdata, rt_value = mean(tipdata),
  model = "OU",
  priors = list(
    alpha = list (df =  function(x, a = 1, b = 125, log_scale = TRUE) {
                            dunif(x, min = a, max = b, log = log_scale)},
                  rf = function(n, a = 1, b = 125) {
                            runif(n, min = a, max = b)} ),
    mu = list (df = function(x, a = 0, b = 20, log_scale = TRUE) {
                            dnorm(x, mean = a, sd = b, log = log_scale)},
               rf = function(n, a = 0, b = 20) {
                            rnorm(n, mean = a, sd = b)}),
    sigma = list (df = function(x, a = 1, b = 225, log_scale = TRUE) {
                            dunif(x, min = a, max = b, log = log_scale) },
                   rf = function(n, a = 1, b = 225) {
                            runif(n, min = a, max = b)})
  ),
  proposals = list(
    alpha = list (df = function(n, alpha, gamma = 0.5, log_scale = TRUE) {
                    dlnorm(n, meanlog = log(alpha), sdlog = gamma,
                    log = log_scale) },
                  rf = function(n, alpha, gamma = 0.5) {
                    rlnorm(n, meanlog = log(alpha), sdlog = gamma) }),
    mu = list (df = function(n, mu, gamma = 0.5, log_scale = TRUE) {
                    dnorm(n, mean = mu, sd = gamma, log = log_scale)},
               rf = function(n, mu, gamma = 0.5) {
                    rnorm(n, mean = mu, sd = gamma)}),
    sigma = list(df = function(n, sigma, gamma = 0.5, log_scale = TRUE) {
                    dlnorm(n, meanlog = log(sigma), sdlog = gamma,
                    log = log_scale)},
                  rf = function(n, sigma, gamma = 0.5) {
                    rlnorm(n, meanlog = log(sigma), sdlog = gamma)})
  ),
  mcmc_type = "tanner-wong", alpha = NULL, mu = NULL, sigma = NULL,
  N=100, init_method="sim", update_method="subtree", iters=5000,
  method = "euler", ...)
}
\arguments{
\item{tr}{single evolutionary tree as an object of the 'phylo' class in
the \code{ape} R package.}

\item{tipdata}{a numeric vector containing values of the trait at the tip.
These must be in the same order as those in the tr$tip.label}

\item{rt_value}{value of the trait at the root.}

\item{model}{either a list containing drift and diffusion coefficients
in quote format as functions of alpha, mu and sigma,
or a string ("OU", "CIR", or "Beta") specifying diffusion process.
See Details.}

\item{priors}{list of lists containing functions for prior distributions
of the model parameters.}

\item{proposals}{list of lists containing functions for proposal distributions
of the model parameters.}

\item{mcmc_type}{Type of MCMC algorithm}

\item{alpha}{NULL if alpha is to be estimated, else either a numeric
  value of a numeric vector specifying the value of the parameter for
  all the branches/edges. In the latter case, the values must be
  specifying in order of the edges in the \code{tr} object.}

\item{mu}{same as \code{alpha}}

\item{sigma}{same as \code{alpha}}

\item{N}{data augmentation frequency.}

\item{init_method}{method for initial data imputation.
Currently only the \code{"sim"} option is available.}

\item{update_method}{method for data imputation during the MCMC.
The \code{"subtree"} will only update a random part of the tree at each
iteration, where as the option \code{"tree"} will update the whole tree. See Details}

\item{iters}{number of MCMC iterations.}

\item{method}{Numerical approximation method to use.}

\item{...}{further arguments for future extensions.}
}
\description{
This function estimates posterior distributions for evolutionary models of continuous traits on a phylogeny. The evolutionary processes considered here belong to a class of diffusion processes which are typically given as solutions to the stochastic differential equations of the form given by
\deqn{dX_t = a(X_t, alpha, mu) dt + b(X_t, sigma) dW_t,~~X_0 = x_0}
where \eqn{X_t} denote the state variable (ie the trait), \code{t} the time,
the drift function \code{a} and
the diffusion function \code{b} are known
in parametric from where
\code{alpha}, \code{mu}, and \code{sigma}
are the parameters, and \eqn{W_t} is Brownian motion.
The value of \eqn{X_t} at time \eqn{t_0}, \eqn{X_0},
is independent of the \eqn{W_t}.
}
\details{
Given the root and tip values, the tree, drift and diffusion functions, the Data Augmentation - Markov Chain Monte Carlo (DA-MCMC) estimates of the parameters are obtained. Parameters may be the same across the tree or allowed to differ in different parts of the tree.

Due to the low frequency nature of the data, we employ Data Augmentation of the evolutionary trajectory (the 'fossil record'), effectively imputing the missing trajectory, and updating trajectories at each MCMC iteration.

If the diffusion process considered is the Ornstein-Uhlenbeck (OU) process, Cox-Ingersoll-Ross (CIR) process or Beta process (Beta) then this can be specified by setting the \code{model} to "OU", "CIR", or "Beta", respectively. In the case of a user-defined diffusion process, the estimation of model parameters can be done by specifying drift and diffusion coefficients in a list assigned to \code{model}. In this case, the list object \code{model} must include functions \code{d} and \code{s} which are functions of \code{t}, \code{x} and \code{theta} which are the time variable, the space variable, and a vector consisting of the parameter values \code{alpha}, \code{mu} and \code{sigma}, and drift coefficient as a list containing a \code{quote} and diffusion coefficient as a list objects which include \code{diffusion} as the diffusion coefficient and \code{x} as the first derivatives of diffusion coefficients.

The Ornstein–Uhlenbeck (OU) model is given by \deqn{dX_t = alpha (mu - X_t) dt + sigma~ dW_t,} with \eqn{X_0 = x_0 > 0}, \eqn{W_t} is Brownian motion, alpha, mu, and sigma are the model parameters where alpha and sigma are positive values.

The Cox-Ingersoll-Ross (CIR) model is given by \deqn{dXt = alpha (mu -
  X_t) dt + sigma~sqrt(X_t) dWt,} with \eqn{X_0 = x_0 > 0}, where Wt is
  Brownian motion, alpha, mu, and sigma are the model parameters which
  are all positive values. If the \code{model == "CIR"} is specified, then the parameter estimation is done by using the transformation \eqn{Y = sqrt(X)} of the Ito diffusion process.

The Beta model is given by \deqn{ dX_t = alpha(mu - X_t) dt + sigma~sqrt(X_t (1-X_t)) dW_t,} with \eqn{X_0 = x_0 > 0}, \eqn{W_t} is Brownian motion, alpha, mu, and sigma are the model parameters where alpha and sigma are positive values. If the \code{model ==} "Beta" is specified, then the parameter estimation is done by using transformation \eqn{Y = 2~sin^{-1}(X)} of the Ito diffusion process.
}

\section{Methods (by class)}{
\itemize{
\item \code{default}: Bayesian Estimator of Diffusion Process
}}
\examples{
set.seed(1)
rpkgs <- c("sde", "ape", "msm")
lapply(rpkgs, require, character.only = TRUE)
# Number of tips
ntips <- 128
# SDE parameters
true.alpha <- 10
true.mu <- 5
true.sigma <- 2
t.root.value <- true.mu
iters <- 200
# Generate tip values
set.seed(1)
tr <-  compute.brlen(stree(n=ntips, type="balanced"))
f_TrCir <- function(x, l)
  rcCIR(n=1, Dt=l, x0=x, theta=c(true.alpha*true.mu, true.alpha, true.sigma))
t.tipdata <- rTraitCont(tr, f_TrCir, ancestor = FALSE, root.value = t.root.value)

set.seed(1)
model.1 <- fit_model(tr=tr, tipdata=t.tipdata, rt_value=t.root.value, iters=iters,
                model = "CIR", alpha = 10,  mu = 15, sigma = NULL,
                N=10, init_method = "sim", update_method = "subtree")

# Look at the MCMC trace of the parameters
# summary(model.1)
model.1

## Use the coda package to analyse the mcmc chain
# require(coda)
# plot(model.1$mcmctrace)
# summary(model.1$mcmctrace)

## The same can be done using the user-defined specification:
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
model.2 <- fit_model(tr=tr, tipdata=tipdata, rt_value=rt_value, iters=iters,
                model = model, alpha = 10,  mu = 15, sigma = NULL,
                N=10, init_method = "sim", update_method = "subtree")

}
