#'fit_model
#'
#'Bayesian Estimator of Parameters of Univariate Diffusion Model for Continuous Trait Evolution
#'
#'  This function estimates posterior distributions for evolutionary
#'  models of continuous traits in a phylogeny. The evolutionary
#'  processes considered here belong to a class of diffusion processes
#'  which are typically given as solutions to the stochastic
#'  differential equations of the form given by
#'  \deqn{dX_t = a(X_t, alpha, mu) dt + b(X_t, sigma) dW_t,~~X_0 = x_0}
#'  where \eqn{X_t} denotes the state variable, \code{t} the time,
#'  \code{a} the drift function and \code{b}
#'  the diffusion function are known
#'  in parametric from which
#'  \code{alpha}, \code{mu}, and \code{sigma}
#'  are the parameters, and \eqn{W_t} is Brownian motion.
#'  The value of \eqn{X_t} at time \eqn{t_0}, \eqn{X_0},
#'  is independent of \eqn{W_t}.
#'
#'  Given root and tip values, the tree,
#'  and drift and diffusion functions,
#'  Markov Chain Monte Carlo (MCMC) estimates
#'  of the parameters are obtained.
#'  The parameter to be estimated is assumed to be the same
#'  for all the branches of the tree; however,
#'  the rest can be allowed to vary.

#'  Due to the low frequency nature of the data, the MCMC
#'  method used assumes intermediate data
#'  are missing. Thus, during each step of the MCMC iterations,
#'  missing data and the model parameters are imputed.
#'
#'  If the diffusion process is
#'  Ornstein-Uhlenbeck (OU),
#'  Cox-Ingersoll-Ross (CIR), or Beta (Beta)
#'  this can be specified by setting the
#'  \code{model} to "OU", "CIR", or
#'  "Beta", respectively.
#'  If the diffusion process is not one of these, the estimation of model parameters
#'  can be done by specifying drift and diffusion coefficients
#'  in a list assigned to \code{model}. In this case,
#'  the list object \code{model} must include : 1)
#'  \code{d} and \code{s} which are
#'  functions of time (\code{t}), space (\code{x}) and a vector of the parameter values (\code{theta})
#'  consisting of \code{alpha}, \code{mu} and \code{sigma},
#'  2) drift coefficients as a list containing
#'  a \code{quote} and 3) a
#'  diffusion coefficient as a list object
#'  which includes \code{diffusion} as
#'  the diffusion coefficients and
#'  \code{x} as the first derivatives of
#'  diffusion coefficients.
#' 
#'  The Ornstein-Uhlenbeck (OU) model is given by
#'  \deqn{dX_t = alpha (mu - X_t) dt + sigma~ dW_t,}
#'  with \eqn{X_0 = x_0 > 0}, \eqn{W_t} is the Brownian process,
#'  alpha, mu, and sigma are the model parameters
#'  where alpha and sigma are positive values.
#'
#'  The Cox-Ingersoll-Ross (CIR) model is given by
#'  \deqn{dXt = alpha (mu - X_t) dt + sigma~sqrt(X_t) dWt,}
#'  with \eqn{X_0 = x_0 > 0}, where Wt is the Brownian process, alpha, mu,
#'  and sigma are the model parameters which are all
#'  positive values. If the \code{model == "CIR"} is specified,
#'  then the parameter estimation is done by
#'  using transformation \eqn{Y = sqrt(X)} of Ito diffusion
#'  process.
#'
#'  The Beta model is given by
#'  \deqn{ dX_t = alpha(mu - X_t) dt + sigma~sqrt(X_t (1-X_t)) dW_t,}
#'  with \eqn{X_0 = x_0 > 0}, \eqn{W_t} is the Brownian process,
#'  alpha, mu, and sigma are the model parameters
#'  where alpha and sigma are positive values.
#'  If the \code{model ==} "Beta" is specified,
#'  then the parameter estimation is done by
#'  using transformation \eqn{Y = 2~sin^{-1}(X)} of Ito diffusion
#'  process.
#'
#'@param tr Single evolutionary tree as an object of the "phylo" class in
#'the \code{ape} package.
#'@param tipdata A numeric vector containing tip values in the same order as the tip labels in \code{tr$tip.label}.
#'@param rt_value Value at the root of \code{tr}.
#'@param model Either a list containing drift and diffusion coefficients
#'in quote format as functions of alpha, mu and sigma,
#'or a string ("OU", "CIR", or "Beta") specifying the diffusion process.
#'See Details.
#'@param priors A list of lists containing functions for prior distributions
#'of the model parameters. Use of the default priors is not recommended.
#'@param proposals A list of lists containing functions for proposal distributions
#'of the model parameters.
#'@param mcmc_type Type of MCMC algorithm, "DA" or "Fuchs"
#'@param alpha Set to NULL if alpha is to be estimated,
#'otherwise set to a numeric value or a numeric vector specifying
#'the value of the parameter for all the branches/edges.
#'In the latter case, the values must be specifying in the same order
#'as the edges in the \code{tr} object.
#'@param mu As \code{alpha}.
#'@param sigma As \code{alpha}.
#'@param N Data augmentation frequency.
#'@param init_method Method for initial data imputation.
#'Currently only the \code{"sim"} option is available.
#'@param update_method Method for data imputation during the MCMC.
#'Option \code{"subtree"} will only update part of the tree at each iteration,
#'where as \code{"tree"} will update the whole tree. See Details.
#'@param iters Number of MCMC iterations.
#'@param method Numerical approximation method, "euler" or "milstein."
#'@param fossils A numeric vector containing the tip values for every fossil added to the tree.
#'@seealso \code{\link{fossil_id}}
#'@param ... Not used.
#'
#'@export
#'@importFrom graphics plot
#'@importFrom ape dist.nodes
#'@importFrom ape compute.brlen
#'@importFrom ape Ntip
#'@importFrom sde EULERloglik
#'@importFrom sde dcEuler
#'@importFrom sde dcElerian
#'@importFrom sde sde.sim
#'@importFrom stats dlnorm
#'@importFrom stats dnorm
#'@importFrom stats dunif
#'@importFrom stats rgamma
#'@importFrom stats rlnorm
#'@importFrom stats rnorm
#'@importFrom stats runif
#'@importFrom stats ts
#'@importFrom stats tsp
#'@importFrom stats window
#'@importFrom stats reorder
#'@importFrom utils setTxtProgressBar txtProgressBar
#'@importFrom stats D time

fit_model <- function(tr, tipdata, rt_value = mean(tipdata),
  model = "OU",
  priors = list(
    alpha = list(df =  function(x, a = 1, b = 125, log_scale = TRUE) {
                            dunif(x, min = a, max = b, log = log_scale)},
                rf = function(n, a = 1, b = 125) {
                            runif(n, min = a, max = b)}),
    mu = list(df = function(x, a = 0, b = 20, log_scale = TRUE) {
                            dnorm(x, mean = a, sd = b, log = log_scale)},
              rf = function(n, a = 0, b = 20) {
                            rnorm(n, mean = a, sd = b)}),
    sigma = list(df = function(x, a = 1, b = 225, log_scale = TRUE) {
                            dunif(x, min = a, max = b, log = log_scale) },
                 rf = function(n, a = 1, b = 225) {
                            runif(n, min = a, max = b)})
  ),
  proposals = list(
    alpha = list(df = function(n, alpha, gamma = 0.5, log_scale = TRUE) {
                        dlnorm(n, meanlog = log(alpha), sdlog = gamma,
                               log = log_scale)},
                rf = function(n, alpha, gamma = 0.5) {
                       rlnorm(n, meanlog = log(alpha), sdlog = gamma) }),
    mu = list(df = function(n, mu, gamma = 0.5, log_scale = TRUE) {
                     dnorm(n, mean = mu, sd = gamma, log = log_scale)},
              rf = function(n, mu, gamma = 0.5) {
                     rnorm(n, mean = mu, sd = gamma)}),
    sigma = list(df = function(n, sigma, gamma = 0.5, log_scale = TRUE) {
                        dlnorm(n, meanlog = log(sigma), sdlog = gamma,
                               log = log_scale)},
                 rf = function(n, sigma, gamma=0.5) {
                        rlnorm(n, meanlog = log(sigma), sdlog = gamma)})
  ),
  mcmc_type = "DA", alpha = NULL, mu = NULL, sigma = NULL,
  N = 1000, init_method = "sim", update_method = "subtree", iters = 5000,
  method = "euler", fossils = NULL, ...)
{

if (class(tr) != "phylo")
  stop("The object tr must be an object of class",
       " \"phylo\" as in the ape R package.", call. = FALSE)

if (!is.vector(tipdata))
  stop("tipdata must be a vector.", call. = FALSE)

if (!is.numeric(tipdata))
  stop("tipdata must be a numeric vector.", call. = FALSE)

if (length(tipdata) != length(tr$tip.label))
  stop("Presently, tipdata must contain all",
       "and only tips values.", call. = FALSE)

if (any(is.na(tipdata[1:length(tr$tip.label)])))
  stop("All tip values must be given.", call. = FALSE)

if ((length(rt_value) > 1) || is.null(rt_value) || !is.numeric(rt_value))
  stop("Root must have a numeric value.", call. = FALSE)

if (!((class(model) == "list") || (class(model) == "character")))
  stop("Specification of the model is not correct. See documentation.",
  call. = FALSE)

if (class(model) == "character") {
  if ((model != "OU") && (model != "CIR") && (model != "Beta"))
    stop("The model must be OU, CIR, Beta, or list of parameters.",
    call. = FALSE)
}

n_para2est <- sum(is.null(alpha) + is.null(mu) + is.null(sigma))
if (n_para2est == 0)
  stop("At least one of alpha, mu or sigma needs to be estimated",
  call. = FALSE)

if ((length(priors) < n_para2est) || !(class(priors) == "list"))
  stop("priors must contain prior distribution",
       " functions for all parameters to be estimated.", call. = FALSE)

if ((length(proposals) < n_para2est) || !(class(proposals) == "list"))
  stop("proposals must be a list of lists containing proposal",
    " distribution functions for parameters to be estimated.", call. = FALSE)

if (!missing(alpha) && !is.null(alpha))
  if ((!is.numeric(alpha)) || is.na(alpha))
    stop("alpha must either be NULL or a numeric value.", call. = FALSE)

if (!missing(mu) && !is.null(mu))
  if ((!is.numeric(mu)) || is.na(mu))
    stop("mu must either be NULL or a numeric value.", call. = FALSE)

if (!missing(sigma) && !is.null(sigma))
  if ((!is.numeric(sigma)) || is.na(sigma))
    stop("Sigma must either be NULL or a numeric value.", call. = FALSE)

if (!is.numeric(N))
  stop("N must be a numeric", call. = FALSE)

if (N <= 1)
  stop("N must be greater than one.", call. = FALSE)

if (!is.numeric(iters))
  stop("iters must be a numeric", call. = FALSE)

if (!((update_method == "tree") || (update_method == "subtree")))
  stop("update_method must only be tree or subtree", call. = FALSE)

if(!((mcmc_type == "DA") || (mcmc_type == "Fuchs")))
  stop("mcmc_type must be only DA or Fuchs", call. = FALSE)

if (!ape::is.binary(tr)) {
  stop("The tree must be rooted.", call. = FALSE)
}
if (min(tr$edge.length) < 1 / N) {
  warning("No data will be imputed to at least ", sum(tr$edge.length < 1 / N),
  " edges in the tree including fossil edges.\nIncrease N to avoid that if these are not fossils.", call. = FALSE)
}
# Current program is only written with no root value updating.
# if ((!updateRoot) && is.null(rt_value))
#  stop ("If root value is not to estimate, it needs to be given")

dots <- list(...)
# if ('constraints' %in% names(dots))
# see what can be done about them
# ideally, priors shoud be specified in a way the that
# the constraint in model parameters are always statisfied

# set the tree to cladewise and reorder parameters

cdwise <- order_tree(tr, alpha, mu, sigma, priors=priors)
tr <- cdwise$tr

theta <- cdwise$theta
para2est <- cdwise$para2est
# set drift and diffusion expressions
sde_comp <- sde_model(model, rt_value, tipdata, ...)
M <- sde_comp$M
tipdata <- sde_comp$tipdata
rt_value <- sde_comp$rt_value

# initialize the path
n_edges <- length(tr$edge.length)
n_tips  <- length(tr$tip.label)
node_len <- ape::node.depth.edgelength(tr)

if (init_method == "sim") {
  lst <- phylo_sde(fossils = fossils, tr = tr, rt_value = rt_value, theta = theta, model = M,
                        N = N, method = method, ...)
} else {
  stop("The init_method specified is not available.", .call = FALSE)
}


    loglike <- numeric(iters)
    loglike[1] <- 0
mcmctrace <- matrix(NA, nrow = iters, ncol = length(para2est))
colnames(mcmctrace) <- para2est
mcmctrace[1, para2est] <- theta[1, para2est]
# for (cvar in para2est) {
#   mcmctrace[1, cvar] <- theta[1, cvar]
# }
theta_star <- theta
# MCMC iterations
n_para_accept <- n_data_accept <- numeric(iters)
n_para_accept[1] <- n_data_accept[1] <- 1
pb <- txtProgressBar(min=1, max=iters)

if (mcmc_type == "DA") {
  for (k in 2:iters) {
    setTxtProgressBar(pb, k)
   
    out_mcmc <- mcmc_steps_DA(fossils = fossils, tr = tr, tipdata = tipdata,
                rt_value = rt_value, lst = lst,
                theta = theta, model = M, para2est = para2est,
                update_method = update_method, proposals = proposals,
                priors = priors, method = method, N = N, ...)
    
    lst <- out_mcmc$lst
    loglike[k] <- out_mcmc$loglike
    theta <- out_mcmc$theta
    mcmctrace[k, para2est] <- theta[1, para2est]
    n_para_accept[k] <- out_mcmc$n_para_accept
    n_data_accept[k] <- out_mcmc$n_data_accept
    ## loglike[k] <- tree_logL(tr = tr, tipdata = tipdata, lst = lst,
    ##                         alpha = theta[, "alpha"], mu = theta[, "mu"],
    ##                         sigma = theta[, "sigma"],
    ##                         model = M,
    ##                         method = method)
  }
  close(pb)
} else if (mcmc_type == "Fuchs"){
  for (k in 2:iters) {
    setTxtProgressBar(pb, k)
    out_mcmc <- mcmc_steps_fuchs(fossils = fossils, tr = tr, tipdata = tipdata,
                rt_value = rt_value, lst = lst,
                theta = theta, model = M, para2est = para2est,
                update_method = update_method, proposals = proposals,
                priors = priors, method = method, N = N, ...)
    lst <- out_mcmc$lst
    loglike[k] <- out_mcmc$loglike
    theta <- out_mcmc$theta
    mcmctrace[k, para2est] <- theta[1, para2est]
    n_para_accept[k] <- out_mcmc$n_para_accept
    n_data_accept[k] <- out_mcmc$n_data_accept
    ## loglike[k] <- tree_logL(tr = tr, tipdata = tipdata, lst = lst,
    ##                         alpha = theta[, "alpha"], mu = theta[, "mu"],
    ##                         sigma = theta[, "sigma"],
    ##                         model = M,
    ##                         method = method)
  }
  close(pb)
}

bt <- back_transform(model = model, tipdata = tipdata, rt_value = rt_value,
      lst = lst)

lst <- bt$lst

attr(mcmctrace, "class") <- "mcmc"
attr(mcmctrace, "mcpar") <- c(Start = 1, End = iters, frequency = 1)
output <- list(mcmctrace = mcmctrace, lst = lst, loglike = loglike,
               n_data_accept = n_data_accept, n_para_accept = n_para_accept)
output$call <- match.call()
output$model_type <- c(model, method, mcmc_type)
output$model_info <- c("Model:" , model , "Method:" , method ,  "Mcmc Type:" , mcmc_type)
class(output) <- c(class(output), "menura")
return(output)
}

