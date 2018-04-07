# ver 0.3.0 - make the acceptace rates vectors
#fit_model <- function(fossils, tr, tipdata, rt_value, model, ...) UseMethod("fit_model")

fit_model.default <- function(fossils, tr, tipdata, rt_value = mean(tipdata),
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
  mcmc_type = "tanner-wong", alpha = NULL, mu = NULL, sigma = NULL,
  N = 100, init_method = "sim", update_method = "subtree", iters = 5000,
  method = "euler", ...)
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

if (!ape::is.binary(tr)) {
  stop("The tree must be rooted.", call. = FALSE)
}
if (min(tr$edge.length) < 1 / N) {
  warning("No data will be imputed to at least ", sum(tr$edge.length < 1 / N),
  " edges in the tree.\nIncreasing N to avoid that.", call. = FALSE)
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
  lst <- phylo_sde_0(fossils = fossils, tr = tr, rt_value = rt_value, theta = theta, model = M,
                        N = N, method = method, ...)
} else {
  stop("The init_method specified is not available.", .call = FALSE)
}


loglike <- numeric(iters)
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

if (mcmc_type == "tanner-wong") {
  for (k in 2:iters) {
    setTxtProgressBar(pb, k)
   
    out_mcmc <- mcmc_steps_tanner_wong(fossils = fossils, tr = tr, tipdata = tipdata,
                rt_value = rt_value, lst = lst,
                theta = theta, model = M, para2est = para2est,
                update_method = update_method, proposals = proposals,
                priors = priors, method = method, N = N, ...)
    
    lst <- out_mcmc$lst
    theta <- out_mcmc$theta
    mcmctrace[k, para2est] <- theta[1, para2est]
    n_para_accept[k] <- out_mcmc$n_para_accept
    n_data_accept[k] <- out_mcmc$n_data_accept
  }
  close(pb)
} else {
  for (k in 2:iters) {
    setTxtProgressBar(pb, k)
    out_mcmc <- mcmc_steps_else(fossils = fossils, tr = tr, tipdata = tipdata,
                rt_value = rt_value, lst = lst,
                theta = theta, model = M, para2est = para2est,
                update_method = update_method, proposals = proposals,
                priors = priors, method = method, N = N, ...)
    lst <- out_mcmc$lst
    theta <- out_mcmc$theta
    mcmctrace[k, para2est] <- theta[1, para2est]
    n_para_accept[k] <- out_mcmc$n_para_accept
    n_data_accept[k] <- out_mcmc$n_data_accept
  }
  close(pb)
}

bt <- back_transform(model = model, tipdata = tipdata, rt_value = rt_value,
      lst = lst)
#tipdata <- bt$tipdata; rt_value <- bt$rt_value
lst <- bt$lst

attr(mcmctrace, "class") <- "mcmc"
attr(mcmctrace, "mcpar") <- c(Start = 1, End = iters, frequency = 1)
output <- list(mcmctrace = mcmctrace, lst = lst, loglike = loglike,
               n_data_accept = n_data_accept, n_para_accept = n_para_accept)
output$call <- match.call()
class(output) <- "menura"
return(output)
}

