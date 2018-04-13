
### This is a test file running through all functions given a tree of type phylo
## A birfucating phylogenetic tree of type phylo is required to run this package
## Polytomies are auto resolved
## fossils must be bound to the tree as tips with a given tipvalue and an edgelength of 0.0, they are inserted
## at secific positions along a branch given a specific node
## 


#clears environment
rm(list=ls())

#set working directory 
setwd("/home/cheyennem/Documents/Menura/menura/R")


# initialize functions
source("fossil_id.R")
source("fit_model.R")
source("phylo_sde_0.R")
source("mcmc_steps_tanner_wong.R")
source("mcmc_steps_else.R")
source("tree_logL.R")
source("order_tree.R")
source("update_subtree.R")
source("update_tree.R")
source("sde_model.R")
source("igamma.R")
source("dc_fn.R")
source("back_transform.R")


set.seed(1)
rpkgs <- c("sde", "ape", "msm")
lapply(rpkgs, require, character.only = TRUE)
# Number of tips
#ntips <- 128
# true values of parameters used to generate tipdata
true.alpha <- 10
true.mu <- 15
true.sigma <- 5


t.root.value <- true.mu

# number of mcmc iterations
iters <- 500



#examples{
#set seed for random num generator to ensure simulation is reproducible
set.seed(1)
# install these packages
rpkgs <- c("sde", "ape", "msm")
lapply(rpkgs, require, character.only = TRUE)

# Random tree with n tips
tr <-  compute.brlen(rtree(n=128))

# plot the tree
plot(tr)
edgelabels()
nodelabels()

#manually adding fossils to tree as new branches with length 0 
#rt_node is shifted up by one after each addition, other nodes may change value depending on insertion point
#be cautious when adding multiple tips at once that they are inserted correctly

#create a tip of class phylo
tip <- list(edge = matrix(c(2,1),1,2), 
            tip.label = "fossil",
            edge.length = 0.0,
            Nnode = 1)
class(tip)<- "phylo"

#bind the tip to the tree at node x and position y
ftr1 <- bind.tree(tr,tip, where = 183, position = 0.1)
plot(ftr1)


tip <- list(edge = matrix(c(2,1),1,2), 
            tip.label = "fossil",
            edge.length = 0.0,
            Nnode = 1)
class(tip)<- "phylo"

ftr2 <- bind.tree(ftr1,tip, where = 164, position = 0.1)
plot(ftr2)

#set the new fossilized tree to the working tree variable tr
tr <- ftr2
#you should now see the same tree with additional fossils
plot(tr)
edgelabels()
nodelabels()


# Generate tip values
set.seed(1)
# tr <-  compute.brlen(stree(n=ntips, type="balanced"))


#originally x and l undefined, generates tipdata for CIR method
# f_TrCir <- function(x, l){
#   x <- 1
#   l <- 0.1
#   rcCIR(n=1, Dt=l, x0=x, theta=c(true.alpha*true.mu, true.alpha, true.sigma))
# }
# 
# t.tipdata <- rTraitCont(tr, f_TrCir, ancestor = FALSE, root.value = t.root.value)

##Generates tipdata for OU method
 t.tipdata <- rTraitCont(tr, "OU", sigma=true.sigma, alpha=true.alpha, theta=true.mu,
                      root.value=t.root.value)


## identify fossils within the tree
fossils <- fossil_id(tr)

## run the simulation and mcmc
set.seed(1)
model.1 <- fit_model.default(fossils = fossils, tr=tr, tipdata=t.tipdata, rt_value=t.root.value, iters=iters,
                     model = "OU", alpha = 10,  mu = 15, sigma = NULL,
                     N=240, init_method = "sim", update_method = "subtree")

# Look at the MCMC trace of the parameters
summary(model.1)
model.1

## Use coda package to analyse the mcmc chain
 require(coda)
 plot(model.1$mcmctrace)
 summary(model.1$mcmctrace)

## The same can be done in the following way
# model <- list()
# model$d <- function (t, x, theta) {
#   ((theta[1]*theta[2] - 0.25* theta[3]^2) / (2 * x)) - theta[2] * x / 2
# }
# model$s <- function(t, x, theta) {
#   0.5 * theta[3]
# }
# model$drift <- quote(((alpha * mu - 0.25 * sigma^2) / (2*x)) -
#                        alpha * x / 2)
# model$diffusion <- quote(sigma/2)
# model$dx_diffusion <- quote(0)
# 
# tipdata  <- sqrt(t.tipdata)
# rt_value <- sqrt(t.root.value)
# set.seed(1)
# model.2 <- fit_model.default(fossils = fossils, tr=tr, tipdata=tipdata, rt_value=rt_value, iters=iters,
#                      model = model, alpha = 10,  mu = 15, sigma = NULL,
#                      N=240, init_method = "sim", update_method = "subtree")
# 
# plot(model.2$mcmctrace)


