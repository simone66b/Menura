

### Code chooses MRCA of 2 random tips for fossils.

## library(ape)
## library(phytools)

## tree.size <- 20
## num.fossils <- 5


## tree <- compute.brlen(rtree(tree.size, rooted=TRUE))
## ##tree <- compute.brlen(stree(tree.size, type="balanced"))
## for (i in 1:num.fossils) {
##     is.fossil <- TRUE
##     while(is.fossil) {
##         to.place  <- sample(tree$tip.label,2)
##         if (length(grep("f", to.place)) > 0) {
##             is.fossil <- TRUE
##         } else {
##             anc <- getMRCA(tree, to.place)
##             desc <- tree$edge[which(tree$edge[,1] == anc),2]
##             desc.tips <- tree$tip.label[desc]
##             if (length(grep("f", desc.tips) > 0)) {
##                 is.fossil <- TRUE
##             } else {
##                 is.fossil <- FALSE
##             }
##         }
##     }
##     tree <- bind.tip(tree, paste("f", i, sep=""), edge.length=0, where=anc)
## }
## plot(tree)



###################################################################3
## Code chooses random nodes for fossils.
library(ape)
## library(menura)
library(phytools)

tree.size <- 10
num.fossils <- 4

tree <- compute.brlen(rtree(tree.size))
## tree <- compute.brlen(stree(tree.size, type="balanced"))
for (i in 1:num.fossils) {
    is.fossil <- TRUE
    while(is.fossil) {
        anc  <- sample(unique(tree$edge[,1]), 1)
        desc <- tree$edge[which(tree$edge[,1] == anc),2]
        desc.tips <- tree$tip.label[desc]
        if (length(grep("f", desc.tips) > 0)) {
            is.fossil <- TRUE
        } else {
            is.fossil <- FALSE
        }
    }
    tree <- bind.tip(tree, paste("f", i, sep=""), edge.length=0, where=anc)
}


plot(tree)
nodelabels()

tr <- multi2di(tree)

traits <- rTraitCont(tr, model="OU", root.value=0, ancestor=TRUE)
fossil.data <- traits[fossil_id(tree)]



fossil_id <- function(tr){
    br_zero <- which(tr$edge.length == 0)
    if (length(br_zero) == 0) return(NULL)
    nodes <- tr$edge[br_zero, 2]
    fossils <- 0
    for(i in 1:length(nodes)) {
      if (nodes[i] <= ape::Ntip(ftr)){
         fossils[i] <-  nodes[i]
      }
    } 
    return(nodes)
  }


test  <- fit_model(tree, rt_value=0, model="OU", alpha=1, mu=0,
                   sigma=NULL, N=100, iters=100, fossils=fossil.data)

set.seed(1)
fossils <- grep("f", tr$tip.label)

sde_edges(tr=tr, node=15, X0 = rt_value, t0 = 0, traits=traits)

sde_edges <- function(tr, node, X0, t0, traits) {
    daughters <- tr$edge[which(tr$edge[, 1] == node), 2]

            for (d_ind in 1:2) {
                edge <- which((tr$edge[, 1] == node) & (tr$edge[, 2] == daughters[d_ind]))
                
                brlen <- tr$edge.length[edge]
                ## print(brlen)
                ## print(edge)
                drift <- as.expression(force(eval(substitute(substitute(e, 
                  list(alpha = theta[edge, "alpha"], mu = theta[edge, 
                    "mu"], sigma = theta[edge, "sigma"])), list(e = model$drift)))))
                diffusion <- as.expression(force(eval(substitute(substitute(e, 
                  list(alpha = theta[edge, "alpha"], mu = theta[edge, 
                    "mu"], sigma = theta[edge, "sigma"])), list(e = model$diffusion)))))
                diffusion_x <- as.expression(force(eval(substitute(substitute(e, 
                  list(alpha = theta[edge, "alpha"], mu = theta[edge, 
                    "mu"], sigma = theta[edge, "sigma"])), list(e = model$dx_diffusion)))))
                n_steps <- brlen * N
                tE <- t0 + brlen

                if (brlen == 0) { ##&& tr$edge[edge, 2] <= Ntip(tr)) { ## fossil edge
                    lst[[edge]] <<- ts(traits[tr$edge[edge,2]], start=t0, end=t0)
                     X0 <- lst[[edge]][n_steps + 1]
                }
               

                if (brlen > 0) {
                    if (t0 == 0) X0 <- 0
                        lst[[edge]] <<- sde::sde.sim(X0 = X0, t0 = t0, 
                        T = tE, N = n_steps, method = method, drift = drift, 
                        sigma = diffusion, sigma.x = diffusion_x, pred.corr = pred.corr)
                        tE <- tsp(lst[[edge]])[2]
                        X0 <- lst[[edge]][n_steps + 1]
                }
                if (daughters[d_ind] > Ntip(tr)) {
                    sde_edges(tr, daughters[d_ind], X0= X0, tE, traits)
                }
            }
}
          


    
