

### Code chooses MRCA of 2 random tips for fossils.

library(ape)
library(phytools)

tree.size <- 20
num.fossils <- 5


tree <- compute.brlen(rtree(tree.size, rooted=TRUE))
##tree <- compute.brlen(stree(tree.size, type="balanced"))
for (i in 1:num.fossils) {
    is.fossil <- TRUE
    while(is.fossil) {
        to.place  <- sample(tree$tip.label,2)
        if (length(grep("f", to.place)) > 0) {
            is.fossil <- TRUE
        } else {
            anc <- getMRCA(tree, to.place)
            desc <- tree$edge[which(tree$edge[,1] == anc),2]
            desc.tips <- tree$tip.label[desc]
            if (length(grep("f", desc.tips) > 0)) {
                is.fossil <- TRUE
            } else {
                is.fossil <- FALSE
            }
        }
    }
    tree <- bind.tip(tree, paste("f", i, sep=""), edge.length=0, where=anc)
}
plot(tree)



###################################################################3
## Code chooses random nodes for fossils.
library(ape)
library(menura)
library(phytools)

tree.size <- 20
num.fossils <- 5

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


traits <- rTraitCont(tree, model="OU", root.value=0)
fossil.data <- traits[fossil_id(tree)]

tree <- multi2di(tree)

fossil_id <- function(ftr){
    br_zero <- which(ftr$edge.length[] == 0)
    nodes <- ftr$edge[br_zero, 2]
    fossils <- 0
    for(i in 1:length(nodes)){
      if (nodes[i] <= ape::Ntip(ftr)){
         fossils[i] <- nodes[i]
      }
    } 
    return(na.omit(fossils))
  }


test  <- fit_model(tree, traits, rt_value=0, model="OU", alpha=1, mu=0,
                   sigma=NULL, N=1000, iters=100, fossils=fossil.data)

phylo_sde <-  function (tr, rt_value, N, theta, model, method, traitdata, fossils = NULL, ...) {
    lst <- list()
    n_tips <- length(tr$tip.label)
    rt_node <- n_tips + 1
    dotslist <- list(...)
    if ("pred.corr" %in% names(dotslist)) {
        pred.corr <- dotslist$pred.corr
    }
    else {
        pred.corr <- FALSE
    }
    if (method == "milstein") {
        pred.corr <- TRUE
    }
    if (pred.corr) {
        if (!exists("dx_diffusion", model)) {
            model$dx_diffusion <- D(model$diffusion, "x")
        }
    }
    else {
        model$dx_diffusion <- quote(NULL)
    }
    
    sde_edges <- function(fossils, tr, node, X0, t0, traits, fossil.data) {
        daughters <- tr$edge[which(tr$edge[, 1] == node), 2]
        if (length(daughters) == 0) return()
        if (any(daughters %in% fossils)) {
            edge <- which((tr$edge[, 1] == node) & !(tr$edge[, 2] %in% fossils))
            f_edge <- which((tr$edge[, 1] == node) & (tr$edge[, 2] %in% fossils))
           
            root <- tr$edge[edge, 2]

            drift <- as.expression(force(eval(substitute(substitute(e, 
                list(alpha = theta[edge, "alpha"], mu = theta[edge, 
                  "mu"], sigma = theta[edge, "sigma"])), list(e = model$drift)))))
            diffusion <- as.expression(force(eval(substitute(substitute(e, 
                list(alpha = theta[edge, "alpha"], mu = theta[edge, 
                  "mu"], sigma = theta[edge, "sigma"])), list(e = model$diffusion)))))
            diffusion_x <- as.expression(force(eval(substitute(substitute(e, 
                list(alpha = theta[edge, "alpha"], mu = theta[edge, 
                  "mu"], sigma = theta[edge, "sigma"])), list(e = model$dx_diffusion)))))
            n_steps <- tr$edge.length[edge] * N
            tE <- t0 + tr$edge.length[edge]
            if (tr$edge.length[edge] == 0) {
                lst[[edge]] <<-list(NULL)
            } else {
                lst[[edge]] <<- sde::sde.sim(X0 =traits[which(names(traits) %in% names(fossil.data))],
                                            t0 = t0, T = tE, 
                                            N = n_steps, method = method, drift = drift, 
                                            sigma = diffusion, sigma.x = diffusion_x,
                                            pred.corr = pred.corr)
            }
            tE <- tsp(lst[[edge]])[2]
            sde_edges(fossils, tr, root, lst[[edge]][n_steps + 1], t0=tE, traits, fossil.data)
        } else {
            for (d_ind in 1:2) {
                edge <- which((tr$edge[, 1] == node) & (tr$edge[, 2] == daughters[d_ind]))
                drift <- as.expression(force(eval(substitute(substitute(e, 
                  list(alpha = theta[edge, "alpha"], mu = theta[edge, 
                    "mu"], sigma = theta[edge, "sigma"])), list(e = model$drift)))))
                diffusion <- as.expression(force(eval(substitute(substitute(e, 
                  list(alpha = theta[edge, "alpha"], mu = theta[edge, 
                    "mu"], sigma = theta[edge, "sigma"])), list(e = model$diffusion)))))
                diffusion_x <- as.expression(force(eval(substitute(substitute(e, 
                  list(alpha = theta[edge, "alpha"], mu = theta[edge, 
                    "mu"], sigma = theta[edge, "sigma"])), list(e = model$dx_diffusion)))))
                n_steps <- tr$edge.length[edge] * N
                tE <- t0 + tr$edge.length[edge]

                if (tr$edge.length[edge] == 0) {
                    lst[[edge]] <<- list(NULL)
                } else {
                    lst[[edge]] <<- sde::sde.sim(X0 = X0, t0 = t0, 
                                                 T = tE, N = n_steps, method = method, drift = drift, 
                                                sigma = diffusion, sigma.x = diffusion_x,
                                                pred.corr = pred.corr)
                    tE <- tsp(lst[[edge]])[2]
                }
                if (lengths(lst[[edge]]) <= 1) t0 <- t0 else t0 <- lst[[edge]][n_steps+1]
                sde_edges(fossils, tr, daughters[d_ind], X0=X0,  t0=t0, traits=traits,
                          fossil.data=fossil.data)
                
            }
        }
     }   
    sde_edges(fossils, tr, rt_node, X0 = rt_value, t0 = 0, traits=traits, fossil.data=fossil.data)
    return(lst)

    }


    
  node_len <- ape::node.depth.edgelength(tr)
    for (nthtip in 1:n_tips) {
        nEdge <- which(tr$edge[, 2] == nthtip)
        ntsp <- tsp(lst[[nEdge]])
        if (nthtip %in% fossils) {
            ntsp[2] <- 0
            lst[[nEdge]] <- 0
        }
        if (ntsp[2] > (node_len[nthtip] - 1/N)) {
            if (length(lst[[nEdge]]) > 1) {
                if (length(lst[[nEdge]]) == 2) {
                  lst[[nEdge]] <- ts(lst[[nEdge]][1], start = ntsp[1], 
                    end = ntsp[1], frequency = N)
                }
                else {
                  lst[[nEdge]] <- window(lst[[nEdge]], start = ntsp[1], 
                    end = ntsp[2] - 1/N)
                }
            }
        }
    }
    return(lst)
}

lst <- lst2
 plot(NA, xlim=c(0,1), ylim=c(-2,2), type="n")
lst2 <- lst[lengths(lst) > 1]
lapply(lst2, lines)
X11();plot(tree)


alpha <- 1
mu <- 0
theta <- cbind(alpha=alpha, mu=mu, sigma=1)
rt_value <- 0
N <- 200

model="OU"
method="euler"
 
