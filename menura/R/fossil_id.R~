fossil_id <- function(ftr){
    br_zero <- which(ftr$edge.length[] == 0)
    nodes <- ftr$edge[br_zero, 2]
    fossils <- 0
    for(i in 1:length(nodes)){
      if (nodes[i] <= ape::Ntip(ftr)){
         fossils[i] <- nodes[i]
      }
    } 
    return(fossils)
  }



