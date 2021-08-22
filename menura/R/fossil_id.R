#'fossil_id
#'
#'Create a vector of the location of the fossils in a tree
#'
#'This function locates the tip values of fossils in the tree.  
#'It is primarily used in other functions of menura to handle 
#'fossils data but can be useful to locate the fossils in the tree.
#' 
#'@param ftr Single evolutionary tree as an object of the "phylo" class in
#'the \code{ape} package.
#'
#'@export
#'
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



