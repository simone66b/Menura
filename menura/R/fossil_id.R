fossil_id <- function(ftr, tr){
  # identify which branches have length 0
  br_zero <-  which(ftr$edge.length[] == 0.0) 
  br_zero
  #from this set, need to identify which ones are connected to tips
  tip <- ftr$edge[br_zero,2]
  #provides a set of node values with edge 0, potential tips
  
  fossils <- 0
  #fossils will be added to the tree separately, this means they will have tip values greater than the original tips
  #and less than the internal notes, the difference in tips in the range they are in
  if(all.equal(ftr,tr)){
    fossils <- NULL
    return (fossils)
  }else{
      for(i in 1:length(tip)){
          if ((tip[i] <= Ntip(ftr)) && (tip[i] > Ntip(tr))){
             fossils[i] <- tip[i]
          }
      } 
      return (fossils)
  }
}  

