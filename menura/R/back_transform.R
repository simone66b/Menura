back_transform <- function(model, rt_value, tipdata, lst, ...) {

  if (any(class(model) %in% "list"))
    return(list(tipdata = tipdata, rt_value = rt_value, lst = lst))

  if (model == "CIR") {
    tipdata  <- tipdata ^ 2
    rt_value <- rt_value ^ 2
    lst <- lapply(1:length(lst), function(i) lst[[i]] ^2)
  }

  if (model == "Beta") {
    tipdata  <-  sin(tipdata / 2)  #2 * asin(tipdata)
    rt_value <-  sin(rt_value / 2) #2 * asin(rt_value)
    lst <- lapply(1:length(lst), function(i)
                            unlist(lapply(lst[[i]]/2, FUN=sin)) )
  }

  return(list(tipdata = tipdata, rt_value = rt_value, lst = lst))
}
