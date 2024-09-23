orderList <- function(X){

    # Sort each inner layer
    X <- lapply(X, function(Z){lapply(Z, sort)})


    .sortWithinLayer <- function(Z){
      # get min element of each vector in list
      minElem <- sapply(Z, function(Y){min(Y)})
      return(Z[order(minElem)])
    }

    X <- lapply(X, .sortWithinLayer)

    return(X)
}


compareOrders <- function(X, Y){
  if(length(X) != length(Y)){
    return(FALSE)
  }
  all(unlist(mapply(function(x,y){mapply(setequal, x, y)}, X, Y)))
}

L <- list(list(c(2), c(3), c(4) ), list(c(5,6,7)), list(c(8)), list(c(1)))
orderList(L)
