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

#' Compare whether two estimated orders are the same
#'
#' @param sMat the covariance matrix for the variables in a root cycle.
#'    If the cycle is not a root cycle, it should be the covariance of the variables after adjusting
#'    for nodes in preceeding layers
#' @param tMat the third moment tensor for the variables in a root cycle.
#'    If the cycle is not a root cycle, it should be the third moments of the variables after adjusting
#'    for nodes in preceeding layers
#' @return
#' \itemize{
#' \item adjMat the estimated adjacency matrix where adjMat[i,j] == 1 indicates the edge i -> j
#' \item Lambda the estimated edge weights where Lambda[i,j] indicates the linear coefficient of i onto j
#' }
#'
compareOrders <- function(X, Y){

  if(length(X) * length(Y) == 0){
    print("Length 0 Error")
    print(X)
    print(Y)
  }

  if(length(X) != length(Y)){
    return(FALSE)
  }

  all(unlist(mapply(function(x,y){mapply(setequal, x, y)}, X, Y)))
}


#' Compare whether two estimated orders are the same
#'
#' @param sMat the covariance matrix for the variables in a root cycle.
#'    If the cycle is not a root cycle, it should be the covariance of the variables after adjusting
#'    for nodes in preceeding layers
#' @param tMat the third moment tensor for the variables in a root cycle.
#'    If the cycle is not a root cycle, it should be the third moments of the variables after adjusting
#'    for nodes in preceeding layers
#' @return
#' \itemize{
#' \item adjMat the estimated adjacency matrix where adjMat[i,j] == 1 indicates the edge i -> j
#' \item Lambda the estimated edge weights where Lambda[i,j] indicates the linear coefficient of i onto j
#' }
#'
compareOrders <- function(X, Y){

  if(length(X) * length(Y) == 0){
    print("Length 0 Error")
    print(X)
    print(Y)
  }

  if(length(X) != length(Y)){
    return(FALSE)
  }

  all(unlist(mapply(function(x,y){mapply(setequal, x, y)}, X, Y)))
}


anySingletonRoots <- function(X){
  any(sapply(X, function(x){any(sapply(x, length) == 1)}))
}



