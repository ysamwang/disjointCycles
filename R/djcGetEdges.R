#' Estimate edges when given a topological ordering fo cyclic graphs with disjoint cycles
#'
#' @param topOrdering a list of lists containing the topological layers estimated by djcGetOrder
#' @param Y n x p matrix of observations
#' @param alpha the level for the hypothesis test when pruning parent edges
#' @param rescaleData whether to rescale data after each step
#' @param verbose print out results for each step
#' @return
#' \itemize{
#' \item adjMat the estimated adjacency matrix where adjMat[i,j] == 1 indicates the edge i -> j
#' \item Lambda the estimated edge weights where Lambda[i,j] indicates the linear coefficient of i onto j
#' \item pvals_adj the adjusted p-values for potential parents to be pruned. If the value is -1 this indicates
#'       that there was no test for the parent in this stage. This is because it is part of a cycle so the
#'       edge was forced to exist or not exist based on the greedy algorithm or because one does not
#'       preceed the other in the topological layers.
#' }
#'
djcGetEdges <- function(topOrdering, Y, alpha, verbose = F, rescaleData = T, pvalAdjMethod = "holm"){

  ## Iterates through the topological layers
  # and orients edges in a cycle and estimates edgeweights
  # for non-root layers we start with an edge from each node in a preceeding
  # layer to every edge in the current layer and then prune edges to estimate a graph
  # uses helper functions constructCycle and pruneParents


  p <- ncol(Y)

  if(rescaleData){
    Y <- scale(Y)
  }


  moments <- disjointCycles::calcSandT(Y)

  Lambda <- adjMat <- matrix(0, p, p)
  pvals <-matrix(-1, p, p)
  # Processing the first topological layer
  for(k in 1:length(topOrdering[[1]])){

    # if the length of the vector is greater than 1, then it corresponds to a cycle
    # so orient the cycle and estimate the edge weights
    if(length(topOrdering[[1]][[k]]) > 1){

      D <- topOrdering[[1]][[k]]
      out_D <- disjointCycles::constructCycle(moments$sMat[D, D, drop = F], moments$tMat[D, D, D, drop = F])
      Lambda[D, D] <- out_D$Lambda
      adjMat[D, D] <- out_D$adjMat

    }

  }


  # If there is only 1 layer in the topological ordering then everything is done
  # No pruning is necessary aside from orienting edges in cycles
  if(length(topOrdering) == 1){

    return(list(adjMat = adjMat, Lambda = Lambda))

  }


  # if there are layers remaining recur through those
  for(i in 2:length(topOrdering)){

    # C contains all nodes in preceeding layers
    # these are the set of potential parents for nodes in the current
    # topological layer
    C <- unlist(topOrdering[1:(i-1)])

    # iterates through units in the current topological layer
    for(k in 1:length(topOrdering[i])){

      # D contains a single set of units from the topological layer
      # can be either a singleton or a cycle
      D <- topOrdering[[i]][[k]]
      out_D <- disjointCycles::pruneParents(C, D, moments$sMat, moments$tMat, Y)

      # update matrices which will be returned with estimates
      Lambda[c(C,D), D] <- out_D$Lambda_CD
      adjMat[c(C,D), D] <- out_D$adj_CD
      pvals[c(C,D), D] <- out_D$pvals

    }

  }

  # Adjust all raw p-values for pruning a parent for multiple testing
  # p-values which correspond to -1 are not actually tested, so shouldn't be counted in the
  # multiplicity adjustment
  pvals[which(pvals > -1)] <- p.adjust(pvals[which(pvals > -1)], method = pvalAdjMethod)

  # remove support from Lambda and adjMat if adjusted p-value is larger
  adjMat <- adjMat * (pvals <= alpha)
  Lambda <- Lambda * (pvals <= alpha)

  return(list(adjMat = adjMat, Lambda = Lambda, pvals_adj = pvals))


}
