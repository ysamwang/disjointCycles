#' Estimate topological ordering for cyclic graphs with disjoint cycles
#'
#' @param Y n x p matrix of observations
#' @param alpha2 the level for the hypothesis test when determining singleton roots
#' @param alpha3 the level for the hypothesis test when determining
#' @param alphaR the level for the hypothesis test when
#' @param methodRoot which method to use for testing singleton roots;
#'  i.e., for fixed u \eqn{H_0: det_{u,v}^{2x2} == 0} for v!=u. Options are indDelta which tests each value individually using
#'  asymptotic normality and the adjusting for multiple testing with holm or jointDelta which uses asymptotic normality of entire vector
#'   and does a single chi-squared test
#' @param methodDet2 which method to use for testing \eqn{H_0: det_{u,v}^{2x2} != 0} or \eqn{det_{v,u}^{2x2} != 0}
#' Options are indDelta which tests each value individually using
#'  asymptotic normality and the adjusting for multiple testing with holm or jointDelta which tests whether the product of the two quantities
#'  are non-zero and uses asymptotic normality
#' @param rescaleData whether to rescale data after each step
#' @param verbose print out results for each step
#' @return
#' A list of lists where each internal list contains vectors of nodes which are either of length 1 (indicating
#' a root singleton) or of length greater than 1 (indicating a root cycle)
#'
djcGetOrderNew <- function(Y, alpha2 = .01, alpha3 = .01, alphaR = .01,
                        methodRoot = "indDelta", methodDet2 = "indDelta",
                        rescaleData = T, verbose = F, pvalAdjMethod = "holm",
                        methodPR = "infFunc", sigmaPop = NULL, rescaleResid = T){

  # Uses djc_recur to iteratively estimate root layers

  # p is number of nodes
  # unordered is the nodes not yet placed in topological ordering
  # topOrder is a list of list where each interior list is a layer of topological orders
  p <- ncol(Y)
  unordered <- 1:p
  topOrder <- list()


  # First pass
  #
  roots <- disjointCycles::djc_recurNew(Y, alpha2 = alpha2, alpha3 = alpha3, alphaR = alphaR,
                                     methodRoot = methodRoot, methodDet2 = methodDet2,
                                     rescaleData = rescaleData,  pvalAdjMethod =  pvalAdjMethod,
                                     methodPR = methodPR, sigmaPop = sigmaPop, rescaleResid = rescaleResid)

  # First layer
  layer <- 1
  topOrder[[layer]] <- roots$roots

  # Remove nodes placed in first topological layer from unordered nodes
  unordered <- setdiff(unordered, unlist(roots$roots))

  if(verbose){

    cat("=== Layer: ")
    cat(layer)
    cat(" ===\n")
    cat("Roots: ")
    print(roots$roots)
    cat("\n")
  }


  # While unordered nodes still remain
  while(length(unordered) > 0){
    layer <- layer + 1


    # estimate next root layer
    roots <- disjointCycles::djc_recurNew(roots$newY, alpha2 = alpha2, alpha3 = alpha3,
                                       alphaR = alphaR,
                                       methodRoot = methodRoot,
                                       methodDet2 = methodDet2,
                                       rescaleData = rescaleData,
                                       pvalAdjMethod = pvalAdjMethod,
                                       methodPR = methodPR, sigmaPop = sigmaPop, rescaleResid= rescaleResid)

    # djc_recur returns indices corresponding to 1:ncol(newY)
    # indexBack maps those back to the original indices
    indexBack <- lapply(roots$roots, function(x){unordered[x]})

    # Put nodes in estimated layer into the topological ordering
    topOrder[[layer]] <- indexBack

    # remove nodes from estimated layer from unordered nodes
    unordered <- setdiff(unordered, unlist(indexBack))

    if(verbose){

      cat("=== Layer: ")
      cat(layer)
      cat(" ===\n")
      cat("Roots: ")
      print(indexBack)
      cat("\n")
    }

  }

  ## Order each layer (and cycles within layer) according
  # to original input ordering to ease comparison
  topOrder <- disjointCycles::orderList(topOrder)

  return(topOrder)


}


### Internal helper function which estimates a root layer
# Returns two objects
#   roots: vectors of nodes which are either of length 1 (indicating
#'        a root singleton) or of length greater than 1 (indicating a root cycle)
#'  newY: an n x p matrix which contains Y for the remaining unordered nodes
djc_recurNew <- function(Y, alpha2 = .01, alpha3 = .01, alphaR = .01,
                      methodRoot = "indDelta", methodDet2 = "indDelta", rescaleData = T,
                      pvalAdjMethod, methodPR, sigmaPop = NULL, rescaleResid = T){



  # p is the number of nodes
  # n is the number of observations
  p <- ncol(Y)
  n <- nrow(Y)


  # If only 1 node remains, it must be the root
  if(p == 1){
    return(list(roots = list(c(1)), newY = NULL))
  }

  # If data should be rescaled before estimating roots
  if(rescaleData){
    Y <- scale(Y)

    if(!is.null(sigmaPop)){
      sigmaPop <- cov2cor(sigmaPop)
    }

  }

  # Calculate 2nd and 3rd moments
  moments <- disjointCycles::calcSandT(Y)


  #### Get singleton Roots ####

  # get the nodes u for which d_{u,v} == 0 for all v
  # i.e., null is not rejected when is pVal >= alpha2
  # If methodRoot = indDelta, each test is done individually  (d_{u,v} == 0 for all v != u) then p-values are
  #     aggregated (for a single u across different v) using holm
  # if methodRoot = jointDelta, a single chi^2 test (for each u) is done
  # In both cases, the p-values (for each u) is adjusted for multiple testing
  # singletonRoots should be a vector


  # Will hold p-values for testing null hypothesis that
  # d_{uv}^{2x2} == 0
  # only fills in upper triangle
  # Use NA so p.adjust does not count diagonals

  det2 <- matrix(NA, p, p)

  for(u in 1:p){
    for(v in 1:p){

      # P-value for checking null that det_{u,v} == 0
      if(u != v){
        det2[u , v] <- .test2ind(u, v, moments$sMat, moments$tMat, Y)
      }
    }
  }

  # adjust for multiplicity
  det2 <- matrix(p.adjust(det2, method = pvalAdjMethod), p, p)

  if(methodRoot == "indDelta"){

    # singletonRoots contains any nodes which do not have rejected null
    singletonRoots <- which(apply(det2, MARGIN = 1, FUN = min, na.rm = T) > alpha2)


  } else if(methodRoot == "jointDelta") {

    # singletonRoots contains any nodes do not have rejected null
    singletonRoots <- which(p.adjust(
      sapply(1:p, function(u){.test2joint(u, moments$sMat, moments$tMat, Y, method = "jointDelta")$pval}),
      method = pvalAdjMethod)  > alpha2)



  } else {
    cat("Incorrect option for methodRoot!")
    stop()
  }



  # If there is at least 1 singleton root certified
  if(length(singletonRoots) > 0){

    # If there are any remaining unordered nodes, then adjust for certified roots
    if(length(singletonRoots) < p){

      # Regress out roots and return
      newY <- matrix(0, n, p)
      for(v in setdiff(1:p, singletonRoots)){

        # regress away all singleton roots
        newY[, v] <- RcppArmadillo::fastLm(Y[, v] ~ Y[, singletonRoots] - 1)$res
      }

      # If there are any singleton nodes, then return the roots and adjusted Y
      return(list(roots = as.list(singletonRoots), newY = newY[, -singletonRoots, drop = F]))

    }  else {

      ## All nodes are roots, so no need to adjust Y
      return(list(roots = as.list(singletonRoots), newY = NULL))
    }

  }


  #### If no singleton roots are found, then test for root cycles ####

  diag(det2) <- 1
  ## Get pairs such that d_{u,v} and d_{v,u} are simultaneously non-zero
  # If no pairs are jointly non-zero at alpha2, then raise threshold
  # to the minimum p-value in det2 so at least 1 pair is jointly certified
  threshold2 <- max(alpha2, min(pmax(det2, t(det2)), na.rm = T))
  nonZero2 <- (pmax(det2, t(det2)) <= threshold2)


  # Will hold p-values for testing if d_{uv}^{3x3} == 0
  # only fill in upper triangle
  # Use NA so p.adjust does not count untested pairs
  det3 <- matrix(NA, p, p)

  for(u in 1:(p-1)){
    for(v in (u+1):p){
      if(nonZero2[u, v] == 1){
        det3[u , v] <- .test3(u, v, moments$sMat, moments$tMat, Y)$pval
      }
    }
  }

  # adjust for multiple testing
  det3[upper.tri(det3, diag = F)] <- p.adjust(det3[upper.tri(det3, diag = F)], method = pvalAdjMethod)



  # Set the threshold to alpha3
  # or decrease to the threshold to the largest p-value if it is smaller
  # than alpha3 so that there is at least 1 pair connected nodes
  threshold3 <- min(alpha3, max(det3, na.rm = T), na.rm = T)

  # create graph which contains the edge u--v
  # if d_{uv}^{3x3} == 0 and d_{uv}^{2x2} != 0 and d_{vu}^{2x2} != 0
  # use upper since det3 is symmetric and we only use upper tri

  G <- igraph::graph_from_adjacency_matrix( det3 >= threshold3 , mode = "upper")

  # compute set of max cliques of at least size 2
  # there must be at least 1 clique since threshold3 is selected to have at least
  # 1 edge in the graph G
  # use as.vector to convert from igraph object to vector
  maxCliques <- lapply(igraph::max_cliques(G, min = 2), function(x){as.vector(x)})

  ### Pruning C to obtain C' ###
  if(length(maxCliques) == 1){

    # if there's only 1 maxclique, then set that as the root no need to prune
    rootCycles <- maxCliques

  } else {

    # if there is more than 1 maxClique found then prune
    # Note that not all maxCliques need be disjoint at this point



    # vector of all nodes in the max cliques list
    C1 <- unique(unlist(maxCliques))

    # will hold p-values for testing if clique k is a root
    pvalRoots <- rep(0, length(maxCliques))

    for (k in 1:length(maxCliques)){

      # Get clique k from maxCliques
      cliqueK <- maxCliques[[k]]

      # get other nodes in maxCliques list and get residuals
      C2 <- setdiff(C1, cliqueK)
      pvalRoots[k] <- pruneHelperRootCycle(cliqueK, C2, Y, methodPR = methodPR, rescaleResid = rescaleResid)
    }

    # Adjust for multiple testing
    pvalRoots <- p.adjust(pvalRoots, method = pvalAdjMethod)

    if(sum(pvalRoots >= alphaR) == 0 ){

      ## if all cycles are rejected, then set root cycle to be all combined
      ## i.e., all put into 1 big cycle
      rootCycles <- list(c(unique(unlist(maxCliques))))


    } else {

      rootCycles <- maxCliques[which(pvalRoots >= alphaR)]

    }





    # Check if any certified cycles are not disjoint
    # if so, combine any connected components together
    if(length(rootCycles) > 1){

      # Create graph with cliques
      temp_Adj <- matrix(0, p,p)
      for(k in 1:length(rootCycles)){
        temp_Adj[t(combn(sort(rootCycles[[k]]),2)) ] <- 1
      }
      G <- igraph::graph_from_adjacency_matrix(temp_Adj, mode = "upper")

      # Get connected components in graph G
      compOut <- igraph::components(G)
      # Set rootCycles to be each connected component
      rootCycles <- lapply(which(compOut$csize > 1), function(x){which(compOut$membership == x)})
    }


  }

  ## if the identified root cycles are all remaining nodes
  ## no need to adjust Y
  if(length(unlist(rootCycles)) == p){

    return(list(roots = rootCycles, newY = NULL))

  } else {





    if(!is.null(sigmaPop)){

      ordered <- unlist(rootCycles)
      unordered <- setdiff(1:p, ordered)
      newY <- Y[, unordered] - Y[, unlist(rootCycles)] %*% solve(sigmaPop[unlist(rootCycles), unlist(rootCycles)], sigmaPop[unlist(rootCycles), unordered])
      newSigmaPop <- sigmaPop[unordered, unordered, drop = F] - sigmaPop[unordered, ordered] %*% solve(sigmaPop[ordered, ordered]) %*% sigmaPop[ordered, unordered]

    } else {

      ordered <- unlist(rootCycles)
      unordered <- setdiff(1:p, ordered)

      newY <- matrix(0, n, p)
      for(v in 1:p){

        if(v %in% unordered){
          # Assumes centered Y so does not include intercept
          newY[, v] <- RcppArmadillo::fastLm(Y[, v] ~ Y[, ordered] - 1)$res
        }
      }

      newY <- newY[, -c(ordered), drop = F]
      newSigmaPop <- NULL
    }


    # rootCycles should be a list of vectors
    return(list(roots = rootCycles, newY = newY, sigmaPop = newSigmaPop))
  }

}





