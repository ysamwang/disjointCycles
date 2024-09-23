#' Estimate edges in a cycle
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
#' }
#'
constructCycle <- function(sMat, tMat){

  p <- ncol(sMat)

  # If cycle is of length 2, then use formula to estimate edges
  # select set of edges which minimize product of the edgweights
  if(p == 2){

    # corresponding to the two "orientations"
    Lambda1 <- Lambda2 <-  matrix(0, 2, 2)

    u <- 1
    v <- 2

    a <- sMat[u, u] * tMat[u, u, v] - sMat[u,v] * tMat[u,u,u]
    b <- sMat[v, v] * tMat[u,u,u] - sMat[u,u] * tMat[u, v, v]
    c <- sMat[u, v] * tMat[u, v, v] - sMat[v,v] * tMat[u, u, v]

    # solution 1 to quadratic formula
    Lambda1[u, v] <- (-b + sqrt(b^2 - 4 * a * c)) / (2  * a)
    Lambda1[v, u] <- (tMat[u, u, v] - Lambda1[u, v] * tMat[u, u, u]) / (tMat[u, v, v] -  Lambda1[u, v] * tMat[u, u, v])

    # solution 2 to quadratic formula
    Lambda2[u, v] <- (-b - sqrt(b^2 - 4 * a * c)) / (2  * a)
    Lambda2[v, u] <- (tMat[u, u, v] - Lambda2[u, v] * tMat[u, u, u]) / (tMat[u, v, v] -  Lambda2[u, v] * tMat[u, u, v])

    # pick the orientation that minimizes the product of the edges
    if(abs(Lambda2[u, v] * Lambda2[v, u]) < abs(Lambda1[u, v] * Lambda1[v, u])){
      Lambda <- Lambda2
    } else {
      Lambda <- Lambda1
    }

    return(list(adjMat = matrix(c(0, 1, 1, 0), 2, 2), Lambda = Lambda))
  }

  # If the cycle is longer than 2 nodes, then use a greedy procedure to add edges
  # Greedily add an edge between the pair with the largest element of the inverse correlation
  # unless it causes some node to have degree greater than 2 or creates a chord in final cycle

  # The (u,v) element of the inverse of the correlation matrix is
  # proportional to det(corY[-u, -v])
  weights <- abs(solve(cov2cor(sMat)))
  diag(weights) <- 0

  # adjacency matrix of single disjont cycle being constructed
  adjMat <- matrix(0, p,p)

  # reachableNodes[i, j] == 1 indicates that i is connected to j
  # in the current adjaceny matrix
  reachableNodes <- diag(p)



  for(i in 1:(p-1)){

    ### Get Valid Edges ###
    # a potential edge is valid if both adjacent vertices have degree less than 2
    # and adding the edge would not create a cycle
    # i.e., make two already connected nodes adjacent

    # Edges for which adjacent nodes have degree less than 2
    unsaturatedNodes <- matrix(0, p,p)
    unsaturatedNodes[which(colSums(adjMat) < 2), which(colSums(adjMat) < 2)] <- 1

    # pair of nodes for which both nodes have degree less than 2
    # and are not currently connected
    # Note that reachableNodes has 1 on diagonal so validEdges is always 0 on diag
    validEdges <- unsaturatedNodes * !reachableNodes

    # Make valid edges upper triangular so that which.max below only returns 1 value
    validEdges[lower.tri(validEdges)] <- 0

    # Find valid edge with largest weight
    selectedPair <- arrayInd(which.max(weights * validEdges), .dim = dim(validEdges))

    # add edge to adjacency matrix
    adjMat[selectedPair[1], selectedPair[2]] <-
      adjMat[selectedPair[2], selectedPair[1]] <- 1


    # if (i,j) is the selected edge, then any node previously connected to i or j
    # is now connected to union of nodes previously connected to i
    # and nodes previously connected to j
    newReachable <- union(which(reachableNodes[selectedPair[1], ] == 1),
                        which(reachableNodes[selectedPair[2], ] == 1))
    reachableNodes[newReachable, newReachable] <- 1
  }

  # After p-1 edges have been added, we can close the cycle with
  # the only remaining unsaturated nodes
  unsaturatedNodes <- matrix(0, p,p)
  unsaturatedNodes[which(colSums(adjMat) < 2), which(colSums(adjMat) < 2)] <- 1

  # pair of nodes for which both nodes have degree less than 2
  # and are not currently connected
  validEdges <- unsaturatedNodes
  diag(validEdges) <- 0
  selectedPair <- arrayInd(which.max(weights * validEdges), .dim = dim(validEdges))

  adjMat[selectedPair[1], selectedPair[2]] <-
    adjMat[selectedPair[2], selectedPair[1]] <- 1

  # adjacency matrix now has a single undirected cycle
  # now we orient the cycle in an arbitrary direction



  dirAdj <- matrix(0, p, p)
  currentNode <- 1
  nextNode <- min(which(adjMat[currentNode, ] == 1))
  dirAdj[currentNode, nextNode] <- 1
  previousNode <- currentNode
  currentNode <- nextNode

  for(i in 1:(p-1)){
    # pick node adjacent to current node which is not the previous node
    nextNode <- setdiff(which(adjMat[currentNode, ] == 1), previousNode)
    dirAdj[currentNode, nextNode] <- 1

    previousNode <- currentNode
    currentNode <- nextNode

  }


  Lambda <- matrix(0, p ,p)

  for(u in 1:p){

      v <- which(dirAdj[u, ] == 1)
      w <- which(dirAdj[v, ] == 1)

      p <- sMat[u,u] * (tMat[u,v,w]^2 - tMat[u, u, w] * tMat[v, v, w]) +
        sMat[u,w] * (tMat[u,u,u] * tMat[v,v,w] - tMat[u,u,v] * tMat[u,v,w]) +
        sMat[v,w] * (tMat[u,u,v] * tMat[u,u,w] - tMat[u,u,u] * tMat[u,v,w])

      q <- -sMat[u,v] * (tMat[u,v,w]^2 - tMat[u,u,w] * tMat[v,v,w]) +
        sMat[u,w] * (tMat[u, v, v] * tMat[u,v,w] - tMat[u,u,v] * tMat[v,v,w]) +
        sMat[v,w] * (tMat[u,u,v] * tMat[u,v,w] - tMat[u,u,w] * tMat[u,v,v])

      Lambda[u,v] <- q / p

  }


  return(list(adjMat = dirAdj, Lambda = -Lambda))
}



pruneParents <- function(C, D, sMat, tMat, Y){


  # initializing objects to return
  adj_CD <- matrix(0, length(C) + length(D), length(D))
  Lambda <- matrix(0, length(C) + length(D), length(D))
  pvals_CD <- matrix(0, length(C) + length(D), length(D))


  # if D is a singleton and not a cycle, only need to prune parents
  if(length(D) == 1){

    # Use least squares coefficients
    Lambda_CD <- solve(sMat[C,C]) %*% sMat[C,D]

    pvals_CD[1:length(C), 1] <- disjointCycles::pruneHelperSingle(C, D, Y, Lambda_CD)


    Lambda[1:length(C), 1:length(D)] <- Lambda_CD
    adj_CD[1:length(C), 1] <- 1

    # no self edge so make (d, d) element of lambda and adjMat 0
    # and set p-value to -1
    Lambda[length(C)+ 1, 1] <- 0
    adj_CD[length(C)+ 1, 1] <- 0
    pvals_CD[length(C) + 1, 1] <- -1


    return(list(adj_CD = adj_CD, Lambda_CD = Lambda, pvals = pvals_CD))

  } else {
    # if length of D is greater than 1, then it is a cycle

    # Step 1: Orient the cycle and estimate internal edges

    # Regression coefficients of C onto D
    R_DC <- solve(sMat[C,C]) %*% sMat[C, D]
    Y_D_adjusted <- Y[, D] - Y[, C] %*% t(R_DC)

    # calculate moments of adjusted data
    adjMoments <- disjointCycles::calcSandT(Y_D_adjusted)
    # orient edges and estimate edgeweights using adjusted data
    Lambda_DD <- disjointCycles::constructCycle(adjMoments$sMat, adjMoments$tMat)

    # estimate edges from C into D
    Lambda_CD <- R_DC %*% (diag(length(D)) - Lambda_DD$Lambda)

    # for dInd in D, get p-values for edge of parent into dInd
    for(dInd in 1:length(D)){

      pvals_CD[1:length(C), dInd] <- disjointCycles::pruneHelper(C, D, dInd, Y, Lambda_CD, Lambda_DD$Lambda)

    }

    # Fill in non-pruned edgeweights and edges
    Lambda[1:length(C), 1:length(D)] <- Lambda_CD
    adj_CD[1:length(C), 1:length(D)] <- 1

    # fill in edgeweights and edges estimated from cycle
    Lambda[length(C)+ 1:length(D), 1:length(D)] <- Lambda_DD$Lambda
    adj_CD[length(C) + 1:length(D), 1:length(D)] <- Lambda_DD$adj

    # Set p-values for (D,D) elements to -1 since they aren't being pruned
    pvals_CD[length(C)+ 1:length(D), 1:length(D)]

    return(list(adj_CD = adj_CD, Lambda_CD = Lambda, pvals = pvals_CD))
  }


}


pruneHelperCycle <- function(C, D, dInd, Y, lambda_CD, lambda_DD){

  # parent of D[dInd]
  d2 <- which(lambda_DD[, dInd] != 0)

  pvals <- rep(0, length(C))

  for(cTest in 1:length(C)){

    # Null hypothesis has Lambda[C, dInd] with 0 for cTest
    lambda_Cd_null <- lambda_CD[, dInd]
    lambda_Cd_null[cTest] <- 0

    # Form errors under null hypothesis
    errsD <- Y[, D[dInd]] - Y[, C] %*% lambda_Cd_null - Y[, D] %*% lambda_DD[, dInd]


    # construct A matrix
    A <- matrix(0, length(C) + 1, length(C) + 1)

    # C_mod simply reorders C so that cTest appears first
    C_mod <- c(cTest, setdiff(C, cTest))
    for(c1 in 1:length(C)){
      for(c2 in 1:length(C)){
        A[c1, c2] <- mean(-Y[ ,C_mod[c2]] * Y[ ,C_mod[c1]]^2)
      }
      A[c1, length(C) + 1] <- mean(-Y[ ,D[d2]] * Y[ ,C_mod[c1]]^2)
    }


    for(c2 in 1:length(C)){
      A[length(C) + 1, c2] <- mean(-2 * errsD * Y[, C_mod[1]] * Y[, C_mod[c2]])
    }
    A[length(C) + 1, length(C) + 1] <- mean(-2 * errsD * Y[, C_mod[1]] * Y[, D[d2]])


    m <- cbind(Y[,C_mod]^2 * c(errsD), errsD^2 * Y[,C_mod[1]])

    g <- m[, 1 , drop = F] - t(A[1, -1, drop = F] %*% solve(A[-1, -1, drop = F]) %*% t(m[, -1, drop = F]))
    pvals[cTest] <- emplik::el.test(g, mu = 0)$Pval

  }
  return(pvals)

}


pruneHelperSingle <- function(C, D, Y, lambda_CD){

  pvals <- matrix(0, length(C))

  for(cTest in 1:length(C)){

    # Null hypothesis has 0 for cTest
    lambda_Cd_null <- lambda_CD
    lambda_Cd_null[cTest] <- 0

    # Form errors under null hypothesis
    errsD <- Y[, D] - Y[, C] %*% lambda_Cd_null


    # construct A matrix
    A <- matrix(0, length(C), length(C))

    # C_mod simply reorders C so that cTest appears first
    C_mod <- c(cTest, setdiff(C, cTest))

    for(c1 in 1:length(C)){
      for(c2 in 1:length(C)){
        A[c1, c2] <- mean(-Y[ ,C_mod[c2]] * Y[ ,C_mod[c1]]^2)
      }

    }

    # form estimating equations m
    m <- Y[, C_mod]^2 * c(errsD)

    # marginalize to only include first
    g <- m[, 1 , drop = F] - t(A[1, -1, drop = F] %*% solve(A[-1, -1, drop = F]) %*% t(m[, -1, drop = F]))
    pvals[cTest] <- emplik::el.test(g, mu = 0)$Pval

  }
  return(pvals)

}


