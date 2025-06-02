#' Estimate edges in a root cycle
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
  # now we orient the cycle in the direction with smallest cycle weight




  ### Try first direction ###
  dirAdj1 <- matrix(0, p, p)
  Lambda1 <- matrix(0, p ,p)

  currentNode <- 1
  # note that this takes min as opposed to below, where max is used
  nextNode <- min(which(adjMat[currentNode, ] == 1))
  dirAdj1[currentNode, nextNode] <- 1
  previousNode <- currentNode
  currentNode <- nextNode

  for(i in 1:(p-1)){
    # pick node adjacent to current node which is not the previous node
    nextNode <- setdiff(which(adjMat[currentNode, ] == 1), previousNode)
    dirAdj1[currentNode, nextNode] <- 1

    previousNode <- currentNode
    currentNode <- nextNode
  }



  for(u in 1:p){

      v <- which(dirAdj1[u, ] == 1)
      w <- which(dirAdj1[v, ] == 1)

      pTerm <- sMat[u,u] * (tMat[u,v,w]^2 - tMat[u, u, w] * tMat[v, v, w]) +
        sMat[u,w] * (tMat[u,u,u] * tMat[v,v,w] - tMat[u,u,v] * tMat[u,v,w]) +
        sMat[v,w] * (tMat[u,u,v] * tMat[u,u,w] - tMat[u,u,u] * tMat[u,v,w])

      qTerm <- -sMat[u,v] * (tMat[u,v,w]^2 - tMat[u,u,w] * tMat[v,v,w]) +
        sMat[u,w] * (tMat[u, v, v] * tMat[u,v,w] - tMat[u,u,v] * tMat[v,v,w]) +
        sMat[v,w] * (tMat[u,u,v] * tMat[u,v,w] - tMat[u,u,w] * tMat[u,v,v])

      Lambda1[u,v] <- qTerm / pTerm

  }

  ### Try Other Direction ###
  dirAdj2 <- matrix(0, p, p)
  Lambda2 <- matrix(0, p ,p)

  currentNode <- 1
  # note that this takes max as opposed to before, we took min
  nextNode <- max(which(adjMat[currentNode, ] == 1))
  dirAdj2[currentNode, nextNode] <- 1
  previousNode <- currentNode
  currentNode <- nextNode

  for(i in 1:(p-1)){
    # pick node adjacent to current node which is not the previous node
    nextNode <- setdiff(which(adjMat[currentNode, ] == 1), previousNode)
    dirAdj2[currentNode, nextNode] <- 1

    previousNode <- currentNode
    currentNode <- nextNode

  }



  for(u in 1:p){

    v <- which(dirAdj2[u, ] == 1)
    w <- which(dirAdj2[v, ] == 1)


    pTerm <- sMat[u,u] * (tMat[u,v,w]^2 - tMat[u, u, w] * tMat[v, v, w]) +
      sMat[u,w] * (tMat[u,u,u] * tMat[v,v,w] - tMat[u,u,v] * tMat[u,v,w]) +
      sMat[v,w] * (tMat[u,u,v] * tMat[u,u,w] - tMat[u,u,u] * tMat[u,v,w])

    qTerm <- -sMat[u,v] * (tMat[u,v,w]^2 - tMat[u,u,w] * tMat[v,v,w]) +
      sMat[u,w] * (tMat[u, v, v] * tMat[u,v,w] - tMat[u,u,v] * tMat[v,v,w]) +
      sMat[v,w] * (tMat[u,u,v] * tMat[u,v,w] - tMat[u,u,w] * tMat[u,v,v])

    Lambda2[u,v] <- qTerm / pTerm

  }


  ### Choose direction with smaller cycle weight
  if(abs(prod(Lambda1[which(Lambda1 != 0)])) <   abs(prod(Lambda2[which(Lambda2 != 0)]))){

    Lambda <- Lambda1
    dirAdj <- dirAdj1

  } else {

    Lambda <- Lambda2
    dirAdj <- dirAdj2

  }



  return(list(adjMat = dirAdj, Lambda = -Lambda))
}


#' Test existence of edges from nodes in preceeding layer
#'
#' @param C the set of nodes in the preceeding layer
#' @param D the set of nodes in the current unit (can be a singleton or cycle)
#' @return
#' \itemize{
#' \item adjMat the estimated adjacency matrix where adjMat[i,j] == 1 indicates the edge i -> j
#' \item Lambda the estimated edge weights where Lambda[i,j] indicates the linear coefficient of i onto j
#' \item pvals the pvalues (unadjusted for multiple testing) for any edge to be tested. If the value is -1 this indicates
#'       that there was no test for the parent in this stage. This is because it is part of a cycle so the
#'       edge was forced to exist or not exist based on the greedy algorithm or because one does not
#'       preceed the other in the topological layers.
#' }
#'
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
    Y_D_adjusted <- Y[, D] - Y[, C] %*% R_DC

    # calculate moments of adjusted data
    adjMoments <- disjointCycles::calcSandT(Y_D_adjusted)
    # orient edges and estimate edgeweights using adjusted data
    Lambda_DD <- disjointCycles::constructCycle(adjMoments$sMat, adjMoments$tMat)

    # estimate edges from C into D
    Lambda_CD <- R_DC %*% (diag(length(D)) - Lambda_DD$Lambda)

    # for dInd in D, get p-values for edge of parent into dInd
    for(dInd in 1:length(D)){

      pvals_CD[1:length(C), dInd] <- disjointCycles::pruneHelperCycle(C, D, dInd, Y, Lambda_CD, Lambda_DD$Lambda)

    }

    # Fill in non-pruned edgeweights and edges
    Lambda[1:length(C), 1:length(D)] <- Lambda_CD
    adj_CD[1:length(C), 1:length(D)] <- 1

    # fill in edgeweights and edges estimated from cycle
    Lambda[length(C)+ 1:length(D), 1:length(D)] <- Lambda_DD$Lambda
    adj_CD[length(C) + 1:length(D), 1:length(D)] <- Lambda_DD$adj

    # Set p-values for (D,D) elements to -1 since they aren't being pruned
    pvals_CD[length(C)+ 1:length(D), 1:length(D)] <- -1

    return(list(adj_CD = adj_CD, Lambda_CD = Lambda, pvals = pvals_CD))
  }


}


#' Helper function to test existence of edges from nodes in preceeding layer
#' used when D is cycle
#'
#' @param C the set of nodes in the preceeding layer
#' @param D the set of nodes in the current unit (must be cycle)
#' @param dInd the specific index in D to test
#' @return the pvalues (unadjusted for multiple testing) for any edge to be tested. If the value is -1 this indicates
#'       that there was no test for the parent in this stage. This is because it is part of a cycle so the
#'       edge was forced to exist or not exist based on the greedy algorithm or because one does not
#'       preceed the other in the topological layers.
#'
pruneHelperCycle <- function(C, D, dInd, Y, lambda_CD, lambda_DD){

  # parent of D[dInd]
  d2 <- which(lambda_DD[, dInd] != 0)

  pvals <- rep(0, length(C))

  # cat("C: ")
  # cat(C)
  # cat("; D: ")
  # cat(D)
  # cat("; dInd: ")
  # cat(dInd)
  # cat("lambda_CD: ")
  # cat(lambda_CD)
  # cat("\n")

  for(cTest in 1:length(C)){

    # Null hypothesis has Lambda[C, dInd] with 0 for cTest
    lambda_Cd_null <- lambda_CD[, dInd]
    lambda_Cd_null[cTest] <- 0

    # Form errors under null hypothesis
    # scaling errors to prevent singular A
    errsD <- scale(Y[, D[dInd], drop = F] - Y[, C, drop = F] %*% lambda_Cd_null -
      Y[, D, drop = F] %*% lambda_DD[, dInd])


    # construct A matrix
    A <- matrix(0, length(C) + 1, length(C) + 1)

    # C_mod simply reorders C so that cTest appears first
    C_mod <- c(C[cTest], setdiff(C, C[cTest]))

    for(c1 in 1:length(C)){
      for(c2 in 1:length(C)){
        A[c1, c2] <- mean(-Y[ ,C_mod[c2]] * Y[ , C_mod[c1]]^2)
      }
      A[c1, length(C) + 1] <- mean(-Y[ ,D[d2]] * Y[ ,C_mod[c1]]^2)
    }


    for(c2 in 1:length(C)){
      A[length(C) + 1, c2] <- mean(-2 * errsD * Y[, C_mod[1]] * Y[, C_mod[c2]])
    }
    A[length(C) + 1, length(C) + 1] <- mean(-2 * errsD * Y[, C_mod[1]] * Y[, D[d2]])


    m <- cbind(Y[,C_mod, drop = F]^2 * c(errsD), errsD^2 * Y[,C_mod[1], drop = F])

    g <- m[, 1 , drop = F] - t(A[1, -1, drop = F] %*% solve(A[-1, -1, drop = F]) %*% t(m[, -1, drop = F]))
    pvals[cTest] <- emplik::el.test(g, mu = 0)$Pval

  }
  return(pvals)

}

#' Helper function to test existence of edges from nodes in preceeding layer
#' used when D is cycle
#'
#' @param C the set of nodes in the preceeding layer
#' @param D the set of nodes in the current unit (must be cycle)
#' @param dInd the specific index in D to test
#' @return the pvalues (unadjusted for multiple testing) for any edge to be tested. If the value is -1 this indicates
#'       that there was no test for the parent in this stage. This is because it is part of a cycle so the
#'       edge was forced to exist or not exist based on the greedy algorithm or because one does not
#'       preceed the other in the topological layers.
#'
pruneHelperSingle <- function(C, D, Y, lambda_CD, pr = "infFunc"){

  # cat("C: ")
  # cat(C)
  # cat("\n D: ")
  # cat(D)
  # cat("\n")

  pvals <- matrix(0, length(C))






      if(pr == "infFunc"){


        # if C is length 1, then there are no nuissance parameters to consider
        if(length(C) == 1){

          pvals[1] <- emplik::el.test(Y[, C] *Y[, D] , mu = 0)$Pval

        } else {


          for(cTest in 1:length(C)){

            # Null hypothesis has 0 for cTest
            lambda_Cd_null <- lambda_CD
            lambda_Cd_null[cTest] <- 0

            # Form errors under null hypothesis
            errsD <- scale(Y[, D, drop = F] - Y[, C, drop = F] %*% lambda_Cd_null)


            # construct A matrixC <
            A <- matrix(0, length(C), length(C))

            # C_mod simply reorders C so that cTest appears first
            C_mod <- c(C[cTest], setdiff(C, C[cTest]))

            for(c1 in 1:length(C)){
              for(c2 in 1:length(C)){
                A[c1, c2] <- mean(-Y[ ,C_mod[c2]] * Y[ ,C_mod[c1]]^2)
              }

            }

            # form estimating equations m
            m <- Y[, C_mod, drop = F]^2 * c(errsD)

            # cat("C_mod: ")
            # cat(C_mod)
            # cat("; M: ")
            # cat(dim(m))
            # cat("; A: ")
            # cat(dim(A))
            # cat("\n")

            # marginalize to only include first
            g <- m[, 1, drop = F] - t(A[1, -1, drop = F] %*% solve(A[-1, -1, drop = F]) %*% t(m[, -1, drop = F]))
            pvals[cTest] <- emplik::el.test(g, mu = 0)$Pval

          }



        }





      } else if(pr == "chisq"){


        # D is singleton

          z <- colMeans(Y[, C] *Y[, D])
          pvals[1] <- pchisq(z %*% solve(cov(Y[, C] * Y[, D])) %*% t(z), lower.tail = F, df = length(D))


      }



  return(pvals)

}


