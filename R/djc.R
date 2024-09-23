djcGetOrdering <- function(Y, alpha2 = .01, alpha3 = .01, alphaR = .01, methodRoot = "indDelta", methodDet2 = "indDelta",
                           rescaleData = F, verbose = F){


  # CLEANUP
  # alpha2 = .01; alpha3 = .01; alphaR = .01; methodRoot = "indDelta"; methodDet2 = "indDelta"; rescaleData = T

  p <- ncol(Y)
  unordered <- 1:p
  topOrder <- list()


  # First pass
  roots <- djc_recur(Y, alpha2 = alpha2, alpha3 = alpha3, alphaR = alphaR,
                     methodRoot = methodRoot, methodDet2 = methodDet2,
                     rescaleData = rescaleData)

  # First layer
  layer <- 1
  topOrder[[layer]] <- roots$roots

  # Get remaining unordered nodes
  unordered <- setdiff(unordered, unlist(roots$roots))

  if(verbose){

    cat("=== Layer: ")
    cat(layer)
    cat(" ===\n")
    cat("Roots: ")
    print(roots$roots)
    cat("\n")
  }


  while(length(unordered) > 0){
    layer <- layer + 1



    roots <- djc_recur(roots$newY, alpha2 = alpha2, alpha3 = alpha3,
                       alphaR = alphaR, methodRoot = methodRoot,
                       methodDet2 = methodDet2)



    # djc_recur returns indices corresponding to 1:ncol(newY)
    # indexBack maps those back to the original indices
    indexBack <- lapply(roots$roots, function(x){unordered[x]})

    # Put those into the topological ordering
    topOrder[[layer]] <- indexBack

    # remove from
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
  topOrder <- orderList(topOrder)

  return(topOrder)


}


### Part 1
djc_recur <- function(Y, alpha2 = .01, alpha3 = .01, alphaR = .01, methodRoot = "indDelta", methodDet2 = "indDelta",
                       rescaleData = F){




  p <- ncol(Y)
  n <- nrow(Y)

  # If only 1 node remains, it must be the root
  if(p == 1){
    return(list(roots = list(c(1)), newY = NULL))
  }

  # If data should be rescaled before estimating roots
  if(rescaleData){
    Y <- scale(Y)
  }

  # Calculate 2nd and 3rd moments
  moments <- .calcSandT(Y)


  #### Get singleton Roots ####

  # get the nodes u for which d_{u,v} == 0 for all v
  # i.e., null is not rejected when is pVal >= alpha2
  # If methodRoot = indDelta, each test is done individually then p-values are aggregated using holm
  # if methodRoot = jointDelta, a single chi^2 test is done
  # singletonRoots should be a vector

  singletonRoots <- which(p.adjust(sapply(1:p, function(u){.test2joint(u, moments$sMat, moments$tMat, Y, method = methodRoot)$pval}), method = "holm") > alpha2)


  # If there is at least 1 singleton root certified
  if(length(singletonRoots) > 0){

    # If there are any remaining unordered nodes, then adjust for certified roots
    if(length(singletonRoots) < p){

      # Regress out roots and return
      newY <- matrix(0, n, p)
      for(v in setdiff(1:p, singletonRoots)){

        # regress away all singleton roots
        newY[, v] <- lm(Y[, v] ~ Y[, singletonRoots] - 1)$res
      }

      # If there are any singleton nodes, then return the roots and adjusted Y
      return(list(roots = as.list(singletonRoots), newY = newY[, -singletonRoots, drop = F]))

    }  else {

      ## All nodes are roots, so no need to adjust Y
      return(list(roots = as.list(singletonRoots), newY = NULL))
    }

  }


  #### If no singleton roots are found, then test for root cycles ####

  # Will hold p-values for testing null hypothesis that
  #   d_{uv}^{2x2} == 0 or d_{vu}^{2x2} == 0
  # only fills in upper triangle
  det2 <- matrix(1, p, p)

  for(u in 1:(p-1)){
    for(v in (u+1):p){

      # P-value for checking null that det_{u,v} == 0 or det_{v,u} == 0
      # if methodDet2 = jointDelta, then tests if det_{u,v} x det_{v,u} == 0
      # if methodDet2 = indDelta, then tests det_{u,v} == 0 and det_{v,u} == 0 individually
      #   and returns maximum p-value
      if(u != v){
        det2[u , v] <- .test2both(u, v, moments$sMat, moments$tMat, Y, method = methodDet2)$p
      }
    }
  }
  # adjust using holm
  det2[upper.tri(det2)] <- p.adjust(det2[upper.tri(det2)], method = "holm")




  ## Get pairs such that d_{u,v} and d_{v,u} are simultaneously non-zero
  # If no pairs are jointly non-zero at alpha2, then raise threshold
  # to the minimum p-value in det2 so at least 1 pair is jointly certified
  threshold2 <- max(alpha2, min(det2))
  nonZero2 <- (det2 <= threshold2)


  # Will hold p-values for testing if d_{uv}^{3x3} == 0
  # only fill in upper triangle
  det3 <- matrix(0, p, p)
  for(u in 1:(p-1)){
    for(v in (u+1):p){
      if(nonZero2[u, v] == 1){
        det3[u , v] <- .test3(u, v, moments$sMat, moments$tMat, Y)$pval
      }
    }
  }
  # adjust using holm
  det3[upper.tri(det3)] <- p.adjust(det3[upper.tri(det3)], method = "holm")


  # Set the threshold to alpha3
  # or decrease to the threshold to the largest p-value if it is smaller
  # than alpha3 so that there is at least 1 pair connected nodes
  threshold3 <- min(alpha3, max(det3))

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

      # Note: Does not include intercept
      res <- Y[, C2] - Y[, cliqueK] %*% solve(t(Y[, cliqueK]) %*% Y[, cliqueK],  t(Y[, cliqueK]) %*% Y[, C2])

      ## Test that moments are 0 using empirical likelihood
      el.mat <- matrix(0, n, length(cliqueK) * length(C2))

      for(s in 1:length(cliqueK)){
        el.mat[, (length(C2) * (s-1) + 1):(s * length(C2)) ] <- Y[, cliqueK[s]]^2 * res
      }

      # Adjusted EL pseudo observation
      el.mat <- rbind(el.mat, -colMeans(el.mat) * 1/2 * log(n) )
      pvalRoots[k] <- emplik::el.test(el.mat, mu = rep(0, length(cliqueK) * length(C2)))$Pval

    }

    # Adjust with holm
    pvalRoots <- p.adjust(pvalRoots, method = "holm")



    ## if all cycles are rejected, then set root cycle to be all combined
    ## i.e., all put into 1 big cycle
    if(sum(pvalRoots >= alphaR) == 0 ){

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



    newY <- matrix(0, n, p)

    for(v in 1:p){
        if(! v %in% unlist(rootCycles)){
          # Assumes centered Y so does not include intercept
          newY[, v] <- lm(Y[, v] ~ Y[, unlist(rootCycles)] - 1)$res
        }
    }

    if(testMode){
      cat("Done adjusting for maxCliques")
      cat("\n")
    }

    # rootCycles should be a list of vectors
    return(list(roots = rootCycles, newY = newY[, -unlist(rootCycles), drop = F]))
  }

}





