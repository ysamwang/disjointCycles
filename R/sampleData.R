#' Generate data from a mixture of normals
#'
#' @param n the sample size
#' @param mu each component is centered at +/- mu
#' @param sd the standard deviation of each component
#' @param piMix the proportion of the sample for one of the components. This is randomly assigned to
#'  either the positive or negative component so the skew may be right or left.
#' @return a vector of length n
#'
#'
rMixNorm <- function(n, mu, sd, piMix){
  mu <- mu * sample(c(-1, 1), size = 1)

  totalVar <- sd^2 + mu^2 - mu^2 * ( 2 * piMix - 1)^2
  X <- rnorm(n, mean = mu * (2 * rbinom(n, size = 1, prob = piMix) - 1), sd = sd)

  return((X - (2 * piMix - 1) * mu) / sqrt(totalVar) )
}


#' Generate Lambda which is a chain of cycles
#'
#' @param p the number of variables
#' @param cycleSize the size of each cycle
#' @param lowEdge the lower bound on the absolute value of each edgeweight
#' @param highEdge the upper bound on the absolute value of each edgeweight
#' @param posNeg whether the edgeweights can be positive and negative
#' @param parentProb the probability of an edge between an ancestor and descendant not in the same cycle
#' @param uniqueTop Whether to enforce a unique topological ordering of cycles by adding an edge from one cycle
#' to the next with probability 1
#' @return
#' a p x p matrix with edgeweight
#'
cycleChain <- function(p, cycleSize = 2, lowEdge = .5, highEdge = .9, posNeg = T,
                       parentProb = 1, uniqueTop = T){

  Lambda <- matrix(0, p, p)

  totalCycles <- floor(p/ cycleSize)

  for(k in 1:totalCycles){

    Lambda[matrix(c(1:cycleSize, 2:cycleSize, 1), ncol = 2) + (k-1) * cycleSize] <-
      runif(cycleSize, lowEdge, highEdge) * sample(c(-(2 * posNeg - 1), 1), cycleSize, replace = T)


    currentCycle <- ((k-1) * cycleSize + 1):(k * cycleSize)


    # If not the last cycle, put edges into
    if((k * cycleSize) < p){
      potentialDesc <- (k * cycleSize + 1):p

      Lambda[currentCycle, potentialDesc] <- runif(length(currentCycle) * length(potentialDesc), lowEdge, highEdge) *
        sample(c(-(2 * posNeg - 1), 1), length(currentCycle) * length(potentialDesc), replace = T) *
        rbinom(length(currentCycle) * length(potentialDesc), size = 1, prob = parentProb)
    }


  }

  # If layer topological ordering is unique, then add an edge from u -> u+1
  if(uniqueTop){

    Lambda[matrix(c(1:(p-1), 2:p), ncol = 2)] <- runif(p-1, lowEdge, highEdge) *
      sample(c(-(2 * posNeg - 1), 1), p-1, replace = T)

  }

 return(Lambda)
}



cycleDag <- function(cycleSize, lowEdge = .5, highEdge = .9, posNeg = T, parentProb = 1/2){

  p <- cycleSize * 3
  Lambda <- adjMat <- matrix(0, p, p)

  for(j in 1:3){
    adjMat[matrix(c(1:cycleSize, 2:cycleSize, 1), ncol = 2) + (j-1) *cycleSize] <- 1
  }

  adjMat[1:(2*cycleSize), -c(1:(2*cycleSize))] <- rbinom(n = 2*cycleSize^2, size = 1, prob = parentProb)


  Lambda <- adjMat * runif(p^2, lowEdge, highEdge) * sample(c(-(2 * posNeg - 1), 1), p^2, replace = T)


  return(Lambda)
}






