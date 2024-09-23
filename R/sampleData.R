rMixNorm <- function(n, mu, sd, piMix){
  # Flip mu positive/negative

  mu <- mu * sample(c(-1, 1), size = 1)

  totalVar <- sd^2 + mu^2 - mu^2 * ( 2 * piMix - 1)^2
  X <- rnorm(n, mean = mu * (2 * rbinom(n, size = 1, prob = piMix) - 1), sd = sd)


  return((X - (2 * piMix - 1) * mu) / sqrt(totalVar) )
}



cycleChain <- function(p, cycleSize = 2, lowEdge = .5, highEdge = .9, posNeg = T,
                       parentProb = 1, uniqueTop = T){

  Lambda <- matrix(0, p, p)

  # If layer topological ordering is unique, then add an edge from u -> u+1


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

  if(uniqueTop){

    Lambda[matrix(c(1:(p-1), 2:p), ncol = 2)] <- runif(p-1, lowEdge, highEdge) *
      sample(c(-(2 * posNeg - 1), 1), p-1, replace = T)

  }

 return(Lambda)
}


ex3_2 <- function(lowEdge = .5, highEdge = .9, posNeg = T){

  p <- 8
  Lambda <- adjMat <- matrix(0, p, p)

  adjMat[matrix(c(1,3,
                  2,4,
                  3,4,
                  5,3,
                  4,5,
                  5,6,
                  6,7,
                  7,8,
                  8,6), byrow = T, ncol = 2)] <- 1


  Lambda <- adjMat * runif(p^2, lowEdge, highEdge) * sample(c(-(2 * posNeg - 1), 1), p^2, replace = T)


  return(Lambda)
}






