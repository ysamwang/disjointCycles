#' Helper function used to select a root cycle from the set C
#'
#' @param C the clique of interest
#' @param D all other nodes
#' @param Y the n x p data matrix
#' @param methodPR the method to be used: chisq, infFunc, or naive
#' @return
#' \itemize{
#' \item B the matrix of edge weights
#' \item Y the n x p data
#' \item errs the error realizations
#' \item scales the variance of the errors
#' \item mu the intercept terms
#' }
#'
pruneHelperRootCycle <- function(C, D, Y, methodPR = "infFunc"){

  if(methodPR == "chisq"){
    res <- Y[, D, drop = F] - Y[, C, drop = F] %*% solve(t(Y[, C, drop = F]) %*% Y[, C, drop = F],  t(Y[, C, drop = F]) %*% Y[, D, drop = F])

    m <- matrix(0, n, length(C) * length(D))
    for(dInd in 1:length(D)){
      inds <- (1:length(C)) +(dInd - 1) * length(C)
      m[, inds] <- res[, dInd] * Y[, C]^2
    }

    return(pchisq(n * colMeans(m) %*% solve(cov(m)) %*% colMeans(m), df = ncol(m), lower.tail = F))


  } else if(methodPR == "infFunc"){

    # Note: Does not include intercept
    res <- Y[, D, drop = F] - Y[, C, drop = F] %*% solve(t(Y[, C, drop = F]) %*% Y[, C, drop = F],  t(Y[, C, drop = F]) %*% Y[, D, drop = F])

    A <- matrix(0, 2 * length(C) * length(D), length(C) * length(D))
    Ablock <- matrix(0, length(C), length(C))

    for(cInd1 in 1:length(C)){
      for(cInd2 in 1:length(C)){
        Ablock[cInd1, cInd2] <- mean(-Y[, C[cInd1]]^2 * Y[, C[cInd2]])
      }
    }



    for(dInd in 1:length(D)){
      inds <- (1:length(C)) +(dInd - 1) * length(C)
      A[inds, inds] <- Ablock

      for(cInd1 in 1:length(C)){
        for(cInd2 in 1:length(C)){

          A[length(C) * length(D) + inds, inds] <- -2 * t(Y[, C, drop = F] * res[, dInd]) %*% Y[, C, drop = F] / n

        }
      }

    }
    A <- cbind(rep(1, length(C) * length(D) * 2), A)


    m <- matrix(0, n, 2 * length(C) * length(D))
    for(dInd in 1:length(D)){
      inds <- (1:length(C)) +(dInd - 1) * length(C)
      m[, inds] <- res[, dInd] * Y[, C]^2
      m[, inds + length(C) * length(D)] <- res[, dInd]^2 * Y[, C]
    }

    activeInds <- 1:(length(C) * length(D))
    g <- m[, activeInds] - t(A[activeInds, -1] %*% solve(A[-c(activeInds), -1]) %*% t(m[, -activeInds]))

    return(emplik::el.test(g, mu = rep(0, ncol(g)))$Pval)



  } else {
    # Note: Does not include intercept
    res <-  Y[, D, drop = F] - Y[, C, drop = F] %*% solve(t(Y[, C, drop = F]) %*% Y[, C, drop = F],  t(Y[, C, drop = F]) %*% Y[, D, drop = F])

    ## Test that moments are 0 using empirical likelihood
    el.mat <- matrix(0, n, length(C) * length(D))

    for(dInd in 1:length(D)){

      el.mat[,  (1:length(C)) + (dInd - 1) * length(C)] <- res[, dInd] * Y[, C]^2
    }

    # Adjusted EL pseudo observation
    # el.mat <- rbind(el.mat, -colMeans(el.mat) * 1/2 * log(n) )
    return(emplik::el.test(el.mat, mu = rep(0, length(C) * length(D)))$Pval)


  }

}

