#' Sample data from an additive noise model
#'
#' @param p number of variables
#' @param n number of observations
#' @param lowScale lower bound on variance of error terms
#' @param highScale upper bound on variance of error terms
#' @param lowEdge lower bound on edgeWeights
#' @param highEdge upper bound on edgeWeights
#' @param dist the distribution of the error terms. Choices are: "gauss", "unif", "lognorm", "gamma", "weibull", "laplace", "mixed"
#' @param AdjMat an adjacency matrix which can be passed in to specify the edges instead of randomly selecting a graph.
#'  If A[i,j] == 1, this indicates i --> j
#' @param LambdaIn a matrix of edges which can be passed in instead of randomly drawing edges
#' @param scalesIn a vector of error variances which can be passed in instead of randomly drawing
#' @param posAndNeg whether edges should be both positive and negative or only positive
#' @return
#' \itemize{
#' \item B the matrix of edge weights
#' \item Y the n x p data
#' \item errs the error realizations
#' \item scales the variance of the errors
#' \item mu the intercept terms
#' }
#'
rLSEM <- function(p, n,
                 lowScale = 1, highScale = 1,
                 lowEdge = .3, highEdge = 1,
                 dist = "gauss", AdjMat = NULL,
                 LambdaIn = NULL, scalesIn = NULL, posAndNeg = T) {

  if(posAndNeg){
    signsAllowed <- c(-1, 1)
  } else {
    signsAllowed <- c(1)
  }

  if(is.null(LambdaIn)){

    Lambda <- AdjMat * matrix(sample(signsAllowed, size = p^2, replace = T) * runif(p^2, min = lowEdge, max = highEdge), p, p)
    } else {

      Lambda <- LambdaIn

  }

  if(is.null(scalesIn)){

    scales <- runif(p, lowScale, highScale)

  } else {

    scales <- scalesIn
  }




  if(dist == "gauss") {

    errs <- matrix(rnorm(n * p) * scales, nrow = n, ncol = p, byrow = T)

  } else if (dist == "lognorm") {

    errs <- matrix((exp(rnorm(n * p)) - exp(1/2)) / sqrt(exp(1) * (exp(1) - 1)) * scales, nrow = n, ncol = p , byrow = T)

  } else if (dist == "gamma"){

    errs <- matrix((rgamma(n * p, 1, 1) - 1) * scales, nrow = n, ncol = p, byrow = T)

  } else if (dist == "weibull"){
    a <- 3/4

    errs <- matrix( ((rweibull(n * p, shape = a, 1) - gamma(1 + 1/a)) / sqrt(gamma(1 + 2/a) - gamma(1 + 1/a)^2)) * scales,
                    nrow = n, ncol = p, byrow = T)

  } else if (dist == "mixedNorm"){

    errs <- matrix(0, nrow = n, ncol = p)

    piMix <- .1

    for(v in 1:p){
      errs[, v] <- rMixNorm(n, 2, .1, piMix) * scales[v]
    }

  } else if (dist == "beta"){
    a <- .5
    b <- 6
    errs <- matrix((rbeta(n * p, a, b) - a / (a + b)) * scales, nrow = n, ncol = p, byrow = T)

  }



  Y <- solve(diag(rep(1, p)) - t(Lambda), t(errs))
  Y <- t(Y)

  sigma <- t(solve(diag(p)- Lambda)) %*% diag(scales) %*% solve(diag(p) - Lambda)


  return(list(Lambda = Lambda, Y = Y, errs = errs, scales = scales, sigma = sigma))
}
