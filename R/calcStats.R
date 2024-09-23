calcSandT <- function(Y){
  n <- nrow(Y)
  p <- ncol (Y)

  sMat <- t(Y) %*% Y / n

  tMat <- array(0, dim = c(p, p, p))

  for(i in 1:p ){
    for(j in 1:p ){
      for(k in 1:p ){
        tMat[i, j, k] <- mean(Y[, i] * Y[, j] * Y[, k])
      }
    }
  }

  return(list(sMat = sMat, tMat = tMat))
}




.test2joint <- function(u, sMat, tMat, Y, method = "indDelta"){


  .test2ind <- function(u, v){

    ## asymptotic Normality through delta method
    det2 <- sMat[u, u] * tMat[u, u, v] - sMat[u, v] * tMat[u, u, u]

    momentCov <- cov(cbind(Y[,u] * Y[, u], Y[,u] * Y[, v], Y[,u]^3, Y[,u]^2 * Y[, v]))
    jacobian <- c(tMat[u,u,v], -tMat[u,u,u], -sMat[u,v], sMat[u,u])


    deltaVar <- t(jacobian) %*% momentCov %*% jacobian

    chiDelta <- n *det2^2 / deltaVar
    pValDelta <- pchisq(chiDelta, df = 1, lower.tail = F)


    return(pValDelta)

  }





  n <- nrow(Y)
  p <- ncol(Y)


  # Chi-squared test for vector (d_{u,v}^{2} : v != u) == 0
  if(method == "jointDelta"){

    Sigma <- matrix(0, 2 * p, 2 * p)
    for(v in 1:p){
      for(w in 1:p){

        Sigma[v, w] <- mean(Y[, u]^2 * Y[, v] * Y[, w]) - sMat[u, v] * sMat[u,w]
        Sigma[v + p, w] <- mean(Y[, u]^3 * Y[, v] * Y[, w]) - tMat[u, u, v] * sMat[u, w]
        Sigma[v, w + p] <- mean(Y[, u]^3 * Y[, v] * Y[, w]) - sMat[u, v] * tMat[u, u, w]
        Sigma[v + p, w + p] <- mean(Y[, u]^4 * Y[, v] * Y[, w]) - tMat[u, u, v] * tMat[u, u, w]

      }
    }


    J <- matrix(0, p, 2*p)

    for(v in 1:p){
      J[v, u] <- tMat[u, u, v]
      J[v, v] <- -tMat[u, u, u]
      J[v, p+u] <- -sMat[u, v]
      J[v, p+v] <- sMat[u, u]
    }

    J <- J[-u, , drop = F]


    finalCov <- J %*% Sigma %*% t(J)

    det2 <- rep(0, p)
    for(v in 1:p){
      if (v != u){
        det2[v] <- sMat[u, u] * tMat[u, u, v] - sMat[u, v] * tMat[u, u, u]
      }
    }

    det2 <- det2[-u]
    chiDelta <- n * t(det2) %*% solve(finalCov) %*% det2
    pValDelta <- pchisq(chiDelta, df = p-1, lower.tail = F)

    return(list(det2 = det2, chi2 = chiDelta, pval = pValDelta))

  } else if(method == "indDelta"){
    # Test each d_{u,v}^{2} == 0 and combine using holm

    pvals <- p.adjust(sapply(setdiff(1:p, u), function(v){.test2ind(u, v)}), method = "holm")

    # cat("min: ")
    # cat(u)
    # cat("; p = ")
    # cat(p)
    # cat("\n")
    # cat(pvals)
    # cat("\n")
    return(list(pval = min(pvals)))


  }


}




.test2both <- function(u, v, sMat, tMat, Y, method = "indDelta"){



  .test2ind <- function(u, v){

    ## asymptotic Normality through delta method
    det2 <- sMat[u, u] * tMat[u, u, v] - sMat[u, v] * tMat[u, u, u]

    momentCov <- cov(cbind(Y[,u] * Y[, u], Y[,u] * Y[, v], Y[,u]^3, Y[,u]^2 * Y[, v]))
    jacobian <- c(tMat[u,u,v], -tMat[u,u,u], -sMat[u,v], sMat[u,u])


    deltaVar <- t(jacobian) %*% momentCov %*% jacobian

    chiDelta <- n *det2^2 / deltaVar
    pValDelta <- pchisq(chiDelta, df = 1, lower.tail = F)


    return(pValDelta)

  }

  if (method == "jointDelta"){

  ## asymptotic Normality of product through delta method
  det2uv <- sMat[u, u] * tMat[u, u, v] - sMat[u, v] * tMat[u, u, u]
  det2vu <- sMat[v, v] * tMat[u, v, v] - sMat[v, u] * tMat[v, v, v]


  momentCov <- cov(cbind(Y[, u]^2,
                         Y[, u] * Y[,v],
                         Y[, v]^2,
                         Y[, u]^3,
                         Y[, u]^2*Y[,v],
                         Y[, u]*Y[,v]^2,
                         Y[, v]^3))

  jacobian <- c(tMat[u,u,v] * det2vu,
                -sMat[u, u] * tMat[u, u, v] * tMat[v, v, v] - sMat[v, v] * tMat[u, u, u] * tMat[u, v, v] +
                  2 * sMat[u, v] * tMat[u,u,u] * tMat[v, v, v],
                tMat[u, v, v] * det2uv,
                -sMat[u,v] * det2vu,
                sMat[u, u] * det2vu,
                sMat[v, v] * det2uv,
                -sMat[u,v] * det2uv)


  deltaVar <- t(jacobian) %*% momentCov %*% jacobian

  chiDelta <- n * (det2uv * det2vu)^2 / deltaVar
  pValDelta <- pchisq(chiDelta, df = 1, lower.tail = F)


    return(list(pval = pValDelta))

  } else  if (method == "indDelta"){

    return(list(pval = max(c(.test2ind(u,v), .test2ind(v,u)))))

  }

}


######################################3




.test3 <- function(u, v, sMat, tMat, Y) {

  det3 <- sMat[u, u] * (tMat[u, u, v] * tMat[v, v, v] - tMat[u, v, v] * tMat[u, v, v]) -
    sMat[u, v] * (tMat[u, u, u] * tMat[v, v, v] - tMat[u, v, v] * tMat[u, u , v]) +
    sMat[v, v] * (tMat[u, u, u] * tMat[u, v, v] - tMat[u, u, v] * tMat[u, u, v])



  momentCov <- cov(cbind(Y[, u]^2,
                         Y[, u] * Y[,v],
                         Y[, v]^2,
                         Y[, u]^3,
                         Y[, u]^2*Y[,v],
                         Y[, u]*Y[,v]^2,
                         Y[, v]^3))




  jacobian <- c(tMat[u,u,v] * tMat[v,v,v] - tMat[u,v,v] * tMat[u,v,v],
                -(tMat[u,u,u] * tMat[v,v,v] - tMat[u,v,v] * tMat[u,u,v]),
                tMat[u,u,u] * tMat[u,v,v] - tMat[u,u,v] * tMat[u,u,v],
                -sMat[u,v] * tMat[v,v,v] + sMat[v,v] * tMat[u,v,v],
                sMat[u,u] * tMat[v,v,v] + sMat[u,v] * tMat[u,v,v] - 2 * sMat[v,v] * tMat[u,u,v],
                -2 * sMat[u,u] * tMat[u,v,v] + sMat[u,v] * tMat[u,u,v] + sMat[v,v] * tMat[u,u,u],
                sMat[u,u] * tMat[u,u,v] - sMat[u,v] * tMat[u,u,u])


  deltaVar <- t(jacobian) %*% momentCov %*% jacobian

  chiDelta <- n *det3^2 / deltaVar
  pValDelta <- pchisq(chiDelta, df = 1, lower.tail = F)


  return(list(det3 = det3, chi2 = chiDelta, pval = pValDelta))


}


# .d3 <- function(theta) {
#   sMat <- matrix(c(theta[1], theta[2],
#                    theta[2], theta[3]), byrow = T, 2, 2)
#
#   tMat <- array(0, dim = c( 2, 2,2))
#   tMat[1, 1, 1] <- theta[4]
#   tMat[1, 1, 2] <- theta[5]
#   tMat[1, 2, 2] <- theta[6]
#   tMat[2, 2, 2] <- theta[7]
#
#
#
#   det3 <- sMat[1, 1] * (tMat[1, 1, 2] * tMat[2, 2, 2] - tMat[1, 2, 2] * tMat[1, 2, 2]) -
#     sMat[1, 2] * (tMat[1, 1, 1] * tMat[2, 2, 2] - tMat[1, 2, 2] * tMat[1, 1 , 2]) +
#     sMat[2, 2] * (tMat[1, 1, 1] * tMat[1, 2, 2] - tMat[1, 1, 2] * tMat[1, 1, 2])
#
#   return(det3)
#
# }
#
# .jacob <- function(theta) {
#   sMat <- matrix(c(theta[1], theta[2],
#                    theta[2], theta[3]), byrow = T, 2, 2)
#
#   tMat <- array(0, dim = c( 2, 2,2))
#   tMat[1, 1, 1] <- theta[4]
#   tMat[1, 1, 2] <- theta[5]
#   tMat[1, 2, 2] <- theta[6]
#   tMat[2, 2, 2] <- theta[7]
#
#
#
#
#   jacobian <- c(tMat[u,u,v] * tMat[v,v,v] - tMat[u,v,v] * tMat[u,v,v],
#                 -(tMat[u,u,u] * tMat[v,v,v] - tMat[u,v,v] * tMat[u,u,v]),
#                 tMat[u,u,u] * tMat[u,v,v] - tMat[u,u,v] * tMat[u,u,v],
#                 -sMat[u,v] * tMat[v,v,v] + sMat[v,v] * tMat[u,v,v],
#                 sMat[u,u] * tMat[v,v,v] + sMat[u,v] * tMat[u,v,v] - 2 * sMat[v,v] * tMat[u,u,v],
#                 -2 * sMat[u,u] * tMat[u,v,v] + sMat[u,v] * tMat[u,u,v] + sMat[v,v] * tMat[u,u,u],
#                 sMat[u,u] * tMat[u,u,v] - sMat[u,v] * tMat[u,u,u])
#
#   return(jacobian)
#
#
# }





