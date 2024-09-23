p <- 2
n <- 50000
Lambda <- t(disjointCycles::cycleChain(p, cycleSize = p, lowEdge = .5, highEdge = .9))

data <- disjointCycles::rLSEM(p, n, dist = "mixedNorm", LambdaIn = Lambda, lowScale = .8)

out <- cycle2(data$Y)

out$Lambda
out$Lambda1
out$LambdaMa
out$Lambda1Ma
Lambda

svd(out$Lambda)$d
svd(out$Lambda1)$d
svd(out$LambdaMa)$d
svd(out$Lambda1Ma)$d
svd(Lambda)$d

L1 <- t(diag(2) - out$Lambda)
L2 <- t(diag(2) - out$Lambda1)
L3 <- t(diag(2) - rbind(out$Lambda[1, ], out$Lambda1[2, ]))
L4 <- t(diag(2) - rbind(out$Lambda[2, ], out$Lambda1[1, ]))
L5 <- t(diag(2) - Lambda)
L6 <- diag()

e1 <- t(L1 %*% t(data$Y))
checkIndEL(e1[,1, drop = F], scale(e1[, 2, drop = F]^2))
e1 <- t(L2 %*% t(data$Y))
checkIndEL(e1[,1, drop = F], scale(e1[, 2, drop = F]^2))
e1 <- t(L3 %*% t(data$Y))
checkIndEL(e1[,1, drop = F], scale(e1[, 2, drop = F]^2))
e1 <- t(L4 %*% t(data$Y))
checkIndEL(e1[,1, drop = F], scale(e1[, 2, drop = F]^2))
e1 <- t(L5 %*% t(data$Y))
checkIndEL(e1[,1, drop = F], scale(e1[, 2, drop = F]^2))

L_oracle <- matrix(0, 2, 2)
L_oracle[1, 2] <- c(out$Lambda[1,2], out$Lambda[1,2])[which.min()]




est <- disjointCycles::constructCycle(calcSandT(data$Y)$s, calcSandT(data$Y)$t)

est$adjMat
round(abs(Lambda + est$Lambda), 3)


Lambda <- ex3_2(lowEdge = .5, highEdge = .9)
adjMat <- (Lambda != 0) + 0
n <- 50000
data <- rLSEM(p = 8, n, dist = "gamma", LambdaIn = Lambda, lowScale = .8)
Y <- scale(data$Y, scale = F)
alpha <- .1
topOrdering <- djcGetOrder(Y, verbose = F, alpha2 = alpha, alpha3 = alpha, alphaR = alpha)
topOrdering
finalEst <- djcGetEdges(topOrdering, data$Y, alpha = .01)
finalEst$adjMat - adjMat




cycle2 <- function(Y){

  moments <- calcSandT(Y); sMat <- moments$sMat; tMat <- moments$tMat

  Lambda <- Lambda1 <- LambdaMatched <- Lambda1Matched <-  matrix(0, 2, 2)

  u <- 1
  v <- 2


  # -l_(0,1)*l_(1,0)*t_(0,0,1)+l_(0,1)*t_(0,0,0)+l_(1,0)*t_(0,1,1)-t_(0,0,1) = 0,

  a <- sMat[u, u] * tMat[u, u, v] - sMat[u,v] * tMat[u,u,u]
  b <- sMat[v, v] * tMat[u,u,u] - sMat[u,u] * tMat[u, v, v]
  c <- sMat[u, v] * tMat[u, v, v] - sMat[v,v] * tMat[u, u, v]
  # quadratic formula
  Lambda[u, v] <- LambdaMatched[u, v] <-  (-b + sqrt(b^2 - 4 * a * c)) / (2  * a)
  LambdaMatched[v, u] <- (tMat[u, u, v] - LambdaMatched[u, v] * tMat[u, u, u]) / (tMat[u, v, v] -  LambdaMatched[u, v] * tMat[u, u, v])

  Lambda1[u, v] <- Lambda1Matched[u, v] <- (-b - sqrt(b^2 - 4 * a * c)) / (2  * a)
  Lambda1Matched[v, u] <- (tMat[u, u, v] - Lambda1Matched[u, v] * tMat[u, u, u]) / (tMat[u, v, v] -  Lambda1Matched[u, v] * tMat[u, u, v])


  u <- 2
  v <- 1

  a <- sMat[u, u] * tMat[u, u, v] - sMat[u,v] * tMat[u,u,u]
  b <- sMat[v, v] * tMat[u,u,u] - sMat[u,u] * tMat[u, v, v]
  c <- sMat[u, v] * tMat[u, v, v] - sMat[v,v] * tMat[u, u, v]


  # quadratic formula
  Lambda[u, v] <- (-b + sqrt(b^2 - 4 * a * c)) / (2  * a)
  Lambda1[u, v] <- (-b - sqrt(b^2 - 4 * a * c)) / (2  * a)
  return(list(adjMat = matrix(c(0, 1, 1, 0), 2, 2), Lambda = Lambda, Lambda1 = Lambda1, LambdaMatched = LambdaMatched, Lambda1Matched = Lambda1Matched))

}




