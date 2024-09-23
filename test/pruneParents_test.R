runInd <- 1
args <- commandArgs(TRUE)
for(i in 1:length(args)){
  eval(parse(text = args[[i]]))
}
print(runInd)


pruneParents1 <- function(C, D, sMat, tMat, Y, alpha, oracleDat = NULL,
                          oracleLambdaD = NULL, infFunc = F){

  testAll <- F
  adj_CD <- matrix(0, length(C) + length(D), length(D))
  adj_CD[1:length(C), 1:length(D)] <- 1

  Lambda <- matrix(0, length(C) + length(D), length(D))
  pvals_CD <- matrix(0, length(C), length(D))

  R_DC <- sMat[D,C] %*% solve(sMat[C,C])


  if(!is.null(oracleDat)){

    adjMoments <- calcSandT(oracleDat)
    Lambda_DD <- disjointCycles::constructCycle(adjMoments$sMat, adjMoments$tMat)

  } else if(!is.null(oracleLambdaD)){

    Lambda_DD <- list(adj = (oracleLambdaD != 0), Lambda = oracleLambdaD)

  } else {

    Y_D_adjusted <- Y[, D] - Y[, C] %*% t(R_DC)
    adjMoments <- calcSandT(Y_D_adjusted)
    Lambda_DD <- disjointCycles::constructCycle(adjMoments$sMat, adjMoments$tMat)

  }


  Lambda_CD <- t(R_DC) %*% (diag(length(D)) - Lambda_DD$Lambda)


  if(infFunc){
    for(dInd in 1:length(D)){

      pvals_CD[, dInd] <- pruneHelper(C, D, dInd, Y, Lambda_CD, Lambda_DD$Lambda)

    }

  } else {

    for(d in 1:length(D)){

      for(c in 1:length(C)){

        Lambda_CD_null <- Lambda_CD
        Lambda_CD_null[c, d] <- 0

        errD_null <- Y[, D[d]] - Y[, C] %*% Lambda_CD_null[, d] - Y[, D] %*% Lambda_DD$Lambda[, d]


        ### Test independence from all C, or just specific c
        if(testAll){
          pvals_CD[c, d] <- disjointCycles::checkIndEL(errD_null, scale(Y[, C]^2))

        } else {

          pvals_CD[c, d] <- disjointCycles::checkIndEL(errD_null, scale(Y[, C[c]]^2))

        }


      }

    }

  }

  pvals_CD_raw <- pvals_CD
  pvals_CD <- matrix(p.adjust(pvals_CD, method = "holm"), nrow = length(C))



  Lambda[1:length(C), 1:length(D)] <- Lambda_CD * (pvals_CD <= alpha)
  adj_CD[1:length(C), 1:length(D)] <- (pvals_CD <= alpha) + 0

  Lambda[length(C)+ 1:length(D), 1:length(D)] <- Lambda_DD$Lambda
  adj_CD[length(C) + 1:length(D), 1:length(D)] <- Lambda_DD$adj

  return(list(adj_CD = adj_CD, Lambda_CD = Lambda, pvals_raw = pvals_CD_raw))


}


pruneHelper <- function(C, D, dInd, Y, lambda_CD, lambda_DD){

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

    g <- m[, 1] - t(A[1, -1] %*% solve(A[-1, -1]) %*% t(m[, -1]))
    pvals[cTest] <- emplik::el.test(g, mu = 0)$Pval

  }
  return(pvals)

}


###########################
library("disjointCycles")
sample.size <- 1000
rep.runs <- 10
n.list <- c(10, 25, 50) * 1000
c.list <- c(4)
alpha <- .05
param.grid <- expand.grid(rep(n.list, sample.size / rep.runs), c.list)
p <- 10

cycleSize <- param.grid[runInd, 2]
n <- param.grid[runInd, 1]

set.seed(100)
trueAdj <- matrix(0, p, p)
trueAdj[matrix(c((p-cycleSize + 1):p, c(p-cycleSize + 2):p, p-cycleSize + 1) , ncol = 2)] <- 1
trueAdj[1:(p-cycleSize), (p-cycleSize + 1):p] <- rbinom((p-cycleSize) * cycleSize, size = 1, prob = 1/2)
Lambda <- trueAdj * runif(p^2, .3, .9) * sample(c(-1, 1), size = p^2, replace = T)
trueGraph <- list(as.list(1:(p-cycleSize)),list(c(p-cycleSize + 1):p))
C <- 1:(p-cycleSize)
D <- (p-cycleSize + 1):p
rec <- rec1 <- rec2 <- rec3 <-  matrix(0, nrow = rep.runs, ncol = length(C) * length(D))

tester <- Lambda[C,D]



set.seed(1000 + runInd)

for(i in 1:rep.runs){
  data <- disjointCycles::rLSEM(p, n, dist = "gamma", LambdaIn = Lambda, lowScale = .8)
  oracleReg <- solve(data$sigma[C,C]) %*% data$sigma[C,D]
  Y <- data$Y
  oracleDat <- Y[,D] - Y[, C] %*% oracleReg

  moments <- calcSandT(Y)
  out <- pruneParents1(C, D, moments$sMat, moments$tMat, Y, alpha)
  rec[i, ] <- c(out$pvals_raw)


  out_dat <- pruneParents1(C, D, moments$sMat, moments$tMat, Y, alpha,
                           oracleDat = oracleDat)

  rec1[i, ] <- c(out_dat$pvals_raw)

  out_lambda <- pruneParents1(C, D, moments$sMat, moments$tMat, Y,
                              alpha, oracleLambdaD = Lambda[D,D])

  rec2[i, ] <- c(out_lambda$pvals_raw)

  out_Inf <- pruneParents1(C, D, moments$sMat, moments$tMat, Y,
                              alpha, infFunc = T)

  rec3[i, ] <- c(out_Inf$pvals_raw)

  cat(i)
  cat("\n")
}


out <- data.frame(n = n, cycleSize = cycleSize, type = rep(c("naive", "oracleDat", "oracleLambda", "infFunc"), each = rep.runs),
           rbind(rec, rec1, rec2, rec3))

write.csv(out, paste("../results/pvals/pvals_",runInd, ".csv", sep = ""), row.names = F)







