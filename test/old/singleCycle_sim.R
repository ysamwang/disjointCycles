runInd <- 1
args <- commandArgs(TRUE)
for(i in 1:length(args)){
  eval(parse(text = args[[i]]))
}
print(runInd)

library("disjointCycles")
sample.size <- 200
rep.runs <- 10
n.list <- c(10, 25, 50) * 1000
c.list <- c(2, 3, 4)
a.list <- c(.05, .01, .001, 1e-4, 1e-5)
param.grid <- expand.grid(rep(n.list, sample.size / rep.runs), c.list, a.list)
# 450

p <- 10
n <- param.grid[runInd, 1]
cycleSize <- param.grid[runInd, 2]
alpha <- param.grid[runInd, 3]

rec <- matrix(0, rep.runs, 10)
for(i in 1:rep.runs){

  trueAdj <- matrix(0, p, p)

  trueAdj[matrix(c((p-cycleSize + 1):p, c(p-cycleSize + 2):p, p-cycleSize + 1) , ncol = 2)] <- 1
  trueAdj[1:(p-cycleSize), (p-cycleSize + 1):p] <- rbinom((p-cycleSize) * cycleSize, size = 1, prob = 1/2)

  Lambda <- trueAdj * runif(p^2, .3, .9) * sample(c(-1, 1), size = p^2, replace = T)

  trueGraph <- list(as.list(1:(p-cycleSize)),list(c(p-cycleSize + 1):p))



  data <- rLSEM(p, n, dist = "gamma", LambdaIn = Lambda, lowScale = .8)
  est_ord <- djcGetOrder(scale(data$Y), verbose = F, alpha2 = alpha, alpha3 = alpha, alphaR = alpha)
  rec[i, 1] <- compareOrders(trueGraph, est_ord)

  ## Estimate edges given true layers
  est_Edges <- djcGetEdges(trueGraph, scale(data$Y, scale = T), alpha = alpha)
  rec[i, 2] <- sum(est_Edges$adjMat[1:(p-cycleSize), (p-cycleSize+1):p] *trueAdj[1:(p-cycleSize), (p-cycleSize+1):p]) / sum(trueAdj[1:(p-cycleSize), (p-cycleSize+1):p])
  rec[i, 3] <- sum((est_Edges$adjMat[1:(p-cycleSize), (p-cycleSize+1):p] + trueAdj[1:(p-cycleSize), (p-cycleSize+1):p]) == 0) / sum(trueAdj[1:(p-cycleSize), (p-cycleSize+1):p] == 0)


  # moments <- calcSandT(data$Y)
  # out <- pruneParents1(C = 1:(p-cycleSize), D = (p- cycleSize + 1):p,
  #                           moments$sMat, moments$tMat, data$Y, alpha = .001)
  ## Estimate edges given estimated layers
  final_ord <- djcGetEdges(est_ord, scale(data$Y, scale = T), alpha = alpha)
  rec[i, 4] <- sum(final_ord$adjMat *trueAdj) / sum(trueAdj)
  rec[i, 5] <- sum((final_ord$adjMat + trueAdj) == 0) / sum(trueAdj == 0)


  data <- rLSEM(p, n, dist = "mixedNorm", LambdaIn = Lambda, lowScale = .8)
  est_ord <- djcGetOrder(scale(data$Y), verbose = F, alpha2 = alpha, alpha3 = alpha, alphaR = alpha)
  rec[i, 6] <- compareOrders(trueGraph, est_ord)

  ## Estimate edges given true layers
  est_Edges <- djcGetEdges(trueGraph, scale(data$Y, scale = T), alpha = alpha)
  est_Edges <- djcGetEdges(trueGraph, scale(data$Y, scale = T), alpha = alpha)
  rec[i, 7] <- sum(est_Edges$adjMat[1:(p-cycleSize), (p-cycleSize+1):p] *trueAdj[1:(p-cycleSize), (p-cycleSize+1):p]) / sum(trueAdj[1:(p-cycleSize), (p-cycleSize+1):p])
  rec[i, 8] <- sum((est_Edges$adjMat[1:(p-cycleSize), (p-cycleSize+1):p] + trueAdj[1:(p-cycleSize), (p-cycleSize+1):p]) == 0) / sum(trueAdj[1:(p-cycleSize), (p-cycleSize+1):p] == 0)

  ## Estimate edges given estimated layers
  final_ord <- djcGetEdges(est_ord, scale(data$Y, scale = T), alpha = alpha)
  rec[i, 9] <- sum(final_ord$adjMat *trueAdj) / sum(trueAdj)
  rec[i, 10] <- sum((final_ord$adjMat + trueAdj) == 0) / sum(trueAdj == 0)


}

colnames(rec) <- paste(rep(c("gamma", "mixNorm"), each = 5),
                       rep(c("order", "oracle_sen", "oracle_spc", "est_sen", "est_spc"), times= 2),
                       sep = "_")
outTab <- data.frame(cycleSize, alpha, n, rec)

write.csv(outTab, paste("../results/singleCycle/singleCycle_",runInd, ".csv", sep = ""), row.names = F)

