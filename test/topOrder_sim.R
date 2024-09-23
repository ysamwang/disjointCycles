runInd <- 1
args <- commandArgs(TRUE)
for(i in 1:length(args)){
  eval(parse(text = args[[i]]))
}
print(runInd)

library("disjointCycles")
sample.size <- 500
rep.runs <- 25
n.list <- c(10, 25, 50) * 1000
c.list <- c(2, 3, 4)
a.list <- c(.05, .01, .001, 1e-4, 1e-5)
param.grid <- expand.grid(rep(n.list, sample.size / rep.runs), c.list, a.list)
# 450

p <- 12
n <- param.grid[runInd, 1]
cycleSize <- param.grid[runInd, 2]
alpha <- param.grid[runInd, 3]

rec <- matrix(0, rep.runs, 10)
for(i in 1:rep.runs){
  Lambda <- cycleChain(p, cycleSize = cycleSize, lowEdge = .5, highEdge = .9, parentProb = 1/2)

  trueAdj <- (Lambda != 0) +0
  trueGraph <- lapply(unname(split(1:p, rep(1:ceiling(p/cycleSize),
                                                          each = cycleSize)[1:p])),
                                    list)



  data <- rLSEM(p, n, dist = "gamma", LambdaIn = Lambda, lowScale = .8)
  est_ord <- djcGetOrder(scale(data$Y), verbose = F, alpha2 = alpha, alpha3 = alpha, alphaR = alpha)
  rec[i, 1] <- compareOrders(trueGraph, est_ord)

  ## Estimate edges given true layers
  est_Edges <- djcGetEdges(trueGraph, scale(data$Y, scale = T), alpha = alpha)
  rec[i, 2] <- sum(est_Edges$adjMat *trueAdj) / sum(trueAdj)
  rec[i, 3] <- sum((est_Edges$adjMat + trueAdj) == 0) / sum(trueAdj == 0)



  ## Estimate edges given estimated layers
  final_ord <- djcGetEdges(est_ord, scale(data$Y, scale = T), alpha = alpha)
  rec[i, 4] <- sum(final_ord$adjMat *trueAdj) / sum(trueAdj)
  rec[i, 5] <- sum((final_ord$adjMat + trueAdj) == 0) / sum(trueAdj == 0)


  data <- rLSEM(p, n, dist = "mixedNorm", LambdaIn = Lambda, lowScale = .8)
  est_ord <- djcGetOrder(scale(data$Y), verbose = F, alpha2 = alpha, alpha3 = alpha, alphaR = alpha)
  rec[i, 6] <- compareOrders(trueGraph, est_ord)

  ## Estimate edges given true layers
  est_Edges <- djcGetEdges(trueGraph, scale(data$Y, scale = T), alpha = alpha)
  rec[i, 7] <- sum(est_Edges$adjMat *trueAdj) / sum(trueAdj)
  rec[i, 8] <- sum((est_Edges$adjMat + trueAdj == 0)) / sum(trueAdj == 0)

  ## Estimate edges given estimated layers
  final_ord <- djcGetEdges(est_ord, scale(data$Y, scale = T), alpha = alpha)
  rec[i, 9] <- sum(final_ord$adjMat *trueAdj) / sum(trueAdj)
  rec[i, 10] <- sum((final_ord$adjMat + trueAdj) == 0) / sum(trueAdj == 0)


}

colnames(rec) <- paste(rep(c("gamma", "mixNorm"), each = 5),
                       rep(c("order", "oracle_sen", "oracle_spc", "est_sen", "est_spc"), times= 2),
                       sep = "_")
outTab <- data.frame(cycleSize, alpha, n, rec)

write.csv(outTab, paste("../results/cycleChain/cycleChain_",runInd, ".csv", sep = ""), row.names = F)

