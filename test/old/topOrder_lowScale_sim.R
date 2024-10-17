runInd <- 1
args <- commandArgs(TRUE)
for(i in 1:length(args)){
  eval(parse(text = args[[i]]))
}
print(runInd)

library("disjointCycles")
sample.size <- 400
rep.runs <- 80
n.list <- c(10, 25, 50) * 1000
a.list <- c(2e-1, 1e-1, 1e-2, 1e-3, 1e-4)
mt.list <- c("BH", "holm")
pr.list <- c("naive")
c.list <- c(2, 3, 4)
res.list <- c(T, F)
param.grid <- expand.grid(rep(n.list, sample.size / rep.runs), a.list, mt.list, pr.list, c.list, res.list)
dim(param.grid)
# 900

p <- 12
n <- param.grid[runInd, 1]
alpha <- param.grid[runInd, 2]
mt <- as.character(param.grid[runInd, 3])
pr <- as.character(param.grid[runInd, 4])
cycleSize <- param.grid[runInd, 5]
rescaleData <- param.grid[runInd, 6]

rec <- matrix(0, rep.runs, 14)
for(i in 1:rep.runs){
  Lambda <- disjointCycles::cycleChain(p, cycleSize = cycleSize, lowEdge = .5, highEdge = .9,
                                       parentProb = 1/2)

  trueAdj <- (Lambda != 0) +0
  trueGraph <- lapply(unname(split(1:p, rep(1:ceiling(p/cycleSize),
                                                          each = cycleSize)[1:p])),
                                    list)
  cat("Iter ")
  cat(i)
  cat(" : ")
  data <- rLSEM(p, n, dist = "gamma", LambdaIn = Lambda, lowScale = .6, highScale = .8)

  est_ord <- djcGetOrderNew(data$Y, verbose = F, alpha2 = alpha,
                         alpha3 = alpha, alphaR = alpha, pvalAdjMethod = mt,
                         methodPR = pr,
                         rescaleData = rescaleData)
  rec[i, 1] <- compareOrders(trueGraph, est_ord)
  cat(" 1 ")

  # ## Estimate edges given true layers
  est_Edges <- djcGetEdges(trueGraph, data$Y, alpha = alpha,pvalAdjMethod = mt)
  rec[i, 2] <- sum(est_Edges$adjMat *trueAdj) / sum(trueAdj)
  rec[i, 3] <- sum((est_Edges$adjMat + trueAdj) == 0) / sum(trueAdj == 0)
  rec[i, 4] <- sum(est_Edges$adjMat == trueAdj)
  cat(" 4 ")

  ## Estimate edges given estimated layers
  final_ord <- djcGetEdges(est_ord, data$Y, alpha = alpha,pvalAdjMethod = mt)
  rec[i, 5] <- sum(final_ord$adjMat *trueAdj) / sum(trueAdj)
  rec[i, 6] <- sum((final_ord$adjMat + trueAdj) == 0) / sum(trueAdj == 0)
  rec[i, 7] <- sum(final_ord$adjMat == trueAdj)
  cat(" 7 ")


  data <- rLSEM(p, n, dist = "mixedNorm", LambdaIn = Lambda, lowScale = .6, highScale = .8)
  est_ord <- djcGetOrderNew(scale(data$Y), verbose = F, alpha2 = alpha,
                         alpha3 = alpha, alphaR = alpha, pvalAdjMethod = mt, methodPR = pr,
                         rescaleData = rescaleData)
  rec[i, 8] <- compareOrders(trueGraph, est_ord)
  cat(" 8 ")

  # ## Estimate edges given true layers
  est_Edges <- djcGetEdges(trueGraph, data$Y,
                           alpha = alpha, pvalAdjMethod = mt)

  rec[i, 9] <- sum(est_Edges$adjMat *trueAdj) / sum(trueAdj)
  rec[i, 10] <- sum((est_Edges$adjMat + trueAdj) == 0) / sum(trueAdj == 0)
  rec[i, 11] <- sum(est_Edges$adjMat == trueAdj)
  cat(" 11 ")
  #
  #
  # ## Estimate edges given true layers
  final_ord <- djcGetEdges(est_ord, data$Y,
                           alpha = alpha, pvalAdjMethod = mt)

  rec[i, 12] <- sum(final_ord$adjMat *trueAdj) / sum(trueAdj)
  rec[i, 13] <- sum((final_ord$adjMat + trueAdj) == 0) / sum(trueAdj == 0)
  rec[i, 14] <- sum(final_ord$adjMat == trueAdj)
  cat(" 14 ")
  cat("\n")

}

colnames(rec) <- paste(rep(c("gamma", "mixNorm"), each = 7),
                       rep(c("order", "oracle_sen", "oracle_spc", "oracle_edges", "est_sen", "est_spc", "est_edges"), times= 2),
                       sep = "_")

outTab <- data.frame(p, cycleSize, alpha, n, mt, pr, rescaleData, rec)

write.csv(outTab, paste("../results/cycleChain/cycleChain_lowScale_",runInd, ".csv", sep = ""), row.names = F)

