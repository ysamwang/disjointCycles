runInd <- 1
args <- commandArgs(TRUE)
for(i in 1:length(args)){
  eval(parse(text = args[[i]]))
}
print(runInd)


set.seed(runInd + 3000)
library("disjointCycles")
sample.size <- 400
rep.runs <- 20
n.list <- c(10, 25, 50) * 1000
a.list <- c(2.5e-1, 1e-1, 1e-2, 1e-3, 1e-4)
mt.list <- c("BH", "holm")
pr.list <- c("naive")
res.list <- c(T)
param.grid <- expand.grid(rep(n.list, sample.size / rep.runs), a.list, mt.list, pr.list, res.list)
dim(param.grid)
# 600



p <- 8
n <- param.grid[runInd, 1]
alpha <- param.grid[runInd, 2]
mt <- as.character(param.grid[runInd, 3])
pr <- as.character(param.grid[runInd, 4])
rescaleData <- param.grid[runInd, 5]

rec <- matrix(0, rep.runs, 21)

for(i in 1:rep.runs){
  Lambda <- disjointCycles::ex3_2(lowEdge = .5, highEdge = .9)
  trueAdj <- (Lambda != 0) + 0

  trueGraph <- list(list(c(1), c(2)), list(c(3,4,5)), list(c(6,7,8)))

  cat("Iter ")
  cat(i)
  data <- disjointCycles::rLSEM(p, n, dist = "gamma", LambdaIn = Lambda, lowScale = .8)

  est_ord <- djcGetOrderNew(data$Y, verbose = F, alpha2 = alpha,
                            alpha3 = alpha, alphaR = alpha, pvalAdjMethod = mt,
                            methodPR = pr,
                            rescaleData = rescaleData)
  rec[i, 1] <- compareOrders(trueGraph, est_ord)
  # ## Estimate edges given true layers
  # est_Edges <- djcGetEdges(trueGraph, data$Y, alpha = alpha,pvalAdjMethod = mt)
  # rec[i, 2] <- sum(est_Edges$adjMat *trueAdj) / sum(trueAdj)
  # rec[i, 3] <- sum((est_Edges$adjMat + trueAdj) == 0) / sum(trueAdj == 0)
  # rec[i, 4] <- sum(est_Edges$adjMat == trueAdj)

  ## Estimate edges given estimated layers
  final_ord <- djcGetEdges(est_ord, data$Y, alpha = alpha, pvalAdjMethod = mt, rescaleData = T)
  rec[i, 5] <- sum(final_ord$adjMat *trueAdj) / sum(trueAdj)
  rec[i, 6] <- sum((final_ord$adjMat + trueAdj) == 0) / sum(trueAdj == 0)
  rec[i, 7] <- sum(final_ord$adjMat == trueAdj)

  cat(": gamma done; ")

  data <- rLSEM(p, n, dist = "mixedNorm", LambdaIn = Lambda, lowScale = .8)
  est_ord <- djcGetOrderNew(data$Y, verbose = F, alpha2 = alpha,
                            alpha3 = alpha, alphaR = alpha, pvalAdjMethod = mt, methodPR = pr,
                            rescaleData = rescaleData)
  rec[i, 8] <- compareOrders(trueGraph, est_ord)

  # # ## Estimate edges given true layers
  # est_Edges <- djcGetEdges(trueGraph, data$Y,
  #                          alpha = alpha, pvalAdjMethod = mt)
  #
  # rec[i, 9] <- sum(est_Edges$adjMat *trueAdj) / sum(trueAdj)
  # rec[i, 10] <- sum((est_Edges$adjMat + trueAdj) == 0) / sum(trueAdj == 0)
  # rec[i, 11] <- sum(est_Edges$adjMat == trueAdj)
  #
  #
  # ## Estimate edges given true layers
  final_ord <- djcGetEdges(est_ord, data$Y,
                           alpha = alpha, pvalAdjMethod = mt, rescaleData = T)

  rec[i, 12] <- sum(final_ord$adjMat *trueAdj) / sum(trueAdj)
  rec[i, 13] <- sum((final_ord$adjMat + trueAdj) == 0) / sum(trueAdj == 0)
  rec[i, 14] <- sum(final_ord$adjMat == trueAdj)
  cat(": mixNorm done; ")

  data <- rLSEM(p, n, dist = "lognorm", LambdaIn = Lambda, lowScale = .8, highScale = 1)
  est_ord <- djcGetOrderNew(data$Y, verbose = F, alpha2 = alpha,
                            alpha3 = alpha, alphaR = alpha, pvalAdjMethod = mt, methodPR = pr,
                            rescaleData = rescaleData)
  rec[i, 15] <- compareOrders(trueGraph, est_ord)

  # # ## Estimate edges given true layers
  # est_Edges <- djcGetEdges(trueGraph, data$Y,
  #                          alpha = alpha, pvalAdjMethod = mt)
  #
  # rec[i, 16] <- sum(est_Edges$adjMat *trueAdj) / sum(trueAdj)
  # rec[i, 17] <- sum((est_Edges$adjMat + trueAdj) == 0) / sum(trueAdj == 0)
  # rec[i, 18] <- sum(est_Edges$adjMat == trueAdj)
  #
  #
  # ## Estimate edges given true layers
  final_ord <- djcGetEdges(est_ord, data$Y,
                           alpha = alpha, pvalAdjMethod = mt, rescaleData = T)

  rec[i, 19] <- sum(final_ord$adjMat *trueAdj) / sum(trueAdj)
  rec[i, 20] <- sum((final_ord$adjMat + trueAdj) == 0) / sum(trueAdj == 0)
  rec[i, 21] <- sum(final_ord$adjMat == trueAdj)
  cat(": lognormal done")


}

colnames(rec) <- paste(rep(c("gamma", "mixNorm", "ln"), each = 7),
                       rep(c("order", "oracle_sen", "oracle_spc", "oracle_edges",
                             "est_sen", "est_spc", "est_edges"), times= 3),
                       sep = "_")
outTab <- data.frame(p, alpha, n, mt, pr, rescaleData, rec)

write.csv(outTab, paste("../results/ex32/ex32_",runInd, ".csv", sep = ""), row.names = F)

