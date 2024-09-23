runInd <- 1
args <- commandArgs(TRUE)
for(i in 1:length(args)){
  eval(parse(text = args[[i]]))
}
print(runInd)


library("disjointCycles")
sample.size <- 400
rep.runs <- 10
n.list <- c(10, 25, 50) * 1000
a.list <- c(.2, .1, .05, .01, .001)
param.grid <- expand.grid(rep(n.list, sample.size / rep.runs), a.list)
# 600



p <- 8
n <- param.grid[runInd, 1]
alpha <- param.grid[runInd, 2]

rec <- matrix(0, rep.runs, 6)

for(i in 1:rep.runs){
  Lambda <- ex3_2(lowEdge = .5, highEdge = .9)
  trueAdj <- (Lambda != 0) + 0

  trueGraph <- list(list(c(1), c(2)), list(c(3,4,5)), list(c(6,7,8)))

  data <- rLSEM(p, n, dist = "gamma", LambdaIn = Lambda, lowScale = .8)
  est_ord <- djcGetOrder(scale(data$Y), verbose = F, alpha2 = alpha, alpha3 = alpha, alphaR = alpha)
  rec[i, 1] <- compareOrders(trueGraph, est_ord)

  ## Estimate edges given true layers
  est_Edges <- djcGetEdges(trueGraph, data$Y, alpha = alpha)
  rec[i, 2] <- all(est_Edges$adjMat == trueAdj)

  ## Estimate edges given estimated layers
  final_ord <- djcGetEdges(est_ord, data$Y, alpha = alpha)
  rec[i, 3] <- all(final_ord$adjMat == trueAdj)


  data <- rLSEM(p, n, dist = "mixedNorm", LambdaIn = Lambda, lowScale = .8)
  est_ord <- djcGetOrder(scale(data$Y), verbose = F, alpha2 = alpha, alpha3 = alpha, alphaR = alpha)
  rec[i, 4] <- compareOrders(trueGraph, est_ord)

  ## Estimate edges given true layers
  est_Edges <- djcGetEdges(trueGraph, data$Y, alpha = alpha)
  rec[i, 5] <- all(est_Edges$adjMat == trueAdj)

  ## Estimate edges given true layers
  final_ord <- djcGetEdges(est_ord, data$Y, alpha = alpha)
  rec[i, 6] <- all(final_ord$adjMat == trueAdj)

}

colnames(rec) <- paste(rep(c("gamma", "mixNorm"), each = 3), rep(c("order", "edges", "total"), times= 2), sep = "_")
outTab <- data.frame(alpha, n, rec)

write.csv(outTab, paste("../results/ex32/ex32_",runInd, ".csv", sep = ""))

