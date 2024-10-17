runInd <- 1
args <- commandArgs(TRUE)
for(i in 1:length(args)){
  eval(parse(text = args[[i]]))
}
print(runInd)

library("disjointCycles")
sample.size <- 100
rep.runs <- 10
n.list <- c(10, 25, 50) * 1000
a.list <- c(1e-1, 1e-2, 1e-3)
mt.list <- c("BH", "holm")
pr.list <- c("naive")
c.list <- c(2, 3, 4)
res.list <- c(T)
param.grid <- expand.grid(rep(n.list, sample.size / rep.runs), a.list, mt.list, pr.list, c.list, res.list)
dim(param.grid)
# 960

p <- 12
n <- param.grid[runInd, 1]
alpha <- param.grid[runInd, 2]
mt <- as.character(param.grid[runInd, 3])
pr <- as.character(param.grid[runInd, 4])
cycleSize <- param.grid[runInd, 5]
rescaleData <- param.grid[runInd, 6]

rec <- matrix(0, rep.runs, 4)

for(i in 1:rep.runs){
  Lambda <- disjointCycles::cycleChain(p, cycleSize = cycleSize, lowEdge = .3,
                                       highEdge = .9, parentProb = 1/3)

  trueAdj <- (Lambda != 0) +0
  trueGraph <- lapply(unname(split(1:p, rep(1:ceiling(p/cycleSize),
                                                          each = cycleSize)[1:p])), list)

  data <- disjointCycles::rLSEM(p, n, dist = "gamma", LambdaIn = Lambda, lowScale = .6, highScale = .95)
  data$sigma

  est_ord <- djcGetOrder(data$Y, verbose = F, alpha2 = alpha,
                         alpha3 = alpha, alphaR = alpha, pvalAdjMethod = mt,
                         methodPR = pr, rescaleData = rescaleData)
  rec[i, 1] <- compareOrders(trueGraph, est_ord)

  est_ord <- djcGetOrder(data$Y, verbose = F, alpha2 = alpha,
                         alpha3 = alpha, alphaR = alpha, pvalAdjMethod = mt,
                         methodPR = pr, rescaleData = rescaleData, sigmaPop = data$sigma)
  rec[i, 2] <- compareOrders(trueGraph, est_ord)


  data <- rLSEM(p, n, dist = "mixedNorm", LambdaIn = Lambda, lowScale = .8)
  est_ord <- djcGetOrder(data$Y, verbose = F, alpha2 = alpha,
                         alpha3 = alpha, alphaR = alpha, pvalAdjMethod = mt,
                         methodPR = pr, rescaleData = rescaleData)
  rec[i, 3] <- compareOrders(trueGraph, est_ord)

  est_ord <- djcGetOrder(data$Y, verbose = F, alpha2 = alpha,
                         alpha3 = alpha, alphaR = alpha, pvalAdjMethod = mt,
                         methodPR = pr, rescaleData = rescaleData, sigmaPop = data$sigma)
  rec[i, 4] <- compareOrders(trueGraph, est_ord)




}

colnames(rec) <- paste(rep(c("gamma", "mixNorm"), each = 2), rep(c("_order","_order_oracle"), times = 2), sep = "")
# colnames(rec) <- paste(c("gamma", "mixNorm"), "_order", sep = "")

outTab <- data.frame(cycleSize, alpha, n, mt, pr, rescaleData, rec)

write.csv(outTab, paste("../results/cycleChain/cycleChain_oracle_",runInd, ".csv", sep = ""), row.names = F)

