runInd <- 1
args <- commandArgs(TRUE)
for(i in 1:length(args)){
  eval(parse(text = args[[i]]))
}
print(runInd)


runInd <- runInd + 3000

set.seed(runInd + 2000)
library("disjointCycles")
sample.size <- 50
rep.runs <- 1
# p.list <- c(20, 25, 30)
p.list <- seq(30, 45, by = 15)
n.list <- c(50, 100) * 1000
a.list <- c(3e-1, 5e-2, 1e-2, 1e-3)
mt.list <- c("BH", "holm")
pr.list <- c("naive", "chisq")
c.list <- c(3)
res.list <- c(T)
param.grid <- expand.grid(rep(n.list, sample.size / rep.runs), a.list, mt.list, pr.list, c.list, res.list, p.list)
dim(param.grid)
# 900

lowEdge <- .5
highEdge <- .8
lowScale <- .8
highScale <- 1
parentProb <- 3/4

n <- param.grid[runInd, 1]
alpha <- param.grid[runInd, 2]
mt <- as.character(param.grid[runInd, 3])
pr <- as.character(param.grid[runInd, 4])
cycleSize <- param.grid[runInd, 5]
rescaleData <- param.grid[runInd, 6]
p <- param.grid[runInd, 7]

rec <- matrix(0, rep.runs, 18)
for(i in 1:rep.runs){
  Lambda <- disjointCycles::cycleChain(p, cycleSize = cycleSize, lowEdge = lowEdge, highEdge = highEdge,
                                       parentProb = parentProb)

  trueAdj <- (Lambda != 0) +0
  trueGraph <- lapply(unname(split(1:p, rep(1:ceiling(p/cycleSize),
                                            each = cycleSize)[1:p])),
                      list)
  cat("Iter ")
  cat(i)
  data <- rLSEM(p, n, dist = "gamma", LambdaIn = Lambda, lowScale = lowScale, highScale = highScale)

  time1 <- system.time(est_ord <- djcGetOrderNew(data$Y, verbose = F, alpha2 = alpha,
                                                 alpha3 = alpha, alphaR = alpha, pvalAdjMethod = mt,
                                                 methodPR = pr,
                                                 rescaleData = rescaleData))

  rec[i, 1] <- compareOrders(trueGraph, est_ord)
  # ## Estimate edges given true layers
  est_Edges <- djcGetEdges(trueGraph, data$Y, alpha = alpha,pvalAdjMethod = mt)
  rec[i, 2] <- sum(est_Edges$adjMat *trueAdj) / sum(trueAdj)
  rec[i, 3] <- sum((est_Edges$adjMat + trueAdj) == 0) / sum(trueAdj == 0)
  rec[i, 4] <- sum(est_Edges$adjMat == trueAdj)

  ## Estimate edges given estimated layers
  time2 <- system.time(final_ord <- djcGetEdges(est_ord, data$Y, alpha = alpha,pvalAdjMethod = mt))
  rec[i, 5] <- sum(final_ord$adjMat *trueAdj) / sum(trueAdj)
  rec[i, 6] <- sum((final_ord$adjMat + trueAdj) == 0) / sum(trueAdj == 0)
  rec[i, 7] <- sum(final_ord$adjMat == trueAdj)

  rec[i, 8] <- time1[3]
  rec[i, 9] <- time2[3]


  cat(": gamma done; ")
#
#   data <- rLSEM(p, n, dist = "mixedNorm", LambdaIn = Lambda, lowScale = lowScale, highScale = highScale)
#   time1 <- system.time(est_ord <- djcGetOrderNew(data$Y, verbose = F, alpha2 = alpha,
#                                                  alpha3 = alpha, alphaR = alpha, pvalAdjMethod = mt, methodPR = pr,
#                                                  rescaleData = rescaleData))
#
#   rec[i, 10] <- compareOrders(trueGraph, est_ord)
#
#   # ## Estimate edges given true layers
#   est_Edges <- djcGetEdges(trueGraph, data$Y,
#                            alpha = alpha, pvalAdjMethod = mt)
#
#   rec[i, 11] <- sum(est_Edges$adjMat *trueAdj) / sum(trueAdj)
#   rec[i, 12] <- sum((est_Edges$adjMat + trueAdj) == 0) / sum(trueAdj == 0)
#   rec[i, 13] <- sum(est_Edges$adjMat == trueAdj)
#   #
#   #
#   # ## Estimate edges given true layers
#   time2 <- system.time(final_ord <- djcGetEdges(est_ord, data$Y,
#                                                 alpha = alpha, pvalAdjMethod = mt))
#
#   rec[i, 14] <- sum(final_ord$adjMat *trueAdj) / sum(trueAdj)
#   rec[i, 15] <- sum((final_ord$adjMat + trueAdj) == 0) / sum(trueAdj == 0)
#   rec[i, 16] <- sum(final_ord$adjMat == trueAdj)
#
#   rec[i, 17] <- time1[3]
#   rec[i, 18] <- time2[3]
#
#   cat(": mixNorm done; ")

  # data <- rLSEM(p, n, dist = "lognorm", LambdaIn = Lambda, lowScale = .8, highScale = 1)
  # est_ord <- djcGetOrderNew(data$Y, verbose = F, alpha2 = alpha,
  #                           alpha3 = alpha, alphaR = alpha, pvalAdjMethod = mt, methodPR = pr,
  #                           rescaleData = rescaleData)
  # rec[i, 15] <- compareOrders(trueGraph, est_ord)
  #
  # # ## Estimate edges given true layers
  # est_Edges <- djcGetEdges(trueGraph, data$Y,
  #                          alpha = alpha, pvalAdjMethod = mt)
  #
  # rec[i, 16] <- sum(est_Edges$adjMat *trueAdj) / sum(trueAdj)
  # rec[i, 17] <- sum((est_Edges$adjMat + trueAdj) == 0) / sum(trueAdj == 0)
  # rec[i, 18] <- sum(est_Edges$adjMat == trueAdj)
  # #
  # #
  # # ## Estimate edges given true layers
  # final_ord <- djcGetEdges(est_ord, data$Y,
  #                          alpha = alpha, pvalAdjMethod = mt)
  #
  # rec[i, 19] <- sum(final_ord$adjMat *trueAdj) / sum(trueAdj)
  # rec[i, 20] <- sum((final_ord$adjMat + trueAdj) == 0) / sum(trueAdj == 0)
  # rec[i, 21] <- sum(final_ord$adjMat == trueAdj)
  # cat(": lognormal done")



}

colnames(rec) <- paste(rep(c("gamma", "mixNorm"), each = 9),
                       rep(c("order", "oracle_sen", "oracle_spc", "oracle_edges",
                             "est_sen", "est_spc", "est_edges", "timeOrd", "timeEdges"), times= 2),
                       sep = "_")

outTab <- data.frame(p, cycleSize, alpha, n, mt, pr, rescaleData, rec)

write.csv(outTab, paste("../results/large/large1_",runInd, ".csv", sep = ""), row.names = F)

#
#
#
# library("disjointCycles")
# library("tidyverse")
# sample.size <- 200
# rep.runs <- 10
# n.list <- c(25, 50) * 1000
# a.list <- c(2.5e-1, 1e-1, 1e-2, 1e-3, 1e-4)
# mt.list <- c("BH", "holm")
# pr.list <- c("naive", "chisq")
# c.list <- c(6)
# res.list <- c(T)
# param.grid <- expand.grid(rep(n.list, sample.size / rep.runs), a.list, mt.list, pr.list, c.list, res.list)
# dim(param.grid)
#
# runInd <- 1
# outTabComp <- read.csv(paste("~/Dropbox/disjointCycles/simResults/large/large_",runInd, ".csv", sep = ""),)
# missing <- c()
# for(runInd in 2:nrow(param.grid)){
#   if(!file.exists(paste("~/Dropbox/disjointCycles/simResults/large/large_", runInd, ".csv", sep = ""))){
#     missing <- c(missing, runInd)
#   } else {
#     temp <- read.csv(paste("~/Dropbox/disjointCycles/simResults/large/large_", runInd, ".csv", sep = ""))
#     outTabComp <- rbind(outTabComp, temp)
#   }
#
# }
# missing
#
#
# agTab <- outTabComp %>% group_by(p, n, alpha, mt, pr, rescaleData) %>%
#   summarize(g_ord = mean(gamma_order), m_ord = mean(mixNorm_order),
#             ln_ord = mean(ln_order),
#             g_edge_est = mean(gamma_est_edges), m_edge_est = mean(mixNorm_est_edges),
#             ln_edge_est = mean(ln_est_edges),
#             g_edge_orc = mean(gamma_oracle_edges), m_edge_orc = mean(mixNorm_oracle_edges),
#             ln_edge_orc = mean(ln_oracle_edges))
#
# agTab %>% filter(n == 50000, mt == "holm")
#
# agTabMax <- agTab %>% filter(rescaleData == T) %>%group_by(p, n, pr) %>%
#   summarize(g_ord = max(g_ord), m_ord = max(m_ord), ln_ord = max(ln_ord),
#             g_edge_est = max(g_edge_est), m_edge_est = max(m_edge_est), ln_edge_est = max(ln_edge_est),
#             g_edge_orc = max(g_edge_orc), m_edge_orc = max(m_edge_orc), ln_edge_orc = max(ln_edge_orc)) %>% print(n = 50)
#
library("disjointCycles")
library("tidyverse")
sample.size <- 20
rep.runs <- 1
# p.list <- c(20, 25, 30)
p.list <- seq(30, 50, by = 10)
n.list <- c(50, 100) * 1000
a.list <- c(2.5e-1, 1e-1, 1e-2, 1e-4)
mt.list <- c("BH", "holm")
pr.list <- c("naive")
c.list <- c(5)
res.list <- c(T)
param.grid <- expand.grid(rep(n.list, sample.size / rep.runs), a.list, mt.list, pr.list, c.list, res.list, p.list)
dim(param.grid)


runInd <- 1
outTabComp1 <- read.csv(paste("~/Dropbox/disjointCycles/simResults/large/large1_",runInd, ".csv", sep = ""),)
missing <- c()
for(runInd in 2:nrow(param.grid)){
  if(!file.exists(paste("~/Dropbox/disjointCycles/simResults/large/large1_", runInd, ".csv", sep = ""))){
    missing <- c(missing, runInd)
  } else {
    temp <- read.csv(paste("~/Dropbox/disjointCycles/simResults/large/large1_", runInd, ".csv", sep = ""))
    outTabComp1 <- rbind(outTabComp1, temp)
  }

}
missing



agTab1 <- outTabComp1 %>% group_by(p, n, mt, alpha, pr, rescaleData) %>%
  summarize(g_ord = mean(gamma_order),
            g_edge_est = mean((gamma_est_edges -p)/(p*(p-1))))




agTabMax1 <- agTab1 %>%  group_by(p, n, mt) %>%
  summarize(g_ord = max(g_ord),
            g_edge_est = max(g_edge_est)) %>% print(n = 50)


agTabMax1

