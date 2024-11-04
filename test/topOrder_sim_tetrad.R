runInd <- 1
args <- commandArgs(TRUE)
for(i in 1:length(args)){
  eval(parse(text = args[[i]]))
}
print(runInd)


set.seed(runInd + 1000)
library("disjointCycles")
library("reticulate")
reticulate::use_python(Sys.which("python3"))
source_python("tools/TetradSearch.py")


sample.size <- 400
rep.runs <- 20
n.list <- c(10, 25, 50) * 1000
a.list <- c(2.5e-1, 1e-1, 1e-2, 1e-3, 1e-4)
mt.list <- c("BH", "holm")
pr.list <- c("naive")
c.list <- c(2,3, 4)
res.list <- c(T)
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

rec <- matrix(0, rep.runs, 3)
for(i in 1:rep.runs){
  Lambda <- disjointCycles::cycleChain(p, cycleSize = cycleSize, lowEdge = .5, highEdge = .9,
                                       parentProb = 1/2)

  trueAdj <- (Lambda != 0) +0
  trueGraph <- lapply(unname(split(1:p, rep(1:ceiling(p/cycleSize),
                                                          each = cycleSize)[1:p])),
                                    list)
  cat("Iter ")
  cat(i)
  data <- disjointCycles::rLSEM(p, n, dist = "gamma", LambdaIn = Lambda, lowScale = .8)
  Y <- data.frame(scale(data$Y))
  ts <- TetradSearch(Y)
  ts$run_ica_lingd()

  ## Estimate edges given estimated layers
  rec[i, 1] <- sum(t(ts$get_stable_bhats()[[1]]) == trueAdj)

  cat(": gamma done; ")

  data <- rLSEM(p, n, dist = "mixedNorm", LambdaIn = Lambda, lowScale = .8)
  Y <- data.frame(scale(data$Y))
  ts <- TetradSearch(Y)
  ts$run_ica_lingd()

  ## Estimate edges given estimated layers
  rec[i, 2] <- sum(t(ts$get_stable_bhats()[[1]]) == trueAdj)

  data <- rLSEM(p, n, dist = "lognorm", LambdaIn = Lambda, lowScale = .8, highScale = 1)
  Y <- data.frame(scale(data$Y))
  ts <- TetradSearch(Y)
  ts$run_ica_lingd()

  ## Estimate edges given estimated layers
  rec[i, 3] <- sum(t(ts$get_stable_bhats()[[1]]) == trueAdj)
  cat(": lognormal done\n")
}

colnames(rec) <- paste(c("gamma", "mixNorm", "ln"), "est_edges", sep = "_")


outTab <- data.frame(p, cycleSize, n, rec)

write.csv(outTab, paste("../results/cycleChain/cycleChain_ica_",runInd, ".csv", sep = ""), row.names = F)

