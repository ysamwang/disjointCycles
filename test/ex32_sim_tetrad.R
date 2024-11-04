runInd <- 1
args <- commandArgs(TRUE)
for(i in 1:length(args)){
  eval(parse(text = args[[i]]))
}
print(runInd)

set.seed(runInd + 3000)
setwd("~/py-tetrad/pytetrad")
library("disjointCycles")
library("reticulate")
reticulate::use_python(Sys.which("python3"))
source_python("tools/TetradSearch.py")


sample.size <- 400
rep.runs <- 1
n.list <- c(10, 25, 50) * 1000
param.grid <- expand.grid(rep(n.list, sample.size / rep.runs))
dim(param.grid)
# 600


p <- 8
n <- param.grid[runInd, 1]

rec <- matrix(0, rep.runs, 3)

for(i in 1:rep.runs){

  Lambda <- disjointCycles::ex3_2(lowEdge = .5, highEdge = .9)
  trueAdj <- (Lambda != 0) + 0

  trueGraph <- list(list(c(1), c(2)), list(c(3,4,5)), list(c(6,7,8)))

  cat("Iter ")
  cat(i)

  data <- disjointCycles::rLSEM(p, n, dist = "gamma",
                                LambdaIn = Lambda, lowScale = .8)

  Y <- data.frame(scale(data$Y))
  ts <- TetradSearch(Y)
  ts$run_ica_lingd()

  ## Estimate edges given estimated layers
  rec[i, 1] <- sum(t(ts$get_stable_bhats()[[1]]!=0) == trueAdj)

  cat(": gamma done; ")

  data <- rLSEM(p, n, dist = "mixedNorm", LambdaIn = Lambda, lowScale = .8)
  Y <- data.frame(scale(data$Y))
  ts <- TetradSearch(Y)
  ts$run_ica_lingd()

  ## Estimate edges given estimated layers
  rec[i, 2] <- sum(t(ts$get_stable_bhats()[[1]]!=0) == trueAdj)

  data <- rLSEM(p, n, dist = "lognorm", LambdaIn = Lambda, lowScale = .8, highScale = 1)
  Y <- data.frame(scale(data$Y))
  ts <- TetradSearch(Y)
  ts$run_ica_lingd()

  ## Estimate edges given estimated layers
  rec[i, 3] <- sum(t(ts$get_stable_bhats()[[1]]!=0) == trueAdj)
  cat(": lognormal done\n")


}

colnames(rec) <- paste(c("gamma", "mixNorm", "ln"), "est_edges", sep = "_")
outTab <- data.frame(p, n, rec)

write.csv(outTab, paste("../results/ex32/ex32_ica_",runInd, ".csv", sep = ""), row.names = F)

