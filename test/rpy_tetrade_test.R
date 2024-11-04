runInd <- 1
args <- commandArgs(TRUE)
for(i in 1:length(args)){
  eval(parse(text = args[[i]]))
}
print(runInd)


setwd("~/py-tetrad/pytetrad")
library("reticulate")
reticulate::use_python(Sys.which("python3"))

source_python("tools/TetradSearch.py")


p <- 12
n <- 50000
cycleSize <- 6

set.seed(runInd + 1000)
Lambda <- disjointCycles::cycleChain(p = p, cycleSize = cycleSize, lowEdge = .5, highEdge = .9,
                                     parentProb = 1/2)
trueAdj <- (Lambda != 0) +0

data <- disjointCycles::rLSEM(p, n, dist = "gamma",
                              LambdaIn = Lambda, lowScale = .8)
Y <- data.frame(scale(data$Y))

ts <- TetradSearch(Y)
ts$set_verbose(F)
ts$run_ica_lingd()

sum(t(ts$get_stable_bhats()[[1]] != 0) != trueAdj)


rec <- matrix(0, sim.size, 3)
for(i in 1:sim.size){
  Lambda <- disjointCycles::cycleChain(p = p, cycleSize = cycleSize, lowEdge = .1, highEdge = .5,
                                       parentProb = 1)
  trueAdj <- (Lambda != 0) +0

  data <- disjointCycles::rLSEM(p, n, dist = "gamma",
                                LambdaIn = Lambda, lowScale = .8)
  Y <- data.frame(scale(data$Y))

  ts <- TetradSearch(Y)
  ts$set_verbose(F)
  ts$run_ica_lingd()

  rec[i, 1] <- sum(t(ts$get_stable_bhats()[[1]] != 0) != trueAdj)


  data <- disjointCycles::rLSEM(p, n, dist = "mixedNorm",
                                LambdaIn = Lambda, lowScale = .8)
  Y <- data.frame(scale(data$Y))

  ts <- TetradSearch(Y)
  ts$set_verbose(FALSE)
  ts$run_ica_lingd()

  rec[i, 2] <- sum(t(ts$get_stable_bhats()[[1]] != 0) != trueAdj)



  data <- disjointCycles::rLSEM(p, n, dist = "beta",
                                LambdaIn = Lambda, lowScale = .8)

  Y <- data.frame(scale(data$Y))
  ts <- TetradSearch(Y)
  ts$set_verbose(FALSE)
  ts$run_ica_lingd()

  rec[i, 3] <- sum(t(ts$get_stable_bhats()[[1]] != 0) != trueAdj)

  cat(i)
  cat(": ")
  cat(colMeans(rec[1:i, , drop = F]))
  cat("\n")
}

colMeans(rec)


