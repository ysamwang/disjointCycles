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
p.list <- seq(24, 30, by = 3)
d.list <- c("gamma", "mixedNorm")
sample.size <- 100
n.list <- c(5, 25) * 1000
param.grid <- expand.grid(rep(n.list, sample.size), p.list, d.list)
dim(param.grid)


cycleSize <- 3
lowEdge <- .5
highEdge <- .8
lowScale <- .8
highScale <- 1
parentProb <- 3/4

n <- param.grid[runInd, 1]
p <- param.grid[runInd, 2]
dist <- as.character(param.grid[runInd, 3])


set.seed(runInd + 1000)
Lambda <- disjointCycles::cycleChain(p, cycleSize = cycleSize, lowEdge = lowEdge, highEdge = highEdge,
                                     parentProb = parentProb)

trueAdj <- (Lambda != 0) +0
trueGraph <- lapply(unname(split(1:p, rep(1:ceiling(p/cycleSize),
                                          each = cycleSize)[1:p])),
                    list)
data <- disjointCycles::rLSEM(p, n, dist = dist, LambdaIn = Lambda, lowScale = lowScale, highScale = highScale)



Y <- data.frame(scale(data$Y))

ts <- TetradSearch(Y)
ts$set_verbose(F)
ts$run_ica_lingd()

rec <- data.frame(p, n, cycleSize, edges = sum(t(ts$get_stable_bhats()[[1]] != 0) != trueAdj))

write.csv(rec, paste("~/cycles/results/icaComp/icaComp_",runInd, ".csv", sep = ""), row.names = F)



