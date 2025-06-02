### Simulations from Section 5 and Appendix B examining performance of the ICA based procedure
#

# runInd is an argument which is passed in. It goes from 1-800 and sets the seed
# and simulation settings below
runInd <- 1
args <- commandArgs(TRUE)
for(i in 1:length(args)){
  eval(parse(text = args[[i]]))
}
print(runInd)

## Folder containing py-tetrad implimentation
setwd("~/py-tetrad/pytetrad")
library("reticulate")
reticulate::use_python(Sys.which("python3"))
source_python("tools/TetradSearch.py")



## the simulation settings considered
p.list <- seq(18, 27, by = 3)
d.list <- c("gamma", "mixedNorm")
sample.size <- 25
n.list <- c(10, 25, 50, 100) * 1000
param.grid <- expand.grid(rep(n.list, sample.size), p.list, d.list)
dim(param.grid)



cycleSize <- 3 # size of cycles
lowEdge <- .5  # lower bound on absolute value of edge
highEdge <- .8 # upper bound on absolute value of edge
lowScale <- .8 # lower bound on sd of idiosyncractic error
highScale <- 1 # upper bound on sd of idiosyncractic  error
parentProb <- 1/2 # probability of edge between ancestor and descendant not in same cycle


n <- param.grid[runInd, 1]
p <- param.grid[runInd, 2]
dist <- as.character(param.grid[runInd, 3])


## Generate Lambda and record adjacency matrix
set.seed(runInd + 1000)
Lambda <- disjointCycles::cycleChain(p, cycleSize = cycleSize, lowEdge = lowEdge, highEdge = highEdge,
                                     parentProb = parentProb)
trueAdj <- (Lambda != 0) +0
trueGraph <- lapply(unname(split(1:p, rep(1:ceiling(p/cycleSize),
                                          each = cycleSize)[1:p])),
                    list)

## Generate data given Lambda
data <- disjointCycles::rLSEM(p, n, dist = dist, LambdaIn = Lambda, lowScale = lowScale, highScale = highScale)


# Scale and center data
Y <- data.frame(scale(data$Y))

# run ICA based procedure from py-tetrad
ts <- TetradSearch(Y)
ts$set_verbose(F)
ts$run_ica_lingd()

# record number of edges (from estimated stable graphs) which do not agree with true adjacency matrix
rec <- data.frame(p, n, cycleSize, dist, edges = sum(t(ts$get_stable_bhats()[[1]] != 0) != trueAdj))

# write.csv(rec, paste("~/cycles/results/icaComp/icaComp_",runInd, ".csv", sep = ""), row.names = F)




#
p.list <- seq(18, 27, by = 3)
d.list <- c("gamma", "mixedNorm")
sample.size <- 50
n.list <- c(10, 25, 50, 1000) * 1000
param.grid <- expand.grid(rep(n.list, sample.size), p.list, d.list)
dim(param.grid)
outTabComp1 <- NULL
missing <- c()

for(runInd in 1:nrow(param.grid)){
  if(!file.exists(paste("~/Dropbox/disjointCycles/simResults/icaComp/icaComp_", runInd, ".csv", sep = ""))){
    missing <- c(missing, runInd)
  } else {
    temp <- read.csv(paste("~/Dropbox/disjointCycles/simResults/icaComp/icaComp_", runInd, ".csv", sep = ""))

    if(is.null(outTabComp1)){
      outTabComp1 <- temp
    } else {
      outTabComp1 <- rbind(outTabComp1, temp)
    }


  }

}
missing
# length(missing)
#
library("tidyverse")
agTab <- outTabComp1 %>% group_by(p, n, dist) %>% summarize(num = length(cycleSize))
agTab[which(agTab$dist == "gamma"), c("p", "n", "num")]%>% xtable::xtable() %>% print(include.rownames = F)
