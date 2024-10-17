runInd <- 1
args <- commandArgs(TRUE)
for(i in 1:length(args)){
  eval(parse(text = args[[i]]))
}
print(runInd)

library("disjointCycles")
p <- 9
n.list <- rep(c(10 , 25, 50) * 1000, each = 8)
rep.runs <- 50
cycleSize <- 3
n <- n.list[runInd]


trueAdj <- matrix(0, p, p)
for(k in 1:3){
  trueAdj[matrix(c(1:cycleSize, 2:cycleSize, 1), byrow = F, ncol = 2) + (k-1) * cycleSize] <- 1
}
trueAdj[1:(p-cycleSize), (p-cycleSize + 1):p] <- 1



rec <- matrix(0, rep.runs, 6)
for(i in 1:rep.runs){

  Lambda <- trueAdj * runif(p^2, .3, .9) * sample(c(-1, 1), size = p^2, replace = T)
  data <- disjointCycles::rLSEM(p, n, dist = "gamma", LambdaIn = Lambda, lowScale = .8)
  Y <- scale(data$Y)

  maxCliques <- list(c(1:3), c(4:6), c(7:9))

  for(k in 1:length(maxCliques)){
    C <- maxCliques[[k]]
    D <- setdiff(1:p, C)
    rec[i, k] <- pruneHelperRootCycle(C, D, Y, methodPR = "infFunc")
    rec[i, k + 3] <- pruneHelperRootCycle(C, D, Y, methodPR = "none")
  }

  cat(i)
  cat("\n")
}

write.csv(rec, paste("~/cycles/results/pruneRoots/pruneRoots_", runInd, ".csv"))

# alpha <- .05
# colMeans(rec < alpha)
# par(mfcol = c(3, 2))
# hist(rec[, 1], main = "Null Inf Func")
# mtext(paste("Reject:", mean(rec[, 1] < alpha), sep = " "))
# hist(rec[, 2], main = "Null Inf Func")
# mtext(paste("Reject:", mean(rec[, 2] < alpha), sep = " "))
# hist(rec[, 3], main = "Alt Inf Func")
# mtext(paste("Reject:", mean(rec[, 3] < alpha), sep = " "))
#
# hist(rec[, 4], main = "Null Naive")
# mtext(paste("Reject:", mean(rec[, 4] < alpha), sep = " "))
# hist(rec[, 5], main = "Null Naive")
# mtext(paste("Reject:", mean(rec[, 5] < alpha), sep = " "))
# hist(rec[, 6], main = "Alt Naive")
# mtext(paste("Reject:", mean(rec[, 6] < alpha), sep = " "))
#
