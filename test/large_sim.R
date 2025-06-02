### Simulations from Section 5 and Appendix B examining performance of the proposed procedure

# runInd is an argument which is passed in. It goes from 1-800 and sets the seed
# and simulation settings below
runInd <- 1
args <- commandArgs(TRUE)
for(i in 1:length(args)){
  eval(parse(text = args[[i]]))
}
print(runInd)

set.seed(runInd + 2000)
library("disjointCycles")


## the simulation settings considered
sample.size <- 100
rep.runs <- 2
p.list <- c(30, 45, 60)
n.list <- c(10, 25, 50, 100) * 1000
a.list <- c(2e-1, 5e-2, 1e-2, 1e-3) # the nominal level of each test conducted
mt.list <- c("BH", "holm") # the multiple testing adjustment procedure used
pr.list <- c("chisq") # the test used for pruning roots
c.list <- c(3) # size of cycles
# c.list <- c(5)
res.list <- c(T)
ep.list <- c(1/2) #probability of edge from ancestor to descendant not in same cycle
param.grid <- expand.grid(rep(n.list, sample.size / rep.runs), a.list, mt.list,
                          pr.list, c.list, res.list, ep.list)
dim(param.grid)
# 900


lowEdge <- .5  # lower bound on absolute value of edge
highEdge <- .8 # upper bound on absolute value of edge
lowScale <- .8 # lower bound on sd of idiosyncractic error
highScale <- 1 # upper bound on sd of idiosyncractic  error


n <- param.grid[runInd, 1]
alpha <- param.grid[runInd, 2] # the nominal level of each test conducted
mt <- as.character(param.grid[runInd, 3]) # the multiple testing adjustment procedure used
pr <- as.character(param.grid[runInd, 4]) # the test used for pruning roots
cycleSize <- param.grid[runInd, 5] # size of each cycle
rescaleData <- param.grid[runInd, 6] # rescale and center the data
parentProb <- param.grid[runInd, 7] #probability of edge from ancestor to descendant not in same cycle

# distributions to consider
dist.list <- c("mixedNorm", "gamma")

rec <- data.frame(matrix(0, rep.runs * length(p.list) * length(dist.list), 7))
recInd <- 0

for(i in 1:rep.runs){

  cat("Iter")
  cat(i)
  cat("\n")

  for(pInd in 1:length(p.list)){

    p <- p.list[pInd]

    cat("p: ")
    cat(p)
    cat("\n")

    ## Generate Lambda and record adjacency matrix
    Lambda <- disjointCycles::cycleChain(p, cycleSize = cycleSize,
                                         lowEdge = lowEdge, highEdge = highEdge,
                                         parentProb = parentProb)
    trueAdj <- (Lambda != 0) +0
    trueGraph <- lapply(unname(split(1:p, rep(1:ceiling(p/cycleSize),
                                              each = cycleSize)[1:p])), list)



    # sample data from each type of distribution in d.list
    for(distInd in 1:length(dist.list)){

        recInd <- recInd + 1
        cat(dist.list[distInd])

        # sample data
        data <- rLSEM(p, n, dist = dist.list[distInd], LambdaIn = Lambda,
                      lowScale = lowScale, highScale = highScale)

        # estimate ordering for the graph
        time1 <- system.time(est_ord <-
                               disjointCycles::djcGetOrderNew(data$Y,
                                                              verbose = T, alpha2 = alpha,
                                                              alpha3 = alpha, alphaR = alpha, pvalAdjMethod = mt,
                                                              methodPR = pr,
                                                              rescaleData = rescaleData))

        rec[recInd, 1] <- p
        rec[recInd, 2] <- as.character(dist.list[distInd])

        rec[recInd, 3] <- compareOrders(trueGraph, est_ord)

        # Estimate edges given estimated layers
        time2 <- system.time(est_Edges <-
                               djcGetEdges(est_ord, data$Y, alpha = alpha,
                                           pvalAdjMethod = mt))

        # True positives
        rec[recInd, 4] <- (sum(est_Edges$adjMat *trueAdj) - p) / (sum(trueAdj) -p)
        # True negatives
        rec[recInd, 5] <- (sum((est_Edges$adjMat + trueAdj) == 0) - p) / (sum(trueAdj == 0) -p)
        # correct edges
        rec[recInd, 6] <- sum(est_Edges$adjMat == trueAdj)

        # time for each step
        rec[recInd, 7] <- time1[3]
        rec[recInd, 8] <- time2[3]

        }
    }
}

colnames(rec) <-c("p", "dist", "order", "est_sen", "est_spc", "est_edges", "timeOrd", "timeEdges")

outTab <- data.frame(cycleSize, alpha, n, mt, pr, rescaleData, rec)
#
# write.csv(outTab, paste("../results/large/large_",runInd, ".csv", sep = ""), row.names = F)

## Create Plots
library("disjointCycles")
library("tidyverse")
sample.size <- 100
rep.runs <- 2
p.list <- c(30, 45, 60)
n.list <- c(10, 25, 50, 100) * 1000
a.list <- c(2e-1, 5e-2, 1e-2, 1e-3)
mt.list <- c("BH", "holm")
pr.list <- c("chisq")
c.list <- c(3)
res.list <- c(T)
ep.list <- c(1/2)
param.grid <- expand.grid(rep(n.list, sample.size / rep.runs), a.list, mt.list, pr.list, c.list, res.list, ep.list)
dim(param.grid)

runInd <- 1
outTabComp1 <- read.csv(paste("~/Dropbox/disjointCycles/simResults/large/large_",runInd, ".csv", sep = ""),)
missing <- c()
for(runInd in 2:nrow(param.grid)){
  if(!file.exists(paste("~/Dropbox/disjointCycles/simResults/large/large_", runInd, ".csv", sep = ""))){
    missing <- c(missing, runInd)
  } else {
    temp <- read.csv(paste("~/Dropbox/disjointCycles/simResults/large/large_", runInd, ".csv", sep = ""))

    if(any(temp$agree < 1) ){
      cat(runInd)
      cat("\n")
    }
    outTabComp1 <- rbind(outTabComp1, temp)
  }

}
missing

agTab1 <- outTabComp1 %>% group_by(p, n, mt, alpha, cycleSize, pr, dist) %>%
  summarize(order = mean(order),
            edges = mean((est_edges -p)/ (p* (p-1))),
            num = length(timeOrd),
            totalTime = mean(timeOrd) + mean(timeEdges))

table(agTab1$num)


agTabMax1 <- agTab1 %>%  group_by(p, dist, n) %>%
  summarize(maxInd = which.max(edges),
            order = max(order),
            edges = max(edges),
            totalTime = totalTime[maxInd])



setEPS()
postscript("~/Dropbox/Apps/Overleaf/Cyclic linear non-Gaussian causal models/figures/gamma3.eps", width = 8, height = 4)
par(mfrow = c(1,2))
plot(1:4, rep(-2, 4), ylim = c(0, .8),
     ylab = "Correct Order", xlab = "", xaxt = "n")
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray92")
abline(h = seq(0, 1, by = .2), lty = 3, col = "white", lwd = 2)
lines(unlist(agTabMax1[which(agTabMax1$p == 30 & agTabMax1$dist == "gamma"), 5]), col = "darkred", type = "b")
lines(unlist(agTabMax1[which(agTabMax1$p == 45 & agTabMax1$dist == "gamma"), 5]), col = "darkgreen", type = "b")
lines(unlist(agTabMax1[which(agTabMax1$p == 60 & agTabMax1$dist == "gamma"), 5]), col = "darkblue", type = "b")
axis(side = 1, labels = n.list/1000, at = 1:4)


plot(1:4, rep(-2, 4), ylim = c(.7, .9),
     ylab = "Correct Edges", xlab = "", xaxt = "n")
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray92")
abline(h = seq(0, 1, by = .05), lty = 3, col = "white", lwd = 2)
axis(side = 1, labels = n.list/1000, at = 1:4)
lines(unlist(agTabMax1[which(agTabMax1$p == 30 & agTabMax1$dist == "gamma"), 6]), col = "darkred", type = "b")
lines(unlist(agTabMax1[which(agTabMax1$p == 45 & agTabMax1$dist == "gamma"), 6]), col = "darkgreen", type = "b")
lines(unlist(agTabMax1[which(agTabMax1$p == 60 & agTabMax1$dist == "gamma"), 6]), col = "darkblue", type = "b")

par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), mar=c(0, 0, 0, 0), new=TRUE)

plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n', main = "")
legend("bottom", legend = c("p =  ",  "30", "45", "60"),
       pch = c(NA, 1, 1, 1), col = c( "white", "darkred", "darkgreen", "darkblue"),
       ncol = 4)
mtext("Gamma Errors", line = -2)
dev.off()



setEPS()
postscript("~/Dropbox/Apps/Overleaf/Cyclic linear non-Gaussian causal models/figures/norm3.eps", width = 8, height = 4)
par(mfrow = c(1,2))
plot(1:4, rep(-2, 4), ylim = c(0, .8),
     ylab = "Correct Order", xlab = "", xaxt = "n")
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray92")
abline(h = seq(0, 1, by = .2), lty = 3, col = "white", lwd = 2)
axis(side = 1, labels = n.list/1000, at = 1:4)
lines(unlist(agTabMax1[which(agTabMax1$p == 30 & agTabMax1$dist == "mixedNorm"), 5]), col = "darkred", type = "b")
lines(unlist(agTabMax1[which(agTabMax1$p == 45 & agTabMax1$dist == "mixedNorm"), 5]), col = "darkgreen", type = "b")
lines(unlist(agTabMax1[which(agTabMax1$p == 60 & agTabMax1$dist == "mixedNorm"), 5]), col = "darkblue", type = "b")


plot(1:4, rep(-2, 4), ylim = c(.7, .9),
     ylab = "Correct Edges", xlab = "", xaxt = "n")
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray92")
abline(h = seq(0, 1, by = .05), lty = 3, col = "white", lwd = 2)
axis(side = 1, labels = n.list/1000, at = 1:4)

lines(unlist(agTabMax1[which(agTabMax1$p == 30 & agTabMax1$dist == "mixedNorm"), 6]), col = "darkred", type = "b")
lines(unlist(agTabMax1[which(agTabMax1$p == 45 & agTabMax1$dist == "mixedNorm"), 6]), col = "darkgreen", type = "b")
lines(unlist(agTabMax1[which(agTabMax1$p == 60 & agTabMax1$dist == "mixedNorm"), 6]), col = "darkblue", type = "b")



par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), mar=c(0, 0, 0, 0), new=TRUE)

plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n', main = "")
legend("bottom", legend = c("p =  ",  "30", "45", "60"),
       pch = c(NA, 1, 1, 1), col = c( "white", "darkred", "darkgreen", "darkblue"),
       ncol = 4)
mtext("Mixture of Normal Errors", line = -2)
dev.off()




setEPS()
postscript("~/Dropbox/Apps/Overleaf/Cyclic linear non-Gaussian causal models/figures/combined.eps",
           width = 5, height = 5)
par(mfrow = c(2,2))
par(mar = c(1,.5, .5, .2), oma = c(5, 4, 2, 1))

plot(1:4, rep(-2, 4), ylim = c(0, .8),
     ylab = "Correct Order", xlab = "", xaxt = "n")
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray92")
abline(h = seq(0, 1, by = .2), lty = 3, col = "white", lwd = 2)
lines(unlist(agTabMax1[which(agTabMax1$p == 30 & agTabMax1$dist == "gamma"), 5]), col = "darkred", type = "b")
lines(unlist(agTabMax1[which(agTabMax1$p == 45 & agTabMax1$dist == "gamma"), 5]), col = "darkgreen", type = "b")
lines(unlist(agTabMax1[which(agTabMax1$p == 60 & agTabMax1$dist == "gamma"), 5]), col = "darkblue", type = "b")
mtext("Gamma Errors")
mtext("% Correct Order", side = 2, line = 3)


plot(1:4, rep(-2, 4), ylim = c(0, .8),
     ylab = "", xlab = "", xaxt = "n", yaxt = "n")
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray92")
abline(h = seq(0, 1, by = .2), lty = 3, col = "white", lwd = 2)
lines(unlist(agTabMax1[which(agTabMax1$p == 30 & agTabMax1$dist == "mixedNorm"), 5]), col = "darkred", type = "b")
lines(unlist(agTabMax1[which(agTabMax1$p == 45 & agTabMax1$dist == "mixedNorm"), 5]), col = "darkgreen", type = "b")
lines(unlist(agTabMax1[which(agTabMax1$p == 60 & agTabMax1$dist == "mixedNorm"), 5]), col = "darkblue", type = "b")
mtext("Mixture of Normal Errors")


plot(1:4, rep(-2, 4), ylim = c(.7, .9),
     ylab = "Correct Edges", xlab = "", xaxt = "n")
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray92")
abline(h = seq(0, 1, by = .05), lty = 3, col = "white", lwd = 2)
axis(side = 1, labels = n.list/1000, at = 1:4)
lines(unlist(agTabMax1[which(agTabMax1$p == 30 & agTabMax1$dist == "gamma"), 6]), col = "darkred", type = "b")
lines(unlist(agTabMax1[which(agTabMax1$p == 45 & agTabMax1$dist == "gamma"), 6]), col = "darkgreen", type = "b")
lines(unlist(agTabMax1[which(agTabMax1$p == 60 & agTabMax1$dist == "gamma"), 6]), col = "darkblue", type = "b")
mtext("% Correct Edges", side = 2, line = 3)

plot(1:4, rep(-2, 4), ylim = c(.7, .9),
     ylab = "", xlab = "", xaxt = "n", yaxt = "n")
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray92")
abline(h = seq(0, 1, by = .05), lty = 3, col = "white", lwd = 2)
axis(side = 1, labels = n.list/1000, at = 1:4)

lines(unlist(agTabMax1[which(agTabMax1$p == 30 & agTabMax1$dist == "mixedNorm"), 6]), col = "darkred", type = "b")
lines(unlist(agTabMax1[which(agTabMax1$p == 45 & agTabMax1$dist == "mixedNorm"), 6]), col = "darkgreen", type = "b")
lines(unlist(agTabMax1[which(agTabMax1$p == 60 & agTabMax1$dist == "mixedNorm"), 6]), col = "darkblue", type = "b")



par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), mar=c(0, 0, 0, 0), new=TRUE)

plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n', main = "")
legend("bottom", legend = c("p =  ",  "30", "45", "60"),
       pch = c(NA, 1, 1, 1), col = c( "white", "darkred", "darkgreen", "darkblue"),
       ncol = 4)
dev.off()






#
#
# setEPS()
# postscript("~/Dropbox/Apps/Overleaf/Cyclic linear non-Gaussian causal models/figures/large_sims_time.eps", width = 8, height = 4)
# plot(1:4, rep(.5, 4), ylim = c(10, 5000),
#      ylab = "Time (sec)", xlab = "n (in 000s)", xaxt = "n", log = "y")
#
#
# abline(h = seq(0, 1, by = .05), lty = 3, col = "white", lwd = 2)
# axis(side = 1, labels = n.list/1000, at = 1:4)
# mtext(side = 2, "Time (sec)", line = 3)
# lines(unlist(agTabMax1[which(agTabMax1$p == 30 & agTabMax1$dist == "gamma"), 7]), col = "darkred", type = "b")
# lines(unlist(agTabMax1[which(agTabMax1$p == 45 & agTabMax1$dist == "gamma"), 7]), col = "darkgreen", type = "b")
# lines(unlist(agTabMax1[which(agTabMax1$p == 60 & agTabMax1$dist == "gamma"), 7]), col = "darkblue", type = "b")
#
# lines(unlist(agTabMax1[which(agTabMax1$p == 30 & agTabMax1$dist == "mixedNorm"), 7]), col = "darkred", type = "b", lty = 2)
# lines(unlist(agTabMax1[which(agTabMax1$p == 45 & agTabMax1$dist == "mixedNorm"), 7]), col = "darkgreen", type = "b", lty = 2)
# lines(unlist(agTabMax1[which(agTabMax1$p == 60 & agTabMax1$dist == "mixedNorm"), 7]), col = "darkblue", type = "b")
# legend("bottomleft", legend = c("p =  ",  "30", "45", "60"),
#        pch = c(NA, 1, 1, 1), col = c( "white", "darkred", "darkgreen", "darkblue"),
#        ncol = 4)
# legend("bottomright", legend = c("gamma", "mix norm"),
#        col = c("black"), lty = c(1, 2),
#        ncol = 2)
# dev.off()



library("disjointCycles")
library("tidyverse")
sample.size <- 100
rep.runs <- 2
p.list <- c(30, 45, 60)
n.list <- c(10, 25, 50, 100) * 1000
a.list <- c(2e-1, 5e-2, 1e-2, 1e-3)
mt.list <- c("BH", "holm")
pr.list <- c("chisq")
c.list <- c(3)
res.list <- c(T)
ep.list <- c(1/2)
param.grid <- expand.grid(rep(n.list, sample.size / rep.runs), a.list, mt.list, pr.list, c.list, res.list, ep.list)
dim(param.grid)

runInd <- 1
outTabComp1 <- read.csv(paste("~/Dropbox/disjointCycles/simResults/large/large5_",runInd, ".csv", sep = ""),)
missing <- c()
for(runInd in 2:nrow(param.grid)){
  if(!file.exists(paste("~/Dropbox/disjointCycles/simResults/large/large5_", runInd, ".csv", sep = ""))){
    missing <- c(missing, runInd)
  } else {
    temp <- read.csv(paste("~/Dropbox/disjointCycles/simResults/large/large5_", runInd, ".csv", sep = ""))

    if(any(temp$agree < 1) ){
      cat(runInd)
      cat("\n")
    }
    outTabComp1 <- rbind(outTabComp1, temp)
  }

}
missing

agTab1 <- outTabComp1 %>% group_by(p, n, mt, alpha, cycleSize, pr, dist) %>%
  summarize(order = mean(order),
            edges = mean((est_edges -p)/ (p* (p-1))),
            num = length(timeOrd),
            totalTime = mean(timeOrd) + mean(timeEdges))

table(agTab1$num)


agTabMax1 <- agTab1 %>%  group_by(p, dist, n) %>%
  summarize(maxInd = which.max(edges),
            order = max(order),
            edges = max(edges),
            totalTime = totalTime[maxInd])


agTabMax1 %>% print(n = 100)

setEPS()
postscript("~/Dropbox/Apps/Overleaf/Cyclic linear non-Gaussian causal models/figures/gamma5.eps", width = 8, height = 4)
par(mfrow = c(1,2))
plot(1:4, rep(-2, 4), ylim = c(0, .8),
     ylab = "Correct Order", xlab = "", xaxt = "n")
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray92")
abline(h = seq(0, 1, by = .2), lty = 3, col = "white", lwd = 2)
lines(unlist(agTabMax1[which(agTabMax1$p == 30 & agTabMax1$dist == "gamma"), 5]), col = "darkred", type = "b")
lines(unlist(agTabMax1[which(agTabMax1$p == 45 & agTabMax1$dist == "gamma"), 5]), col = "darkgreen", type = "b")
lines(unlist(agTabMax1[which(agTabMax1$p == 60 & agTabMax1$dist == "gamma"), 5]), col = "darkblue", type = "b")
axis(side = 1, labels = n.list/1000, at = 1:4)


plot(1:4, rep(-2, 4), ylim = c(.7, .9),
     ylab = "Correct Edges", xlab = "", xaxt = "n")
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray92")
abline(h = seq(0, 1, by = .05), lty = 3, col = "white", lwd = 2)
axis(side = 1, labels = n.list/1000, at = 1:4)
lines(unlist(agTabMax1[which(agTabMax1$p == 30 & agTabMax1$dist == "gamma"), 6]), col = "darkred", type = "b")
lines(unlist(agTabMax1[which(agTabMax1$p == 45 & agTabMax1$dist == "gamma"), 6]), col = "darkgreen", type = "b")
lines(unlist(agTabMax1[which(agTabMax1$p == 60 & agTabMax1$dist == "gamma"), 6]), col = "darkblue", type = "b")

par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), mar=c(0, 0, 0, 0), new=TRUE)

plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n', main = "")
legend("bottom", legend = c("p =  ",  "30", "45", "60"),
       pch = c(NA, 1, 1, 1), col = c( "white", "darkred", "darkgreen", "darkblue"),
       ncol = 4)
mtext("Gamma Errors", line = -2)
mtext("5-cycles", line = -3)
dev.off()



setEPS()
postscript("~/Dropbox/Apps/Overleaf/Cyclic linear non-Gaussian causal models/figures/norm5.eps", width = 8, height = 4)
par(mfrow = c(1,2))
plot(1:4, rep(-2, 4), ylim = c(0, .8),
     ylab = "Correct Order", xlab = "", xaxt = "n")
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray92")
abline(h = seq(0, 1, by = .2), lty = 3, col = "white", lwd = 2)
axis(side = 1, labels = n.list/1000, at = 1:4)
lines(unlist(agTabMax1[which(agTabMax1$p == 30 & agTabMax1$dist == "mixedNorm"), 5]), col = "darkred", type = "b")
lines(unlist(agTabMax1[which(agTabMax1$p == 45 & agTabMax1$dist == "mixedNorm"), 5]), col = "darkgreen", type = "b")
lines(unlist(agTabMax1[which(agTabMax1$p == 60 & agTabMax1$dist == "mixedNorm"), 5]), col = "darkblue", type = "b")


plot(1:4, rep(-2, 4), ylim = c(.7, .9),
     ylab = "Correct Edges", xlab = "", xaxt = "n")
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray92")
abline(h = seq(0, 1, by = .05), lty = 3, col = "white", lwd = 2)
axis(side = 1, labels = n.list/1000, at = 1:4)

lines(unlist(agTabMax1[which(agTabMax1$p == 30 & agTabMax1$dist == "mixedNorm"), 6]), col = "darkred", type = "b")
lines(unlist(agTabMax1[which(agTabMax1$p == 45 & agTabMax1$dist == "mixedNorm"), 6]), col = "darkgreen", type = "b")
lines(unlist(agTabMax1[which(agTabMax1$p == 60 & agTabMax1$dist == "mixedNorm"), 6]), col = "darkblue", type = "b")



par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), mar=c(0, 0, 0, 0), new=TRUE)

plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n', main = "")
legend("bottom", legend = c("p =  ",  "30", "45", "60"),
       pch = c(NA, 1, 1, 1), col = c( "white", "darkred", "darkgreen", "darkblue"),
       ncol = 4)
mtext("Mixture of Normal Errors", line = -2)
mtext("5-cycles", line = -3)
dev.off()


#
# par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), mar=c(0, 0, 0, 0), new=TRUE)
#
# plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
# legend("bottom", legend = c("p =  ",  "30", "45", "60"),
#        pch = c(NA, 1, 1, 1), col = c( "white", "darkred", "darkgreen", "darkblue"),
#        ncol = 4)
# dev.off()
#
#
# setEPS()
# postscript("~/Dropbox/Apps/Overleaf/Cyclic linear non-Gaussian causal models/figures/large_sims_time.eps", width = 8, height = 4)
# plot(1:4, rep(.5, 4), ylim = c(10, 5000),
#      ylab = "Time (sec)", xlab = "n (in 000s)", xaxt = "n", log = "y")
#
#
# abline(h = seq(0, 1, by = .05), lty = 3, col = "white", lwd = 2)
# axis(side = 1, labels = n.list/1000, at = 1:4)
# mtext(side = 2, "Time (sec)", line = 3)
# lines(unlist(agTabMax1[which(agTabMax1$p == 30 & agTabMax1$dist == "gamma"), 7]), col = "darkred", type = "b")
# lines(unlist(agTabMax1[which(agTabMax1$p == 45 & agTabMax1$dist == "gamma"), 7]), col = "darkgreen", type = "b")
# lines(unlist(agTabMax1[which(agTabMax1$p == 60 & agTabMax1$dist == "gamma"), 7]), col = "darkblue", type = "b")
#
# lines(unlist(agTabMax1[which(agTabMax1$p == 30 & agTabMax1$dist == "mixedNorm"), 7]), col = "darkred", type = "b", lty = 2)
# lines(unlist(agTabMax1[which(agTabMax1$p == 45 & agTabMax1$dist == "mixedNorm"), 7]), col = "darkgreen", type = "b", lty = 2)
# lines(unlist(agTabMax1[which(agTabMax1$p == 60 & agTabMax1$dist == "mixedNorm"), 7]), col = "darkblue", type = "b")
# legend("bottomleft", legend = c("p =  ",  "30", "45", "60"),
#        pch = c(NA, 1, 1, 1), col = c( "white", "darkred", "darkgreen", "darkblue"),
#        ncol = 4)
# legend("bottomright", legend = c("gamma", "mix norm"),
#        col = c("black"), lty = c(1, 2),
#        ncol = 2)
# dev.off()

