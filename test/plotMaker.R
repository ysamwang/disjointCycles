library("tidyverse")
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
# 450

runInd <- 1
outTabComp <- read.csv(paste("~/Dropbox/disjointCycles/simResults/cycleChain/cycleChain_", runInd, ".csv", sep = ""))

missing <- c()
for(runInd in 2:nrow(param.grid)){
  if(!file.exists(paste("~/Dropbox/disjointCycles/simResults/cycleChain/cycleChain_", runInd, ".csv", sep = ""))){
    missing <- c(missing, runInd)
  } else {
    temp <- read.csv(paste("~/Dropbox/disjointCycles/simResults/cycleChain/cycleChain_", runInd, ".csv", sep = ""))
    outTabComp <- rbind(outTabComp, temp)
  }

}
missing


agTab <- outTabComp %>% group_by(n, cycleSize, alpha, mt, pr, rescaleData) %>%
  summarize(g_ord = mean(gamma_order), m_ord = mean(mixNorm_order),
            ln_ord = mean(ln_order),
            g_edge = mean(gamma_est_edges), m_edge = mean(mixNorm_est_edges),
            ln_edge = mean(ln_est_edges))


agTabMax <- agTab %>% filter(rescaleData == T) %>%group_by(cycleSize, n) %>%
  summarize(g_ord = max(g_ord), m_ord = max(m_ord), ln_ord = max(ln_ord),
            g_edge = max(g_edge), m_edge = max(m_edge),
            ln_edge = max(ln_edge)) %>% print(n = 50)

setEPS()
postscript("~/Dropbox/disjointCycles/simResults/plots/cycleChain12.eps", width = 6, height = 4)
par(mfrow = c(2,3))
par(mar = c(1, .5, 1, .5))
par(oma = c(4.5, 4, 2,0))
for(cInd in c.list){
  plot(agTabMax$g_ord[which(agTabMax$cycleSize == cInd)], type = "n",
       main = "", xaxt = "n", ylab = "", yaxt = "n",
       xlab = "", ylim = c(0, 1), pch = NA)
  mtext(paste("Cycle Size =", cInd))
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray92")
  abline(h = seq(0, 1, by = .2), lty = 3, col = "white", lwd = 1.5)
  if(cInd == 2){
    mtext("% Cor. Ordering", side = 2, line = 2)
    axis(side = 2, at = seq(0, 1, by = .2))
  }
  lines(agTabMax$g_ord[which(agTabMax$cycleSize == cInd)],
        type= "b", col = "darkred", pch = 1, lwd =2)
  lines(agTabMax$m_ord[which(agTabMax$cycleSize == cInd)],
        type= "b", col = "darkblue", pch = 2, lwd =2)
  lines(agTabMax$ln_ord[which(agTabMax$cycleSize == cInd)],
        col = "darkgreen", pch = 3, type = "b", lwd =2)
}

for(cInd in c.list){
  plot(agTabMax$g_edge[which(agTabMax$cycleSize == cInd)] - 12, type = "n",
       main = "", xaxt = "n", ylab = "", yaxt = "n",
       xlab = "", ylim = c(60, 12*11), pch = NA)
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray92")
  axis(side = 1, at = 1:3, labels = n.list/ 1000)
  abline(h = seq(0, 12 * 11, by = 12), lty = 3, col = "white", lwd = 1.5)
  if(cInd == 2){
    mtext("Num. Cor. Edges", side = 2, line = 2)
    axis(side = 2, at = seq(0, 12 * 11, by = 12))
  }

  lines(agTabMax$g_edge[which(agTabMax$cycleSize == cInd)] - 12,
        type= "b", col = "darkred", pch = 1, lwd =2)
  lines(agTabMax$m_edge[which(agTabMax$cycleSize == cInd)] - 12,
        type= "b", col = "darkblue", pch = 2, lwd =2)
  lines(agTabMax$ln_edge[which(agTabMax$cycleSize == cInd)] - 12,
        col = "darkgreen", pch = 3, type = "b", lwd =2)
}
par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), mar=c(0, 0, 0, 0), new=TRUE)
plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
legend("bottom", legend = c("Mixture", "Gamma", "Lognormal"),
       pch = c(2,1,3), col = c("darkred", "darkblue", "darkgreen"),
       ncol = 3)
dev.off()


###################################################
library("tidyverse")
sample.size <- 400
rep.runs <- 25
n.list <- c(10, 25, 50) * 1000
a.list <- c(2.5e-1, 1e-1, 1e-2, 1e-3, 1e-4)
mt.list <- c("BH", "holm")
pr.list <- c("naive")
c.list <- c(2, 4)
res.list <- c(T)
param.grid <- expand.grid(rep(n.list, sample.size / rep.runs), a.list, mt.list, pr.list, c.list, res.list)
dim(param.grid)


runInd <- 1
outTabComp8 <- read.csv(paste("~/Dropbox/disjointCycles/simResults/cycleChain8/cycleChain8_", runInd, ".csv", sep = ""))

missing <- c()
for(runInd in 2:nrow(param.grid)){
  if(!file.exists(paste("~/Dropbox/disjointCycles/simResults/cycleChain8/cycleChain8_", runInd, ".csv", sep = ""))){
    missing <- c(missing, runInd)
  } else {
    temp <- read.csv(paste("~/Dropbox/disjointCycles/simResults/cycleChain8/cycleChain8_", runInd, ".csv", sep = ""))
    outTabComp8 <- rbind(outTabComp8, temp)
  }
}
missing


agTab8 <- outTabComp8 %>% group_by(n, cycleSize, alpha, mt, pr, rescaleData) %>%
  summarize(g_ord = mean(gamma_order), m_ord = mean(mixNorm_order),
            ln_ord = mean(ln_order),
            g_edge = mean(gamma_est_edges), m_edge = mean(mixNorm_est_edges),
            ln_edge = mean(ln_est_edges))

agTabMax8 <- agTab8 %>% filter(rescaleData == T) %>%group_by(cycleSize, n) %>%
  summarize(g_ord = max(g_ord), m_ord = max(m_ord),
            ln_ord = max(ln_ord),
            g_edge = max(g_edge), m_edge = max(m_edge),
            ln_edge = max(ln_edge)) %>% print(n = 50)

p <- 8
setEPS()
postscript("~/Dropbox/disjointCycles/simResults/plots/cycleChain8.eps", width = 6, height = 4)
par(mfrow = c(2,2))
par(mar = c(1, .5, 1, .5))
par(oma = c(4.5, 4, 2,0))
for(cInd in c.list){
  plot(agTabMax8$g_ord[which(agTabMax8$cycleSize == cInd)], type = "n",
       main = "", xaxt = "n", ylab = "", yaxt = "n",
       xlab = "", ylim = c(0, 1), pch = NA)
  mtext(paste("Cycle Size =", cInd))
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray92")
  abline(h = seq(0, 1, by = .2), lty = 3, col = "white", lwd = 1.5)
  if(cInd == 2){
    mtext("% Cor. Ordering", side = 2, line = 2)
    axis(side = 2, at = seq(0, 1, by = .2))
  }
  lines(agTabMax8$g_ord[which(agTabMax8$cycleSize == cInd)],
        type= "b", col = "darkred", pch = 1, lwd =2)
  lines(agTabMax8$m_ord[which(agTabMax8$cycleSize == cInd)],
        type= "b", col = "darkblue", pch = 2, lwd =2)
  lines(agTabMax8$ln_ord[which(agTabMax8$cycleSize == cInd)],
        col = "darkgreen", pch = 3, type = "b", lwd =2)
}

for(cInd in c.list){
  plot(agTabMax8$g_edge[which(agTabMax8$cycleSize == cInd)] - 12, type = "n",
       main = "", xaxt = "n", ylab = "", yaxt = "n",
       xlab = "", ylim = c(35, p*(p-1)), pch = NA)
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray92")
  axis(side = 1, at = 1:3, labels = n.list/ 1000)
  abline(h = seq(0, p*(p-1), by = p), lty = 3, col = "white", lwd = 1.5)
  if(cInd == 2){
    mtext("Num. Cor. Edges", side = 2, line = 2)
    axis(side = 2, at = seq(0, p*(p-1), by = p))
  }

  lines(agTabMax8$g_edge[which(agTabMax8$cycleSize == cInd)] - p,
        type= "b", col = "darkred", pch = 1, lwd =2)
  lines(agTabMax8$m_edge[which(agTabMax8$cycleSize == cInd)] - p,
        type= "b", col = "darkblue", pch = 2, lwd =2)
  lines(agTabMax8$ln_edge[which(agTabMax8$cycleSize == cInd)] - p,
        col = "darkgreen", pch = 3, type = "b", lwd =2)
}
par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), mar=c(0, 0, 0, 0), new=TRUE)
plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
legend("bottom", legend = c("Mixture", "Gamma", "Lognormal"),
       pch = c(2,1,3), col = c("darkred", "darkblue", "darkgreen"),
       ncol = 3)
dev.off()


#############################
library("tidyverse")
sample.size <- 400
rep.runs <- 20
n.list <- c(10, 25, 50) * 1000
a.list <- c(2.5e-1, 1e-1, 1e-2, 1e-3, 1e-4)
mt.list <- c("BH", "holm")
pr.list <- c("naive")
res.list <- c(T)
param.grid <- expand.grid(rep(n.list, sample.size / rep.runs), a.list, mt.list, pr.list, res.list)
dim(param.grid)



runInd <- 1
outTabCompEx32 <- read.csv(paste("~/Dropbox/disjointCycles/simResults/ex32/ex32_", runInd, ".csv", sep = ""))
for(runInd in 2:nrow(param.grid)){
  temp <- read.csv(paste("~/Dropbox/disjointCycles/simResults/ex32/ex32_", runInd, ".csv", sep = ""))
  outTabCompEx32 <- rbind(outTabCompEx32, temp)
}


agTabEx32 <- outTabCompEx32 %>% group_by(n, alpha, mt, pr, rescaleData) %>%
  summarize(g_ord = mean(gamma_order), m_ord = mean(mixNorm_order),
            ln_ord = mean(ln_order),
            g_edge = mean(gamma_est_edges), m_edge = mean(mixNorm_est_edges),
            ln_edge = mean(ln_est_edges))

agTabMaxEx32 <- agTabEx32 %>% filter(rescaleData == T) %>%group_by(n) %>%
  summarize(g_ord = max(g_ord), m_ord = max(m_ord), ln_ord = max(ln_ord),
            g_edge = max(g_edge), m_edge = max(m_edge), ln_edge = max(ln_edge)) %>% print(n = 50)

agTabEx32 %>% filter(rescaleData == T, n == 50000)
agTabEx32 %>% filter(rescaleData == T) %>%group_by(n) %>%
  summarize(g_ord = which.max(g_ord), m_ord = which.max(m_ord),
            g_edge = which.max(g_edge), m_edge = which.max(m_edge)) %>% print(n = 50)

setEPS()
postscript("~/Dropbox/disjointCycles/simResults/plots/ex32.eps", width = 6, height = 4)
par(mfrow = c(1,2))
par(mar = c(2, 4, 2, .5))
par(oma = c(2,0,0,0))

plot(agTabMaxEx32$g_ord, type = "n",
       main = "", xaxt = "n", ylab = "", yaxt = "n",
       xlab = "", ylim = c(0, 1), pch = NA)
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray92")
  abline(h = seq(0, 1, by = .2), lty = 3, col = "white", lwd = 1.5)
    mtext("% Cor. Ordering", side = 2, line = 2)
    axis(side = 2, at = seq(0, 1, by = .2))
    axis(side = 1, at =1:3, labels = n.list/1000)
  lines(agTabMaxEx32$g_ord,
        type= "b", col = "darkred", pch = 1, lwd =2)
  lines(agTabMaxEx32$m_ord,
        type= "b", col = "darkblue", pch = 2, lwd =2)
  lines(agTabMaxEx32$ln_ord,
        col = "darkgreen", pch = 3, type = "b", lwd =2)

plot(agTabMaxEx32$g_edge, type = "n",
       main = "", xaxt = "n", ylab = "", yaxt = "n",
       xlab = "", ylim = c(44, 56), pch = NA)
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray92")
  abline(h = seq(40, p * (p-1), by = p/2), lty = 3, col = "white", lwd = 1.5)
  mtext("Num. Cor. Edges", side = 2, line = 2)
  axis(side = 2, at = seq(40, 56, by = 4))
  lines(agTabMaxEx32$g_edge - 8,
        type= "b", col = "darkred", pch = 1, lwd =2)
  lines(agTabMaxEx32$m_edge -8,
        type= "b", col = "darkblue", pch = 2, lwd =2)
  lines(agTabMaxEx32$ln_edge -8,
        col = "darkgreen", pch = 3, type = "b", lwd =2)
  axis(side = 1, at =1:3, labels = n.list/1000)


par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), mar=c(0, 0, 0, 0), new=TRUE)
plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
mtext(side = 3, "Example 3.2", line = -1)
legend("bottom", legend = c("Mixture", "Gamma", "Lognormal"),
       pch = c(2,1,3), col = c("darkred", "darkblue", "darkgreen"),
       ncol = 3)
dev.off()

#############################################
sample.size <- 500
rep.runs <- 100
n.list <- c(5, 10, 25, 50) * 1000
a.list <- c(1e-1, 1e-2, 1e-3, 1e-4)
mt.list <- c("BH", "holm")
pr.list <- c("infFunc", "naive")
c.list <- c(2, 3, 4)
param.grid <- expand.grid(rep(n.list, sample.size / rep.runs), a.list, mt.list, pr.list, c.list)
# 960
# 450

runInd <- 1
outTabComp <- read.csv(paste("~/Dropbox/disjointCycles/simResults/cycleChain/cycleChain_",runInd, ".csv", sep = ""),)
missing <- c()
for(runInd in 2:nrow(param.grid)){
  if(!file.exists(paste("~/Dropbox/disjointCycles/simResults/cycleChain/cycleChain_", runInd, ".csv", sep = ""))){
    missing <- c(missing, runInd)
  } else {
    temp <- read.csv(paste("~/Dropbox/disjointCycles/simResults/cycleChain/cycleChain_", runInd, ".csv", sep = ""))
    outTabComp <- rbind(outTabComp, temp)
  }

}
missing

agTab <- aggregate(cbind(gamma_order, mixNorm_order, gamma_oracle_sen, mixNorm_oracle_sen, gamma_oracle_spc , mixNorm_oracle_spc,
                         gamma_graph_oracle = gamma_oracle_sen == 1 & gamma_oracle_spc == 1, mixNorm_graph_oracle =mixNorm_oracle_sen == 1 & mixNorm_oracle_spc == 1,
                         gamma_est_sen, mixNorm_est_sen,gamma_est_spc, mixNorm_est_spc,
                         gamma_graph_est = gamma_est_sen == 1 & gamma_est_spc == 1, mixNorm_graph_est = mixNorm_est_sen == 1 & mixNorm_est_spc == 1)~
                     n + cycleSize + alpha + mt + pr, data = outTabComp, FUN = mean)




agTabMax <-  aggregate(cbind(gamma_order, mixNorm_order, gamma_oracle_sen, mixNorm_oracle_sen, gamma_oracle_spc , mixNorm_oracle_spc,
                             gamma_graph_oracle = gamma_oracle_sen == 1 & gamma_oracle_spc == 1, mixNorm_graph_oracle =mixNorm_oracle_sen == 1 & mixNorm_oracle_spc == 1,
                             gamma_est_sen, mixNorm_est_sen,gamma_est_spc, mixNorm_est_spc,
                             gamma_graph_est = gamma_est_sen == 1 & gamma_est_spc == 1, mixNorm_graph_est = mixNorm_est_sen == 1 & mixNorm_est_spc == 1)~
                         n + cycleSize + mt + pr, data = agTab, FUN = max)


########
sample.size <- 1000
rep.runs <- 10
n.list <- c(10, 25, 50) * 1000
c.list <- c(4)
alpha <- .05
param.grid <- expand.grid(rep(n.list, sample.size / rep.runs), c.list)


runInd <- 1
outTabComp <- read.csv(paste("~/Dropbox/disjointCycles/simResults/pvals/pvals_",runInd, ".csv", sep = ""),)

missing <- c()
for(runInd in 2:nrow(param.grid)){
  if(!file.exists(paste("~/Dropbox/disjointCycles/simResults/pvals/pvals_", runInd, ".csv", sep = ""))){
    missing <- c(missing, runInd)
  } else {
    temp <- read.csv(paste("~/Dropbox/disjointCycles/simResults/pvals/pvals_", runInd, ".csv", sep = ""))
    outTabComp <- rbind(outTabComp, temp)
  }

}
missing


tabNull <- aggregate(cbind(X1, X2, X4, X5, X6, X8, X10, X13, X14, X17, X18, X19) < .05 ~ n + type,
          data = outTabComp, FUN = mean)


tabAlt <- aggregate(cbind(X3, X7 ,X9, X11, X12, X15, X16, X20, X21, X22, X23, X24) < .05 ~ n + type,
                 data = outTabComp, FUN = mean)

#################
p <- 9
n.list <- rep(c(5, 10 , 25, 50) * 1000, each = 100)
rep.runs <- 10
cycleSize <- 3

runInd <- 1
outTabComp <- read.csv(paste("~/Dropbox/disjointCycles/simResults/pruneRoots/pruneRoots_",runInd, ".csv", sep = " "))

for(runInd in 2:length(n.list)){
  temp <- read.csv(paste("~/Dropbox/disjointCycles/simResults/pruneRoots/pruneRoots_",runInd, ".csv", sep = " "))
  outTabComp <- rbind(outTabComp, temp)
}

aggregate(cbind(X1, X2, X3, X4, X5, X6) < .1 ~n, data = outTabComp, FUN = mean)



hist(outTabComp[which(outTabComp$n == 10000), "V1"])
hist(outTabComp[which(outTabComp$n == 25000), "V1"])
hist(outTabComp[which(outTabComp$n == 50000), "V1"])

hist(outTabComp[which(outTabComp$n == 10000), "V2"])
hist(outTabComp[which(outTabComp$n == 25000), "V2"])
hist(outTabComp[which(outTabComp$n == 50000), "V2"])
