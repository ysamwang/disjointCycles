sample.size <- 500
rep.runs <- 25
n.list <- c(10, 25, 50) * 1000
c.list <- c(2, 3, 4)
a.list <- c(.05, .01, .001, 1e-4, 1e-5)
param.grid <- expand.grid(rep(n.list, sample.size / rep.runs), c.list, a.list)
# 450

runInd <- 1
outTabComp <- read.csv(paste("~/Dropbox/disjointCycles/simResults/cycleChain/cycleChain_", runInd, ".csv", sep = ""))
colnames(outTabComp) <- c("cycLen", "alpha", "n", paste(rep(c("gamma", "mixNorm"), each = 5), rep(c("order", "oracle_sen", "oracle_spc", "est_sen", "est_spc"), times= 2), sep = "_"))
missing <- c()
for(runInd in 2:nrow(param.grid)){
  if(!file.exists(paste("~/Dropbox/disjointCycles/simResults/cycleChain/cycleChain_", runInd, ".csv", sep = ""))){
    missing <- c(missing, runInd)
  } else {
    temp <- read.csv(paste("~/Dropbox/disjointCycles/simResults/cycleChain/cycleChain_", runInd, ".csv", sep = ""))
    colnames(temp) <- c("cycLen", "alpha", "n", paste(rep(c("gamma", "mixNorm"), each = 5), rep(c("order", "oracle_sen", "oracle_spc", "est_sen", "est_spc"), times= 2), sep = "_"))
    outTabComp <- rbind(outTabComp, temp)
  }

}
missing

agTab <- aggregate(cbind(gamma_order, mixNorm_order, gamma_oracle_sen, mixNorm_oracle_sen, gamma_oracle_spc , mixNorm_oracle_spc,
                         gamma_graph_oracle = gamma_oracle_sen == 1 & gamma_oracle_spc == 1, mixNorm_graph_oracle =mixNorm_oracle_sen == 1 & mixNorm_oracle_spc == 1,
                         gamma_est_spc, mixNorm_est_spc, gamma_est_sen, mixNorm_est_sen,
                         gamma_graph_est = gamma_est_sen == 1 & gamma_est_spc == 1, mixNorm_graph_est = mixNorm_est_sen == 1 & mixNorm_est_spc == 1
                         )~ n + cycLen + alpha, data = outTabComp, FUN = mean)


par(mfrow = c(length(n.list), length(c.list)))
for(n in n.list){
  for(cyc in c.list){

    plot(agTab[which(agTab$n == n & agTab$cycLen == cyc), "gamma_order" ],xaxt = "n",
         ylim = c(0, 1), type = "b", ylab = "", main = paste("n = ", n, "; cycle = ", cyc, sep = " "))
    lines(agTab[which(agTab$n == n & agTab$cycLen == cyc), "mixNorm_order" ], type = "b", lty = 2, pch = 2)
    axis(side = 1, at = 1:5, labels = rev(a.list))

  }
}


par(mfrow = c(length(n.list), length(c.list)))
for(n in n.list){
  for(cyc in c.list){

    plot(agTab[which(agTab$n == n & agTab$cycLen == cyc), "gamma_graph_oracle" ],xaxt = "n",
         ylim = c(0, 1), type = "b", ylab = "", main = paste("n = ", n, "; cycle = ", cyc, sep = " "))
    lines(agTab[which(agTab$n == n & agTab$cycLen == cyc), "mixNorm_graph_oracle" ], type = "b", lty = 2, pch = 2)
    axis(side = 1, at = 1:5, labels = rev(a.list))

  }
}


par(mfrow = c(length(n.list), length(c.list)))
for(n in n.list){
  for(cyc in c.list){

    plot(1- agTab[which(agTab$n == n & agTab$cycLen == cyc), "gamma_oracle_spc" ],
         agTab[which(agTab$n == n & agTab$cycLen == cyc), "gamma_oracle_sen" ],
         ylim = c(.8, 1), xlim = c(0, .2), type = "b", ylab = "", main = paste("n = ", n, "; cycle = ", cyc, sep = " "))
    lines(1- agTab[which(agTab$n == n & agTab$cycLen == cyc), "mixNorm_oracle_spc" ],
          agTab[which(agTab$n == n & agTab$cycLen == cyc), "mixNorm_graph_sen" ],
          type = "b", lty = 2, pch = 2)
  }
}







agTab1 <- aggregate(cbind(gamma_order, mixNorm_order, gamma_oracle_sen, mixNorm_oracle_sen, gamma_oracle_spc, mixNorm_oracle_spc,
                         gamma_est_spc, mixNorm_est_spc)~ n + cycLen, data = agTab, FUN = max)









#############################
sample.size <- 400
rep.runs <- 10
n.list <- c(10, 25, 50) * 1000
a.list <- c(.01, .001)
param.grid <- expand.grid(rep(n.list, sample.size / rep.runs), a.list)


runInd <- 1
outTabComp <- read.csv(paste("~/Dropbox/disjointCycles/simResults/ex32/ex32_", runInd, ".csv", sep = ""))
for(runInd in 2:nrow(param.grid)){
  temp <- read.csv(paste("~/Dropbox/disjointCycles/simResults/ex32/ex32_", runInd, ".csv", sep = ""))
  outTabComp <- rbind(outTabComp, temp)

}

aggregate(cbind(gamma, mixNorm)~ n + alpha, data = outTabComp, FUN = mean)



#############################################
sample.size <- 200
rep.runs <- 10
n.list <- c(10, 25, 50) * 1000
c.list <- c(2, 3, 4)
a.list <- c(.05, .01, .001, 1e-4, 1e-5)
param.grid <- expand.grid(rep(n.list, sample.size / rep.runs), c.list, a.list)
# 450

runInd <- 1
outTabComp <- read.csv(paste("~/Dropbox/disjointCycles/simResults/singleCycle/singleCycle_",runInd, ".csv", sep = ""),)
colnames(outTabComp) <- c("cycLen", "alpha", "n", paste(rep(c("gamma", "mixNorm"), each = 5), rep(c("order", "oracle_sen", "oracle_spc", "est_sen", "est_spc"), times= 2), sep = "_"))
missing <- c()
for(runInd in 2:nrow(param.grid)){
  if(!file.exists(paste("~/Dropbox/disjointCycles/simResults/singleCycle/singleCycle_", runInd, ".csv", sep = ""))){
    missing <- c(missing, runInd)
  } else {
    temp <- read.csv(paste("~/Dropbox/disjointCycles/simResults/singleCycle/singleCycle_", runInd, ".csv", sep = ""))
    colnames(temp) <- c("cycLen", "alpha", "n", paste(rep(c("gamma", "mixNorm"), each = 5), rep(c("order", "oracle_sen", "oracle_spc", "est_sen", "est_spc"), times= 2), sep = "_"))
    outTabComp <- rbind(outTabComp, temp)
  }

}
missing

agTab <- aggregate(cbind(gamma_order, mixNorm_order, gamma_oracle_sen, mixNorm_oracle_sen, gamma_oracle_spc , mixNorm_oracle_spc,
                         gamma_graph_oracle = gamma_oracle_sen == 1 & gamma_oracle_spc == 1, mixNorm_graph_oracle =mixNorm_oracle_sen == 1 & mixNorm_oracle_spc == 1,
                         gamma_est_spc, mixNorm_est_spc, gamma_est_sen, mixNorm_est_sen,
                         gamma_graph_est = gamma_est_sen == 1 & gamma_est_spc == 1,
                         mixNorm_graph_est = mixNorm_est_sen == 1 & mixNorm_est_spc == 1)~
                     n + cycLen + alpha, data = outTabComp, FUN = mean)




agTab2 <- aggregate(cbind(gamma_order, mixNorm_order, gamma_oracle_sen, mixNorm_oracle_sen, gamma_oracle_spc , mixNorm_oracle_spc,
                         gamma_graph_oracle, mixNorm_graph_oracle,
                         gamma_est_spc, mixNorm_est_spc, gamma_est_sen, mixNorm_est_sen,
                         gamma_graph_est, mixNorm_graph_est
)~ n + cycLen, data = agTab, FUN = max)


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



