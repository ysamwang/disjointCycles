n <- 500
p <- 10

### Compare testing jointly vs
sim.size <- 1000
rec <- matrix(0, sim.size, 8)
for(i in 1:sim.size){
  dat <- cdcs::rDAG(p, n, dist = "gamma", parent_prob = .5, lowEdge = .5, highEdge = .9)
  Y <- scale(dat$Y, center = F)

  moments <- .calcSandT(Y)



  sMat <- moments$sMat; tMat <- moments$tMat

  u <- 1

  rec[i, 1] <- .test2joint(u, moments$sMat, moments$tMat, Y, method = "jointDelta")$p
  rec[i, 2] <- .test2joint(u, moments$sMat, moments$tMat, Y, method = "indDelta")$p

  u <- 2

  rec[i, 3] <- .test2joint(u, moments$sMat, moments$tMat, Y, method = "jointDelta")$p
  rec[i, 4] <- .test2joint(u, moments$sMat, moments$tMat, Y, method = "indDelta")$p

  u <- floor(p/2)

  rec[i, 5] <- .test2joint(u, moments$sMat, moments$tMat, Y, method = "jointDelta")$p
  rec[i, 6] <- .test2joint(u, moments$sMat, moments$tMat, Y, method = "indDelta")$p

  u <- p-1

  rec[i, 7] <- .test2joint(u, moments$sMat, moments$tMat, Y, method = "jointDelta")$p
  rec[i, 8] <- .test2joint(u, moments$sMat, moments$tMat, Y, method = "indDelta")$p

  if(i%%100==0){cat(i);cat("\n")}
}

alpha <- .05

final <- colMeans(rec < alpha)
names(final) <- paste(rep(c("size", "two", "half", "end"), each = 2), rep(c("joint", "Ind"), times = 4), sep = "_")
final


