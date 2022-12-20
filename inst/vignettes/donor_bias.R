library(lute)

ndonor = 5
gindexv = c(1, 2)
mean.offset.pos = 10
mean.offset.neg = 2
seed.num = 0
lambda.pos = 20
lambda.neg = 1
ktotal = 2

set.seed(seed.num)
nmarkers <- length(gindexv)
# draw random offsets from normal dist
offposv <- rnorm(n = ndonor, mean = mean.offset.pos)
offnegv <- rnorm(n = ndonor, mean = mean.offset.neg)
# get value vectors
meanv.pos <- offposv + lambda.pos
meanv.neg <- offnegv + lambda.pos
# convert negative means
meanv.pos[meanv.pos < 0] <- -1*meanv.pos
meanv.neg[meanv.neg < 0] <- -1*meanv.neg
# get matrix of markers (rows) by donors (cols)
md <- do.call(cbind, lapply(seq(ndonor), function(ii){
  unlist(random_lgv(gindexv, num.iter = 1,
                    lambda.pos = meanv.pos[ii],
                    lambda.neg = meanv.neg[ii]))
}))
md <- as.data.frame(md)
colnames(md) <- paste0("donor", seq(ndonor))
md$type <- paste0("type", rep(seq(ktotal), each = nmarkers))
md$marker <- paste0("marker", rep(seq(nmarkers), times = ktotal))
md$marker.type <- paste0("type", gindexv)