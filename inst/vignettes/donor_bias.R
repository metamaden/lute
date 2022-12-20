library(lute)

ndonor <- 2
mean.offset.pos <- 10
mean.offset.neg <- 1
lambda.pos = 25
lambda.neg = 2

gindexv <- c(1,2)
lambda.pos = 25
lambda.neg = 2
ktotal <- 2
num.iter <- 1

# draw random offsets from normal dist
offposv <- rnorm(n = ndonor, lambda = mean.offset.pos)
offnegv <- rnorm(n = ndonor, lambda = mean.offset.neg)
# get value vectors
meanv.pos <- offposv + lambda.pos
meanv.neg <- offnegv + lambda.pos
# convert negative means
meanv.pos[meanv.pos < 0] <- -1*meanv.pos
meanv.neg[meanv.neg < 0] <- -1*meanv.neg


# get matrix of markers (rows) by donors (cols)
md <- do.call(cbind, lapply(seq(ndonor), function(ii){
  unlist(random_lgv(gindexv, num.iter = num.iter,
                    lambda.pos = meanv.pos[ii],
                    lambda.neg = meanv.neg[ii]))
}))
md <- as.data.frame(md)
colnames(md) <- paste0("donor", seq(ndonor))
md$marker <- paste0("marker", rep(seq(ktotal), each = ndonor))
md$type <- paste0("type", rep(seq(ktotal), ndonor))

# get summary stats
# by types
lapply(unique(md$marker), function(mi){
  mdi <- md[md$marker == mi,]
  apply(mdi[,c(1:ndonor)], 1, var)
})

# by markers
lapply(unique(md$type), function(ti){
  mdi <- md[md$type == ti,]
  apply(mdi[,c(1:ndonor)], 1, var)
})

lgv.total <- apply(md, 2, mean)
lgvi <- list(length(gindexv))













ndonor = 5
set.seed(seed.num)
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
  unlist(random_lgv(gindexv, num.iter = num.iter,
                    lambda.pos = meanv.pos[ii],
                    lambda.neg = meanv.neg[ii],
                    ...))
}))
md <- as.data.frame(md)
colnames(md) <- paste0("donor", seq(ndonor))
md$marker <- paste0("marker", rep(seq(ktotal), each = ndonor))
md$type <- paste0("type", rep(seq(ktotal), ndonor))
