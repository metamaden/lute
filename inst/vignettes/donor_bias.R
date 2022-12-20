library(lute)

ndonor <- 2
ktotal <- 2
num.iter <- 1
mean.offset <- 10
gindexv <- c(1,2)
lambda.pos = 25
lambda.neg = 2

# randomize dist means
meanv.pos <- rpois(n = ndonor, lambda = mean.offset + lambda.pos)
meanv.neg <- rpois(n = ndonor, lambda = mean.offset + lambda.neg)
# get matrix of markers (rows) by donors (cols)
md <- do.call(cbind, lapply(seq(ndonor), function(ii){
  unlist(random_lgv(gindexv, num.iter = num.iter,
                    lambda.pos = meanv.pos[ii],
                    lambda.neg = meanv.neg[ii]))
}))
# get summary stats
# across types
mvar.cross <- apply(md, 1, var)
mvar.in <- apply(md, 2, var)
msd.cross <- apply(md, 1, sd)
msd.in <- apply(md, 2, sd)
# by types
# by markers


lgv.total <- apply(md, 2, mean)
lgvi <- list(length(gindexv))

