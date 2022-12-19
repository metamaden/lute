library(lute)

ndonor <- 2
mean.offset <- 10
gindexv <- c(1,2)
lambda.pos = 25
lambda.neg = 2

lgv <- random_lgv(gindexv)
meanv.pos <- rpois(n = ndonor, lambda = mean.offset + lambda.pos)
meanv.neg <- rpois(n = ndonor, lambda = mean.offset + lambda.neg)

llgv <- lapply(seq(ndonor), function(ii){
  random_lgv(gindexv, 
             lambda.pos = meanv.pos[ii], 
             lambda.neg = meanv.neg[ii])
})
