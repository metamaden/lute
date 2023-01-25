#!/usr/bin/env R

# tests for donor.data.frame methods

require(lute)

#--------------------------
# test donordf_from_mexpr()
#--------------------------
mexpr
typev.rd
groupv.cd
verbose

mexpr = assays(scef)$logcounts
typev.rd = NULL
groupv.cd = scef$donor
typev.cd = scef$k2
verbose = T



utypev.cd <- unique(typev.cd)
if(is(typev.rd, "NULL")){ 
  # get marker.type var
  mgk <- do.call(cbind, lapply(utypev.cd, function(ti){
    matrix(rowMeans(mexpr[,typev.cd==ti]), ncol = 1)
  }))
  typev.rd <- unlist(apply(mgk, 1, function(ri){
    utypev.cd[ri==max(ri)]}))
}

df <- do.call(cbind, lapply(udonorv, function(di){
  filt.donor <- which(donorv==di)
  do.call(rbind, lapply(utypev.cd, function(ti){
    filt.type <- filt.donor & which(utypev.cd==ti)
    matrix(rowMeans(mexpr[,filt.type]), ncol = 1)
  }))
}))
colnames(df) <- as.character(udonorv)
donor.mean <- rowMeans(df)
donor.median <- rowMedians(df)
df <- as.data.frame(df)
df$donor.mean <- donor.mean
df$donor.median <- donor.median
df$type <- rep(utypev.cd, each = nrow(mexpr))
df$marker <- rep(rownames(mexpr), length(utypev.cd))
df$marker.type <- typev.rd
df <- as(df, "donor.data.frame")
if(!check_donordf(df)){stop("Error, couldn't make new valid donor.data.frame")}
return(df)



#-----------------------
# test random_donordf()
#-----------------------
ndonor = 2
gindexv = c(1, 2)
method = "nbinom"
lambda.pos = 20
lambda.neg = 2
lambda.sdoff.pos = 0
lambda.sdoff.neg = 0
gamma.pos = 10
gamma.neg = 10
seed.num = 0
verbose = FALSE

set.seed(seed.num)
nmarkers <- length(gindexv); ktotal <- length(unique(gindexv))

# get random hyperparameter offsets from normal dist
# get random offsets for means
offposv <- rnorm(n = ndonor, mean = 0, sd = lambda.sdoff.pos)
offnegv <- rnorm(n = ndonor, mean = 0, sd = lambda.sdoff.neg)

# get new hyperparameter values
# get new means
meanv.pos <- offposv + lambda.pos
meanv.neg <- offnegv + lambda.neg

# convert negative values
# convert means
meanv.pos[meanv.pos < 0] <- -1*meanv.pos[meanv.pos < 0]
meanv.neg[meanv.neg < 0] <- -1*meanv.neg[meanv.neg < 0]

# get matrix of markers (rows) by donors (cols)
md <- do.call(cbind, lapply(seq(ndonor), function(ii){
  unlist(random_lgv(gindexv, num.iter = 1, 
                    lambda.pos = meanv.pos[ii],
                    lambda.neg = meanv.neg[ii], 
                    gamma.size.pos = gamma.pos,
                    gamma.size.neg = gamma.neg, 
                    method = method, seed.num = ii))
}))
md <- as.data.frame(md); colnames(md) <- paste0("donor", seq(ndonor))
if(ndonor > 1){
  if(verbose){message("Getting donor summary columns...")}
  which.cnv.donor <- which(grepl("donor", colnames(md)))
  md$donor.combn.all.mean <- apply(md[,which.cnv.donor], 1, mean)
  md$donor.combn.all.median <- apply(md[,which.cnv.donor], 1, median)
}
md$type <- paste0("type", rep(seq(ktotal), each = nmarkers))
md$marker <- paste0("marker", rep(seq(nmarkers), times = ktotal))
md$marker.type <- paste0("type", gindexv)
# final check
md <- as(md, "donor.data.frame")