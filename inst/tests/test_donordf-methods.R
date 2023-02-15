#!/usr/bin/env R

# tests for donor.data.frame methods

require(lute)

#--------------------------
# test donordf_from_mexpr()
#--------------------------

# args

mexpr
typev.rd
groupv.cd
verbose

mexpr = assays(scef)$logcounts
typev.rd = NULL
groupv.cd = scef$donor
typev.cd = scef$k2
verbose = T

# function

utypev.cd <- unique(typev.cd)
ugroupv.cd <- unique(groupv.cd)

if(is(typev.rd, "NULL")){ 
  # get marker.type var
  mgk <- do.call(cbind, lapply(utypev.cd, function(ti){
    matrix(rowMeans(mexpr[,typev.cd==ti]), ncol = 1)
  }))
  typev.rd <- unlist(apply(mgk, 1, function(ri){
    utypev.cd[ri==max(ri)]}))
}

df <- do.call(cbind, lapply(ugroupv.cd, function(di){
  filt.donor <- groupv.cd==di
  do.call(rbind, lapply(utypev.cd, function(ti){
    filt.type <- filt.donor & typev.cd==ti
    mdat <- rowMeans(mexpr[,filt.type], na.rm = T)
    matrix(mdat, ncol = 1)
  }))
}))

colnames(df) <- as.character(ugroupv.cd)
donor.mean <- rowMeans(df, na.rm = T)
donor.median <- rowMedians(df, na.rm = T)

df <- as.data.frame(df)
df$donor.mean <- donor.mean
df$donor.median <- donor.median
df$type <- rep(utypev.cd, each = nrow(mexpr))
df$marker <- rep(rownames(mexpr), length(utypev.cd))
df$marker.type <- rep(typev.rd, length(utypev.cd))
df <- as(df, "donor.data.frame")
if(!check_donordf(df)){stop("Error, couldn't make new valid donor.data.frame")}


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

function (ndonor = 2, gindexv = c(1, 2), method = "nbinom", lambda.pos = 20, 
          lambda.neg = 2, lambda.sdoff.pos = 0, lambda.sdoff.neg = 0, 
          gamma.pos = 10, gamma.neg = 10, seed.num = 0, verbose = FALSE, 
          ...) 
{
  set.seed(seed.num)
  nmarkers <- length(gindexv)
  ktotal <- length(unique(gindexv))
  offposv <- rnorm(n = ndonor, mean = 0, sd = lambda.sdoff.pos)
  offnegv <- rnorm(n = ndonor, mean = 0, sd = lambda.sdoff.neg)
  meanv.pos <- offposv + lambda.pos
  meanv.neg <- offnegv + lambda.neg
  meanv.pos[meanv.pos < 0] <- -1 * meanv.pos[meanv.pos < 0]
  meanv.neg[meanv.neg < 0] <- -1 * meanv.neg[meanv.neg < 0]
  md <- do.call(cbind, lapply(seq(ndonor), function(ii) {
    unlist(random_lgv(gindexv, num.iter = 1, lambda.pos = meanv.pos[ii], 
                      lambda.neg = meanv.neg[ii], gamma.size.pos = gamma.pos, 
                      gamma.size.neg = gamma.neg, method = method, seed.num = ii))
  }))
  md <- as.data.frame(md)
  colnames(md) <- paste0("donor", seq(ndonor))
  if (ndonor > 1) {
    if (verbose) {
      message("Getting donor summary columns...")
    }
    which.cnv.donor <- which(grepl("donor", colnames(md)))
    md$donor.mean <- apply(md[, which.cnv.donor], 1, mean)
    md$donor.median <- apply(md[, which.cnv.donor], 1, median)
  }
  type.vector <- paste0("type", rep(seq(ktotal), each = nmarkers))
  type.levels <- unique(type.vector)
  type.levels <- marker.levels[order(type.levels)]
  type.factor <- factor(type.vector, levels = type.levels)
  md$type <- type.factor
  
  marker.vector <- paste0("marker", rep(seq(nmarkers), times = ktotal))
  marker.levels <- unique(marker.vector)
  marker.levels <- marker.levels[order(marker.levels)]
  marker.factor <- factor(marker.vector, levels = marker.levels)
  md$marker <- marker.factor
  
  md$marker.type <- paste0("type", gindexv)
  md <- as(md, "donor.data.frame")
  if (check_donordf(md)) {
    return(md)
  }
  else {
    stop("Error, couldn't make new donor.data.frame object.")
  }
  return(NULL)
}