#!/usr/bin/env R

# Author: Sean Maden
#
# Test donor simulations.
#

require(lute)

# simulate donor type signals
df <- donordf <- rand_donor_marker_table(ndonor = 10, gindexv = c(1, 2), 
                              sd.offset.pos = 5,sd.offset.neg = 5)

# perform combat adj
donorv.adj <- donoradj(df = donordf, method = "combat")

# simulate donor signals data
# get ypb
P <- c(0.75, 0.25); ktotal <- length(P)
Z <- matrix(df[,"donor1"], ncol = ktotal)
Ypb <- ypb_fromtypes(Z = Z, P = P)
# get bias expt results
li <- biasexpt(df = donordf, Ypb = Ypb, P = P)

# get bias expt series
lb <- donor_marker_biasexpt(df = df, offsetv = c(5, 100), donor.adj.method = "combat")


# test donor_marker_biasexpt
offsetv = c(1, 10)
P = c(0.25, 0.75)
donor.adj.method = "combat"
plot.biasadj = TRUE
plot.pca = TRUE
cname.donorsummary = "donor.combn.all.mean"
gindexv = c(1, 2)
ndonor = 10
seed.num = 0
verbose = FALSE

set.seed(seed.num); lr <- list()
if(verbose){message("Making pseudobulk sample from types matrix...")}
df <- rand_donor_marker_table(ndonor = 1, gindexv = gindexv,
                              sd.offset.pos = 0, sd.offset.neg = 0)
ktotal <- length(P)
Z <- matrix(df[,"donor1"], ncol = ktotal)
Ypb <- ypb_fromtypes(Z = Z, P = P)
if(verbose){message("Getting randomized donor marker data...")}
ldonordf <- lapply(offsetv, function(offi){
  rand_donor_marker_table(ndonor = ndonor, gindexv = gindexv,
                          sd.offset.pos = offi, sd.offset.neg = offi,
                          seed.num = seed.num)
})
names(ldonordf) <- paste0("offset:", offsetv)
if(verbose){message("Getting type predictions...")}
type.indexv <- seq(ktotal)
lexpt <- lapply(seq(length(ldonordf)), function(ii){
  namei <- names(ldonordf)[ii]
  df <- ldonordf[[namei]] # get full donordf
  offsetv <- rep(gsub(".*:", "", namei), ktotal)
  donor.unadj <- df[,cname.donorsummary] # get donor summary datas
  li <- biasexpt(df = df, Ypb = Ypb, P = P, donor.unadj = donor.unadj,
                 donor.adj.method = donor.adj.method,
                 plot.biasadj = plot.biasadj,
                 verbose = verbose)
  # append offset values
  li$dfi$offset <- rep(offsetv, nrow(li$dfi)/length(offsetv))
  return(li)
})
names(lexpt) <- names(ldonordf)
# get results df
dfres <- do.call(rbind, lapply(lexpt, function(ii){ii$dfi}))
ldonorv <- lapply(lexpt, function(ii){ii[c("donor.unadj", "donor.adj")]})
names(ldonorv) <- names(ldonordf)
# get return object
lmd.adj <- list(donor.adj.method = donor.adj.method)
lmd <- list(offsetv = offsetv, P = P, donor.adj.info = lmd.adj)
lr[["dfres"]] <- dfres
lr[["ldonorv"]] <- ldonorv
lr[["ldonordf"]] <- ldonordf
lr[["Ypb"]] <- Ypb
lr[["metadata"]] <- lmd
# get plot objects
if(plot.pca){
  lpca <- lapply(ldonordf, function(dfi){
    pcaplots_donor(dt = df, title.append = NULL)
  })
  names(lpca) <- names(ldonordf)
  lr[["lpca"]] <- lpca
}
if(plot.biasadj){
  lpt <- lapply(lexpt, function(ii){ii$ggpt.biasadj})
  names(lpt) <- names(lexpt)
  lr[["ggpt.biasadj"]] <- lpt
}