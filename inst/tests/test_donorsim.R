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
donordf = NULL
method = "nbinom"
lambda.pos = 20
lambda.neg = 2
lambda.sdoff.pos = 5
lambda.sdoff.neg = 2
gamma.pos = 20
gamma.neg = 2
P = c(0.25, 0.75)
donor.adj.method = "combat"
plot.biasadj = TRUE
plot.pca = TRUE
cname.donorsummary = "donor.combn.all.mean"
gindexv = c(1, 2)
ndonor = 10
seed.num = 0
verbose = TRUE

set.seed(seed.num); lr <- list()
if(verbose){message("Making pseudobulk sample from types matrix...")}
df <- rand_donor_marker_table(ndonor = 1, gindexv = gindexv,
                              lambda.sdoff.pos = 0, lambda.sdoff.neg = 0)
ktotal <- length(P)
Z <- matrix(df[,"donor1"], ncol = ktotal)
Ypb <- ypb_fromtypes(Z = Z, P = P)

if(is(donordf, "NULL")){
  if(verbose){message("Getting randomized donor marker data...")}  
  donordf <- rand_donor_marker_table(ndonor = ndonor, 
                                     method = method,
                                     gindexv = gindexv,
                                     lambda.sdoff.pos = lambda.sdoff.pos,
                                     lambda.sdoff.neg = lambda.sdoff.neg,
                                     gamma.pos = gamma.pos,
                                     gamma.neg = gamma.neg,
                                     seed.num = seed.num)
} else{
  if(verbose){message("Checking if provided donordf passes checks...")}
  if(!check_donordf(donordf)){
    stop("Error, provided donordf is invalid. ",
         "Check that it contains all required columns.")}
}

if(verbose){message("Getting type predictions...")}
type.indexv <- seq(ktotal)
donor.unadj <- donordf[,cname.donorsummary] # get donor summary datas
lbias <- biasexpt(df = donordf, Ypb = Ypb, P = P, donor.unadj = donor.unadj,
               donor.adj.method = donor.adj.method,
               plot.biasadj = plot.biasadj,
               verbose = verbose)

# get return object
lr[["dfres"]] <- lbias$dfi # get results df
lr[["donordf"]] <- donordf
lr[["Ypb"]] <- Ypb
lr[["adj.method"]] <- method
# make new plots
if(plot.pca){lr[["lpca"]] <- pcaplots_donor(dt = donordf)}
if(plot.biasadj){lr[["ggpt.biasadj"]] <- lbias$ggpt.biasadj}
return(lr)
