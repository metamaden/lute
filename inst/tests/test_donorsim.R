#!/usr/bin/env R

# Author: Sean Maden
#
# Test donor simulations.
#

require(lute)

df <- rand_donor_marker_table(ndonor = 10, gindexv = c(1, 2), 
                              sd.offset.pos = 5,sd.offset.neg = 5)

# perform combat adj
donordf <- rand_donor_marker_table()
donorv <- donordf$donor.combn.all.mean
donorv.adj <- donoradj(donorv, donordf, method = "combat")

# simulate donor signals data
donordf <- rand_donor_marker_table()
# get ypb
P <- c(0.75, 0.25)
ktotal <- length(P)
Z <- matrix(df[,"donor1"], ncol = ktotal)
Ypb <- ypb_fromtypes(Z = Z, P = P)
# get bias expt results
li <- biasexpt(df = donordf, Ypb = Ypb)