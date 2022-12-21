#!/usr/bin/env R

# Author: Sean Maden
#
# Quantify donor bias; show impact on deconvolution outcomes with simulations.
#
#

libv <- c("lute")
sapply(libv, library, character.only = T)

# get simulated marker data
dt <- rand_donor_marker_table(ndonor = 2, lambda.pos = 10)