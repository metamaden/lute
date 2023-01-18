#!/usr/bin/env R

# Author: Sean Maden
#
# Test donor simulations.
#

require(lute)

df <- rand_donor_marker_table(ndonor = 10, gindexv = c(1, 2), 
                              sd.offset.pos = 5,sd.offset.neg = 5)