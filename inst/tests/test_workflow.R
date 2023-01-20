#!/usr/bin/env R

# Author: Sean Maden
#
# Test functions for NextFlow workflow.
#

sce <- random_sce(na.include = T, na.fract = 0.4, 
                  zero.include = T, zero.fract = 0.4)
scef <- filter_value_cells(sce, verbose = T)