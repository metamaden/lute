#!/usr/bin/env R

# Author: Sean Maden
#
# Test functions for NextFlow workflow.
#

# test filter_value_type()
set.seed(0)
sce <- random_sce(zero.include = T, zero.fract = 0.3)
scef <- filter_value_type(sce, filter.term = "zerocount", verbose = T,
                          max.gene.value.freq = 0.15)
metadata(scef)$filter.zerocount.by.type$df.type


