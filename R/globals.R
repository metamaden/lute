#!/usr/bin/env R

###
### Defines global variables used in the package. For instance, names variables 
### used in dplyr functions filter(), arrange(), etc.
###

utils::globalVariables(c("cellType.target", "rank_ratio", 
                         "overlapping.samples", "abs.summary", "marker.table"))