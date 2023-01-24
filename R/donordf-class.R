#!/usr/bin/env R

# Author: Sean Maden
#
# Defines the donor.data.frame class, for tables containing donor-level signals 
# by marker and type. 
#
#

setClass(
  "donor.data.frame",
  representation(type_summary = "character", annotation = "character"),
  contains = "data.frame"
)
