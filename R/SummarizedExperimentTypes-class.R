#!/usr/bin/env R

# Author: Sean Maden
#
# Define set class, validity tests, maker and conversion functions, etc.
#


setClass(
  "SummarizedExperimentSets",
  representation(preprocessMethod = "character", annotation = "character"),
  contains = "SummarizedExperiment"
)
