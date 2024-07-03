#!/usr/bin/env R

#
# Run from active R session to complete lute.yml conda environment setup.
#

install.packages("BiocManager")
BiocManager::install("lute")