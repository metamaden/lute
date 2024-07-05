#!/usr/bin/env R

# lute_supported.R
#
# Author: Sean Maden
#
# Install dependencies for the lute_bioconductor Docker image.
#
#
install.packages("nnls")
install.packages("BisqueRNA")
install.packages("MCMCpack")
install.packages("pkgmaker")
install.packages("openxlsx")
install.packages("rngtools")
install.packages("NMF")
install.packages("Formula")
install.packages("checkmate")
install.packages("rafalib")
install.packages("RefManageR")
install.packages("e1071")
install.packages("reshape2")
install.packages("openxlsx")
BiocManager::install("TOAST")
BiocManager::install("DeconRNASeq")
BiocManager::install("BiocStyle")
BiocManager::install("spatialLIBD")