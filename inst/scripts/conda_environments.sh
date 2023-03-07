#!/usr/bin/env sh

# Author: Sean Maden
#
# Install conda dependencies for intronomer.
#

#-----------------
# make basic r env
#-----------------
conda create --name r_4 r=4
activate r_4
conda install -c conda-forge r-devtools

#-------------------
# Get SCDC conda env
#-------------------
conda create --name scdc --clone r_411
conda activate scdc
conda install -c conda-forge r-devtools
R
devtools::install_github("renozao/xbioc")
devtools::install_github("meichendong/SCDC")

conda env export > scdc.yml

#--------------------
# Get MuSiC conda env
#--------------------
conda create --name music --clone r_4
conda activate music
# conda install -c conda-forge toast
# conda install -c bioconda bioconductor-biobase
# conda install -c bioconda bioconductor-singlecellexperiment
R
install.packages("BiocManager")
BiocManager::install("TOAST")
BiocManager::install("Biobase")
BiocManager::install("SingleCellExperiment")
devtools::install_github("xuranw/MuSiC")

