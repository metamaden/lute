#!/usr/bin/env sh

# Author: Sean Maden
#
# Install conda dependencies for intronomer.
#

# make new env
# conda create --name scdc r=4.2
# conda activate scdc

#-----------------
# make basic r env
#-----------------
conda create --name r_411 r=4.1.1
conda activate music
conda install -c conda-forge r-devtools

#-------------------
# Get SCDC conda env
#-------------------

# make new env
# module load conda_R/4.2
# conda create --name scdc r=4.2.2
#conda create --name scdc
#conda activate scdc
#conda install -c conda-forge r-base=4.1.2
#conda activate scdc

conda create --name scdc
conda activate scdc
conda install -c conda-forge r-devtools
R
devtools::install_github("renozao/xbioc")
devtools::install_github("meichendong/SCDC")

conda env export > scdc.yml

#-------------------
# Get MuSiC conda env
#-------------------
conda create --name music
conda activate music
conda install -c conda-forge r-devtools
R
devtools::install_github("renozao/xbioc")
devtools::install_github("meichendong/SCDC")

module load conda_R/4.2
R
devtools::install_github("xuranw/MuSiC")

