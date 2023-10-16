#!/usr/bin/env sh

### Author: Sean Maden
###
### Install conda dependencies for intronomer.
###

###-----------------
### make basic r env
###-----------------
conda create --name r_4 r=4
activate r_4
conda install -c conda-forge r-devtools

###-------------------
### Get NNLS conda env
###-------------------
conda create --name nnls --clone r_422
conda activate nnls
R
install.packages("nnls", repos = "https://cloud.r-project.org")

conda env export > nnls.yml

###-------------------
### Get SCDC conda env
###-------------------
conda create --name scdc --clone r_411
conda activate scdc
conda install -c conda-forge r-devtools
R
devtools::install_github("renozao/xbioc")
devtools::install_github("meichendong/SCDC")

conda env export > scdc.yml

###--------------------
### Get MuSiC conda env
###--------------------
conda create --name music --clone r_4
conda activate music
## conda install -c conda-forge toast
## conda install -c bioconda bioconductor-biobase
## conda install -c bioconda bioconductor-singlecellexperiment
R
install.packages("BiocManager")
BiocManager::install("TOAST")
BiocManager::install("Biobase")
BiocManager::install("SingleCellExperiment")
devtools::install_github("xuranw/MuSiC")

conda env export > music.yml

###---------------------
### Get MuSiC2 conda env
###---------------------
conda create --name music2 --clone music
conda activate music2
R
devtools::install_github("renozao/xbioc")
devtools::install_github("Jiaxin-Fan/MuSiC2")

conda env export > music2.yml

###-------------------
### Get SCDC conda env
###-------------------
conda create --name scdc --clone r_411
conda activate scdc
R
remotes::install_github("renozao/xbioc")
devtools::install_github("meichendong/SCDC")

conda env export > scdc.yml

###--------------------
### Get EPIC conda env
###--------------------
conda create --name epic --clone r_4
conda activate epic
R
devtools::install_github("GfellerLab/EPIC")

conda env export > epic.yml

###--------------------------
### Get DeconRNASeq conda env
###--------------------------
conda create --name deconrnaseq --clone r_4
conda activate deconrnaseq
R
install.packages("BiocManager")
BiocManager::install("DeconRNASeq")

conda env export > deconrnaseq.yml

###---------------------
### Get Bisque conda env
###---------------------
conda create --name bisque --clone r_4
conda activate bisque
R
devtools::install_github("cozygene/bisque")

conda env export > bisque.yml