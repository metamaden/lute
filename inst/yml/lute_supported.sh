#!/usr/bin/env bash

# lute_supported.sh
#
# Author: Sean Maden
#
# Install dependencies from shell/bash command line calls.
#
#
#

#
# chained calls
#
url=(
	https://github.com/metamaden/musicParam
	https://github.com/metamaden/music2Param
	https://github.com/LieberInstitute/DeconvoBuddies
	https://github.com/metamaden/meanratiosParam
	https://github.com/metamaden/deconrnaseqParam
	https://github.com/metamaden/epicParam
	https://github.com/xuranw/MuSiC
	git clone https://github.com/Jiaxin-Fan/MuSiC2
	https://github.com/GfellerLab/EPIC
	https://github.com/shenorrLabTRDF/bseqsc
	https://github.com/shenorrLabTRDF/csSAM
	https://github.com/renozao/xbioc
)
for url in "${urlarr[@]}" do
	git clone $url
done
pkgArr=(
	csSAM
	xbioc
	bseqsc
	musicParam
	music2Param
	meanratiosParam
	deconrnaseqParam
	EPIC
	epicParam
)	
for pkg in "${pkgArr[@]}" do
	R CMD build $pkg
	R CMD INSTALL $pkg_*
done
pkgArr=(
	bseqsc
	DeconvoBuddies
	MuSiC
	MuSiC2
)
for pkg in "${pkgArr[@]}" do
	R CMD build $pkg --no-build-vignettes
	R CMD INSTALL $pkg_*
done

#
# successive calls
#

git clone https://github.com/metamaden/musicParam
git clone https://github.com/metamaden/music2Param
git clone https://github.com/LieberInstitute/DeconvoBuddies
git clone https://github.com/metamaden/bisqueParam
git clone https://github.com/metamaden/meanratiosParam
git clone https://github.com/metamaden/deconrnaseqParam
git clone https://github.com/metamaden/epicParam
git clone https://github.com/xuranw/MuSiC
git clone https://github.com/Jiaxin-Fan/MuSiC2
git clone https://github.com/GfellerLab/EPIC
git clone https://github.com/shenorrLabTRDF/bseqsc
git clone https://github.com/shenorrLabTRDF/csSAM
git clone https://github.com/renozao/xbioc

R CMD build csSAM
R CMD INSTALL csSAM_*

R CMD build xbioc
R CMD INSTALL xbioc_*

R CMD build bseqsc
R CMD build bseqsc

R CMD build MuSiC --no-build-vignettes
R CMD INSTALL MuSiC_*

R CMD build MuSiC2 --no-build-vignettes
R CMD INSTALL MuSiC2_*

R CMD build musicParam
R CMD INSTALL musicParam_*

R CMD build music2Param
R CMD INSTALL music2Param_*

R CMD build meanratiosParam
R CMD INSTALL meanratiosParam_*

R CMD build deconrnaseqParam
R CMD INSTALL deconrnaseqParam_*

R CMD build EPIC
R CMD INSTALL EPIC_*

R CMD build epicParam
R CMD INSTALL epicParam_*

R CMD build bseqsc --no-build-vignettes
R CMD INSTALL bseqsc_*

R CMD build DeconvoBuddies --no-build-vignettes
R CMD INSTALL DeconvoBuddies_*