
<!-- badges: start -->

[![R build
status](https://github.com/metamaden/lute/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/metamaden/lute/actions)

<!-- badges: end -->


# lute

Authors: Sean Maden, Stephanie Hicks

[<img style="float: right;" src = "inst/png/lute_hexsticker_basic1.png" height="200"/>](https://github.com/metamaden/lute)

`lute` is a framework for deconvolution experiments.

## Installation

### From GitHub

Install `lute` from GitHub by running the following in an R session:

```
devtools::install("metamaden/lute")
```

## Framework overview

The `lute` framework package causes deconvolution with the `lute()` function (see `?lute` for details). Data may be incorporated from the `cellScaleFactor` package and passed to the `s` argument for transformation of cell reference expression data prior to deconvolution.   

[<img style="float: center;" src = "inst/jpg/lute_framework_diagram.jpeg"/>](https://github.com/metamaden/lute)

## Deconvolution methods supported

The following deconvolution algorithms are currently supported with novel methods and classes in `lute`. They are listed with links to their main repos as well as conda YML scripts to install them with dependencies.

* MuSiC : The `musicParam` class supports the `MuSiC::music.basic()` implementation of the MuSiC deconvolution algorithm ([url](https://github.com/xuranw/MuSiC); [yml](https://github.com/metamaden/lute/blob/main/inst/yml/music.yml)).

* MuSiC2 : The `music2Param` class supports both the `MuSiC::music2_prop()` and `MuSiC2::music2_prop()` implementations of the MuSiC2 deconvolution algorithm ([url](https://github.com/Jiaxin-Fan/MuSiC2); [yml](https://github.com/metamaden/lute/blob/main/inst/yml/music2.yml)).

* EPIC : The `epicParam` class supports the `EPIC::EPIC()` implementation of the EPIC deconvolution algorithm ([url](https://github.com/GfellerLab/EPIC); [yml](https://github.com/metamaden/lute/blob/main/inst/yml/epic.yml)).

* DeconRNASeq : The `deconrnaseqParam` class supports the `DeconRNASeq::DeconRNASeq()` implementation of the DeconRNASeq deconvolution algorithm ([url](https://bioconductor.org/packages/release/bioc/html/DeconRNASeq.html); [yml](https://github.com/metamaden/lute/blob/main/inst/yml/deconrnaseq.yml)).

* SCDC : The `scdcParam` class supports the `SCDC::SCDC_prop()` implementation of the SCDC deconvolution algorithm ([url](https://github.com/meichendong/SCDC); [yml](https://github.com/metamaden/lute/blob/main/inst/yml/scdc.yml)).

* NNLS : The `nnlsParam` class supports the `nnls::nnls` implementation of the NNLS deconvolution algorithm ([url](https://cran.r-project.org/web/packages/nnls/index.html); [yml](https://github.com/metamaden/lute/blob/main/inst/yml/nnls.yml)).

* Bisque : The `bisqueParam` class supports the `BisqueRNA::ReferenceBasedDecomposition` implementation of the Bisque deconvolution algorithm ([url](https://github.com/cozygene/bisque); [yml](https://github.com/metamaden/lute/blob/main/inst/yml/bisque.yml)).

## Conda environments

To run specific deconvolution methods and algorithms, `lute` contains conda YML scripts. For example, set up an environement to run the MuSiC method with the following:

```
conda env create -f ./lute/inst/yml/music.yml
```

## Acknowledgements

We acknowledge the following individuals for their helpful feedback and suggestions on improving this project: Louise Huuki-Myers