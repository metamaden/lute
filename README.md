# lute

Authors: Sean Maden, Stephanie Hicks

[<img style="float: right;" src = "inst/png/lute_hexsticker_basic1.png" height="200"/>](https://github.com/metamaden/lute)

The `lute` R package, a framework and utilities for deconvolution experiments.

## Installation

### From GitHub

Install `lute` from GitHub by running the following in an R session:

```
devtools::install("metamaden/lute")
```

### From Bioconductor (pending)

Install `lute` from Bioconductor with the following:

```
BiocManager::install("lute")
```

## Generics, methods, and classes for deconvolution

The `lute` R package introduces several new generics for deconvolution experiments. It defines the new 
`deconvolution` generic, whose methods provide access to multiple state-of-the-science algorithms and maps 
standard inputs to synonymous function-specific arguments. 

Since there are numerous routine tasks shared across many deconvolution experiments, deconvolution algorithms are supported using a new class hierarchy that is useful
for understanding and comparing current methods, and can serve as a starting point for new experiments and
new method development.

The inheritance hierarchy for the numerous novel classes in `lute` for various deconvolution algorithms is depicted in the following flowchart:

[<img style="float: center;" src = "inst/jpg/deconvolutionParam_hierarchy_diagram.jpeg"/>](https://github.com/metamaden/lute)

## Deconvolution methods supported

The following deconvolution algorithms are currently supported with novel methods and classes in `lute`. They are listed with links to their main repos as well as conda YML scripts to install them with dependencies.

* MuSiC : The `musicParam` class supports the `MuSiC::music.basic()` implementation of the MuSiC deconvolution algorithm ([url](https://github.com/xuranw/MuSiC); [yml](https://github.com/metamaden/lute/blob/main/inst/yml/music.yml)).

* MuSiC2 : The `music2Param` class supports both the `MuSiC::music2_prop()` and `MuSiC2::music2_prop()` implementations of the MuSiC2 deconvolution algorithm ([url](https://github.com/Jiaxin-Fan/MuSiC2); [yml](https://github.com/metamaden/lute/blob/main/inst/yml/music2.yml)).

* EPIC : The `epicParam` class supports the `EPIC::EPIC()` implementation of the EPIC deconvolution algorithm ([url](https://github.com/GfellerLab/EPIC); [yml](https://github.com/metamaden/lute/blob/main/inst/yml/epic.yml)).

* DeconRNASeq : The `deconrnaseqParam` class supports the `DeconRNASeq::DeconRNASeq()` implementation of the DeconRNASeq deconvolution algorithm ([url](https://bioconductor.org/packages/release/bioc/html/DeconRNASeq.html); [yml](https://github.com/metamaden/lute/blob/main/inst/yml/deconrnaseq.yml)).

* SCDC : The `scdcParam` class supports the `SCDC::SCDC_prop()` implementation of the SCDC deconvolution algorithm ([url](https://github.com/meichendong/SCDC); [yml](https://github.com/metamaden/lute/blob/main/inst/yml/scdc.yml)).

* NNLS : The `nnlsParam` class supports the `nnls::nnls` implementation of the NNLS deconvolution algorithm ([url](https://cran.r-project.org/web/packages/nnls/index.html); [yml](https://github.com/metamaden/lute/blob/main/inst/yml/nnls.yml)).

* Bisque : The `bisqueParam` class supports the `BisqueRNA::ReferenceBasedDecomposition` implementation of the Bisque deconvolution algorithm ([url](https://github.com/cozygene/bisque); [yml](https://github.com/metamaden/lute/blob/main/inst/yml/bisque.yml)).

## Conda environment support

Set up a new conda environment for a deconvolution method by passing the following code your Terminal:

```
conda env create -f ./lute/inst/yml/music.yml
```

## Further reading

Additional resources to learn more about deconvolution algorithms and their applications.

* [`awesome-deconvolution`](https://github.com/metamaden/awesome-deconvolution) : List of deconvolution algorithms, frameworks, and resources across various disciplines.