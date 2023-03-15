# lute

[<img style="float: right;" src = "inst/png/lute_hexsticker_basic1.png" height="200"/>](https://github.com/metamaden/lute)

R package for deconvolution simulation and optimization.

## Installation

## Generics, methods, and classes for deconvolution

The `lute` R package introduces several new generics for deconvolution experiments. It defines the new 
`deconvolution` generic, whose methods provide access to multiple state-of-the-science algorithm and maps 
standard inputs to synonymous function-specific arguments. Since there are numerous routine tasks shared across many deconvolution experiments, deconvolution algorithms are supported using a new class hierarchy that is useful
for understanding and comparing current methods, and can serve as a starting point for new experiments and
new method development.

The inheritance for the new class hierarchy for deconvolution methods in `lute` is shown in the following 
flowchart:

[<img style="float: right;" src = "inst/png/deconvolution_classes_flowchart.jpeg" height="500"/>](https://github.com/metamaden/lute)

## Deconvolution methods supported

The following methods are currently supported in `lute`:

* MuSiC : The `musicParam` class supports the `MuSiC::music.basic()` implementation of the MuSiC deconvolution algorithm (see also: ).

* MuSiC2 : The `music2Param` class supports both the `MuSiC::music2_prop()` and `MuSiC2::music2_prop()` implementations of the MuSiC2 deconvolution algorithm (see also: ).

* EPIC : The `epicParam` class supports the `EPIC::EPIC()` implementation of the EPIC deconvolution algorithm (see also: ).

* DeconRNASeq : The `deconrnaseqParam` class supports the `DeconRNASeq::DeconRNASeq()` implementation of the DeconRNASeq deconvolution algorithm (see also: ).

* SCDC : The `scdcParam` class supports the `SCDC::SCDC_prop()` implementation of the SCDC deconvolution algorithm (see also: ).

* NNLS : The `nnlsParam` class supports the `nnls::nnls` implementation of the NNLS deconvolution algorithm (see also: ).

* Bisque : The `bisqueParam` class supports the `BisqueRNA::ReferenceBasedDecomposition` implementation of the Bisque deconvolution algorithm (see also: ).

