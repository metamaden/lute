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

[<img style="float: right;" src = "inst/png/deconvolution_classes_flowchart.jpg" height="500"/>](https://github.com/metamaden/lute)