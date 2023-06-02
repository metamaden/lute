#!/usr/bin/env R

# Author: Sean Maden
#
# Tests using standard formats for the nnlsParam class.
#

require(lute)

lparam_from_lexample <- function(lexample){
  
  lparam <- referencebasedParam(
    
    y = lexample[["y"]], z = lexample[["z"]], s = lexample[["s"]]
    
  )
  
  return(lparam)
}

lparam <- lparam_from_lexample(lute:::.get_decon_example_data())

output.from.one.bulk.sample <- function(algo.index.from.mappings){
  mappings.data <- lute:::algo.from.index()
  
  
  
  eval(parse(text = paste0(algorithm.string, "(",z, y,")")))
}

#-------------------
# embedded algorithm
#-------------------
# test algorithm contained by the new param class

require(nnls)
bulk.samples.index.vector <- seq(ncol(y))
standard.output.first.index <- nnls::nnls(A = z, b = lparam[["y"]][,1])




class(standard.output.first.index) == "nnls"
paste0(names(standard.output.first.index), collapse = ";") == "x;deviance;residuals;fitted;mode;passive;bound;nsetp"



result <- lapply(
  
  bulk.samples.index.vector, 
  function(index){
    
    nnls::nnls(A = z, b = y[,index])
    
  }
  
  
)

names(result) <- colnames(y)
predictions <- lapply(result, 
                      function(iter){
                        iter$x
                        }
                      )

lr <- .parse_deconvolution_predictions_results(predictions, 
                                               colnames(z), 
                                               colnames(y))

#------------
# default use
#------------
# run new param class with any applicable generics




param <- nnlsParam(s = lexample[["s"]], y = lexample[["y"]], z = lexample[["z"]])

# return only predicted proportions
result <- deconvolution(param)

# return detailed results
param[["return.info"]] <- T

result <- deconvolution(param)

