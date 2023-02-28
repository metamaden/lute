#!/usr/bin/env r

# Author: Sean Maden
#
# Main functions managing mappings between standard deconvolution interface and
# specific functions.
#

#' run_deconvolution
#'
#' Perform deconvolution to predict cell type proportions in mixed samples by
#' passing standard inputs and getting standard outputs across a number of 
#' different deconvolution methods/functions.
#'
#' @param Z Signature matrix of dimensions G (marker genes) x K (types).
#' @param Y Bulk matrix of dimensions G (marker genes) x J (bulk samples).
#' @param seed.num Random seed for computational reproducibility.
#' @param method Character string of a valid deconvolution method to use (see
#' available methods with `valid_deconvolution_methods()`).
#' @param arguments List of additional valid arguments for the method.
#' @returns List of objects of type `deconvolution.results` containing 
#' predictions, metadata, and benchmarking data, organized by J samples in the
#' provided Y bulk data.
#' 
#' @details This is the main function to estimate cell type proportions using 
#' numerous different reference-based deconvolution methods. It calls several
#' other functions and provides a standard output from supported methods. 
#' Currently, one bulk/mixed signals sample is processed at a time. If multiple
#' bulk samples are provided (e.g. `ncol(Y) > 1`), then results for each
#' sample are returned in a named list of `deconvolution.results` objects.
#' 
#' This function addresses several issues found across bulk deconvolution 
#' methods. First, it provides a standard way of calling the deconvolution 
#' functions, including a single way of specifying the signature matrix (i.e. 
#' argument `Z`) and bulk signals matrix (i.e. argument `Y`). This is needed 
#' because most function arguments use non-standard references to these objects,
#' such as "a", "B", "X", "signatures", etc. 
#' 
#' Second, predictions are returned in standard format including metadata and 
#' the character string of the final function evaluation. Third, basic 
#' benchmarking support is available by default, including a timing of the run 
#' with `Sys.time()` and an assessment of memory using `gc()`. These are also 
#' provided in the default returned data. Lastly, results are provided using the 
#' `deconvolution.results` object class, which provides several convenient 
#' methods for handling the outputs (see `?deconvolution.results` for details)
#' 
#' Currently supported deconvolution methods include:
#' 
#' * nnls : Non-negative least squares (NNLS) function from the `nnls` R package
#' (available on CRAN: https://cran.r-project.org/web/packages/nnls/index.html).
#' 
#' * music : The function `music.basic` from the `MuSiC` R package (available on 
#' GitHub: https://github.com/xuranw/MuSiC).
#' 
#' * DeconRNASeq : The function `DeconRNASeq` from the `DeconRNASeq` 
#' R package (available on Bioconductor: 
#' https://doi.org/doi:10.18129/B9.bioc.DeconRNASeq).
#' 
#' * EPIC : The `EPIC` method from the `EPIC` R package (available on GitHub: 
#' https://github.com/GfellerLab/EPIC)
#' 
#' @examples 
#' sce <- random_sce()
#' typev <- unique(sce[["celltype"]])
#' Z <- do.call(cbind, lapply(typev, function(typei){
#' rowMeans(counts(sce[,sce[["celltype"]]==typei]))
#' }))
#' colnames(Z) <- c('type1', 'type2')
#' Y <- matrix(rowMeans(counts(sce)), ncol = 1)
#' 
#' # run methods
#' # run nnls
#' ldecon <- run_deconvolution(method = "nnls", Y = Y, Z = Z)
#' 
#' # inspect results
#' ldecon
#' 
#' @export
run_deconvolution <- function(Z, Y, seed.num = 0, method = "nnls", arguments = list()){
  set.seed(seed.num)
  message("preparing predictions using method ", method, " for :")
  message("G = ", nrow(Z), " marker genes...")
  message("K = ", ncol(Z), " cell types...")
  message("J = ", ncol(Y), " bulk samples...")
  arguments[["Z"]] <- "as.data.frame(Z)"
  arguments[["Y"]] <- "as.data.frame(Y)"
  if(ncol(Y) > 1){
    message("parsing multiple bulk samples...")
    lr <- lapply(seq(ncol(Y)), function(ii){
      arguments[["Y"]] <- Y[,ii,drop=F]
      command.list <- map_deconvolution_arguments(method = method, 
                                                  arguments = arguments)
      get_deconvolution_predictions(command.list = command.list)
    })
    names(lr) <- paste0("bulk_sample", seq(ncol(Y)), ";id:", colnames(Y))
  } else{
      command.list <- map_deconvolution_arguments(method = method,
                                                  arguments = arguments)
      lr <- get_deconvolution_predictions(command.list = command.list)
  }
  return(lr)
}

#' map_deconvolution_arguments
#'
#' Wrapper to manage deconvolution method mapping calls.
#' 
#' @param method Name of deconvolution method.
#' @param arguments List of valid arguments.
#' @returns Result of calling `mapping_[method](...)`
#' @details See `?valid_deconvolution_methods()` for valid methods and their 
#' arguments.
#' 
#' @export
map_deconvolution_arguments <- function(method, arguments){
  command.string <- paste0("map_", method, "(arguments = arguments)")
  return(eval(parse(text = command.string)))
}

#' get_deconvolution_predictions
#'
#' Parses deconvolution method to get proportions predictions, with timing and
#' memory usage for benchmarking.
#' 
#' @param command.list Valid command list, returned from map_deconvolution_arguments().
#' @param item.vector Vector of special item names in command.list. These aren't
#' passed to the deconvolution prediction call.
#' @returns List containing predictions, metadata, and benchmark data.
#'
#' @export
get_deconvolution_predictions <- function(command.list, 
                                          item.vector = c("command.text", 
                                                          "method", "seed.num")){
  # instantiate list objects in environment
  for(name in names(command.list)){
    eval(parse(text = paste0(name, " <- ", command.list[[name]])))
  }
  
  message("getting predictions...")
  method <- command.list[["method"]]
  command.text <- command.list[["command.str"]]
  command.text <- paste0("try(", command.text, ", silent = T)")
  t1 <- Sys.time()
  m1 <- gc(); predictions <- eval(parse(text = command.text)); m2 <- gc()
  duration <- Sys.time() - t1
  if(is(predictions, "try-error")){
    message("error when evaluating deconvolution method")
  } else{
    message("deconvolution method evaluation successful")
  }
  
  # parse return objects
  # parse timing details
  ltime <- list(duration <- Sys.time()-t1, units = "sec", 
                timestamp = as.character(t1), method = "Sys.time()")
  # parse memory usage details
  m.change <- m2-m1
  lmem <- list(start = m1, end = m2, change = m.change, method = "gc")
  lbench <- list(time = ltime, memory = lmem)
  # parse predictions as proportions
  if(is(predictions, "numeric")){
    if(sum(predictions) > 1){predictions <- predictions/sum(predictions)}
  }
  # parse metadata
  lmd <- list(command.text = command.text, method = method, 
              number.of.markers = nrow(Z), number.of.types = ncol(Z),
              markers = rownames(Z), types = colnames(Z), 
              bulk.samples = colnames(Y))
  
  # return
  lreturn <- list(predictions = predictions, metadata = lmd, benchmark = lbench)
  return(as(lreturn, "deconvolution.results"))
}

#-------------------
# main use functions
#-------------------

#' map_nnls
#'
#' Mapping provided arguments to required arguments for NNLS, using defaults
#' for any required arguments not provided
#'
#' @param arguments List of user-provided arguments.
#' @param method Character string of the method.
#' @param library.name Name of library to call function from.
#' @param method.arguments Arguments required for this method.
#' @returns List of data and command character string to parse.
#' 
#' @export
map_nnls <- function(arguments, method = "nnls", library.name = "nnls",
                     method.arguments = c("A" = "Z", "b" = "Y")){
  require(nnls)
  # parse arguments
  message("validating provided arguments...")
  arg.user <- names(arguments)
  arg.method <- names(method.arguments)
  overlapping.args <- intersect(arg.user, arg.method)
  filter.user <- arg.user %in% overlapping.args
  filter.method <- !arg.method %in% overlapping.args
  af.user <- arguments[filter.user]
  af.method <- method.arguments[filter.method]
  message("the following required arguments were provided: ", 
          paste0(names(af.user), collapse = "; "))
  if(length(af.method) > 0){
    message("the following required arguments were not provided: ", 
            paste0(names(af.method), collapse = "; "))
    message("parsing defaults for required methods not provided...")
    for(ai in af.method){
      if(ai == "A"){
        a <- arguments[["Z"]]
      } else if(ai == "b"){
        b <- arguments[["Y"]]
      } else{}
    }
  }
  
  # get the command string
  final.method.vector <- c(af.user, af.method)
  method.string <- paste0(names(final.method.vector), "=", 
                          final.method.vector, collapse = ",")
  command.string <- paste0(library.name, "::", method, "(", method.string, ")$x")
  
  # get final command string in return list
  lr <- lapply(c(af.user, af.method), function(methodi){methodi})
  lr[["command.str"]] <- command.string
  lr[["method"]] <- method
  return(lr)
}

#' map_music
#' 
#' Mapping provided arguments to required arguments for MuSiC, using defaults
#' for any required arguments not provided
#'
#' @param arguments List of user-provided arguments.
#' @param method Character string of the method.
#' @param library.name Name of library to call function from.
#' @param method.arguments Arguments required for this method.
#' @returns List of data and command character string to parse.
#' 
#' @export
map_music <- function(arguments, method = "music.basic", library.name = "MuSiC",
                      method.arguments = c("X" = "Z", "Y" = "Y", "S" = "S", 
                                           "Sigma" = "Sigma", "nu" = "1e-10", 
                                           "iter.max" = "100", "eps" = "0")){
  require(MuSiC)
  require(nnls)
  # parse arguments
  message("validating provided arguments...")
  arg.user <- names(arguments)
  arg.method <- names(method.arguments)
  overlapping.args <- intersect(arg.user, arg.method)
  filter.user <- arg.user %in% overlapping.args
  filter.method <- !arg.method %in% overlapping.args
  af.user <- arguments[filter.user]
  af.method <- method.arguments[filter.method]
  message("the following required arguments were provided: ", 
          paste0(names(af.user), collapse = "; "))
  if(length(af.method) > 0){
    message("the following required arguments were not provided: ", 
            paste0(names(af.method), collapse = "; "))
    message("parsing defaults for required methods not provided...")
    for(ai in af.method){
      if(ai == "Sigma"){
        af.method["Sigma"] <- paste0("matrix(0, ncol = 1, nrow = nrow(Z))")
      } else if(ai == "S"){
        af.method["S"] = paste0("rep(1, ncol(Z))")
      } else if(ai == "nu"){
        af.method["nu"] = paste0("1e-10")
      } else if(ai == "iter.max"){
        af.method["iter.max"] = paste0("1000")
      } else if(ai == "eps"){
        af.method["eps"] = 0
      } else if(ai == "X"){
        af.method["X"] = "Z"
      } else{}
    }
  }
  
  # get the command string
  final.method.vector <- c(af.user, af.method)
  method.string <- paste0(names(final.method.vector), "=", 
                          final.method.vector, collapse = ",")
  command.string <- paste0(library.name, "::", method, 
                           "(", method.string, ")$p.weight")
  
  # get final command string in return list
  lr <- lapply(c(af.user, af.method), function(methodi){methodi})
  lr[["command.str"]] <- command.string
  lr[["method"]] <- method
  return(lr)
}

#' map_deconrnaseq
#'
#' Mapping provided arguments to required arguments for DeconRNASeq, using 
#' defaults for any required arguments not provided
#'
#' @param arguments List of user-provided arguments.
#' @param method Character string of the method.
#' @param library.name Name of library to call function from.
#' @param method.arguments Arguments required for this method.
#' @returns List of data and command character string to parse.
#' 
#' @export
map_deconrnaseq <- function(arguments, method = "DeconRNASeq", 
                            library.name = "DeconRNASeq",
                            method.arguments = c("signatures" = "as.data.frame(Z)",
                                                 "datasets" = "as.data.frame(cbind(Y,Y))",
                                                 "use.scale" = "FALSE")){
  require(DeconRNASeq)
  # parse arguments
  message("validating provided arguments...")
  arg.user <- names(arguments)
  arg.method <- names(method.arguments)
  overlapping.args <- intersect(arg.user, arg.method)
  filter.user <- arg.user %in% overlapping.args
  filter.method <- !arg.method %in% overlapping.args
  af.user <- arguments[filter.user]
  af.method <- method.arguments[filter.method]
  message("the following required arguments were provided: ", 
          paste0(names(af.user), collapse = "; "))
  if(length(af.method) > 0){
    message("the following required arguments were not provided: ", 
            paste0(names(af.method), collapse = "; "))
    message("parsing defaults for required methods not provided...")
    for(ai in names(af.method)){
      if(ai == "signatures"){
        Z <- as.data.frame(arguments[["Z"]])
      } else if(ai == "datasets"){
        Y <- cbind(arguments[["Y"]], arguments[["Y"]])
        Y <- as.data.frame(Y)
      } else{}
    }
  }
  
  # get the command string
  final.method.vector <- c(af.user, af.method)
  filter.vector <- !names(final.method.vector) %in% c("Y", "Z")
  final.method.vector <- final.method.vector[filter.vector]
  method.string <- paste0(names(final.method.vector), "=", 
                          final.method.vector, collapse = ",")
  command.string <- paste0(library.name, "::", 
                           method, "(", method.string, ")$out.all[1,]")
  
  # get final command string in return list
  lr <- lapply(c(af.user, af.method), function(methodi){methodi})
  lr[["command.str"]] <- command.string
  lr[["method"]] <- method
  return(lr)
}

#' map_epic
#'
#' Mapping provided arguments to required arguments for EPIC, using 
#' defaults for any required arguments not provided
#'
#' @param arguments List of user-provided arguments.
#' @param method Character string of the method.
#' @param library.name Name of library to call function from.
#' @param method.arguments Arguments required for this method.
#' @returns List of data and command character string to parse.
#' 
#' @export
map_epic <- function(arguments, method = "EPIC", library.name = "EPIC",
                     method.arguments = c("bulk" = "as.data.frame(Y)",
                                          "reference" = "NA")){
  require(EPIC)
  # parse arguments
  message("validating provided arguments...")
  arg.user <- names(arguments)
  arg.method <- names(method.arguments)
  overlapping.args <- intersect(arg.user, arg.method)
  filter.user <- arg.user %in% overlapping.args
  filter.method <- !arg.method %in% overlapping.args
  af.user <- arguments[filter.user]
  af.method <- method.arguments[filter.method]
  message("the following required arguments were provided: ", 
          paste0(names(af.user), collapse = "; "))
  if(length(af.method) > 0){
    message("the following required arguments were not provided: ", 
            paste0(names(af.method), collapse = "; "))
    message("parsing defaults for required methods not provided...")
    for(ai in names(af.method)){
      if(ai == "bulk"){
        af.method[ai] <- "as.data.frame(Y)"
      } else if(ai == "reference"){
        af.method[ai] <- paste0("list(refProfiles = Z, ",
                                "sigGenes = rownames(Z),",
                                "refProfiles.z = Z)")
      } else{}
    }
  }
  
  # get the command string
  af.user <- af.user[!names(af.user) %in% c("Y", "Z")]
  final.method.vector <- c(af.user, af.method)
  method.string <- paste0(names(final.method.vector), "=", 
                          final.method.vector, collapse = ",")
  command.string <- paste0(library.name, "::", 
                           method, "(", method.string, ")$mRNAProportions")
  
  # get final command string in return list
  lr <- lapply(c(af.user, af.method), function(methodi){methodi})
  lr[["command.str"]] <- command.string
  lr[["method"]] <- method
  return(lr)
}

#' map_bisque
#'
#'
map_bisque <- function(){
  
}

