#!/usr/bin/env R

# Author: Sean Maden
#
# Utilities to run deconvolution simulations
#
#
#

#-------------
# dependencies
#-------------
# libv <- c("nnls", "ggplot2", "ComplexHeatmap")
# sapply(libv, library, character.only = T)

#-------------------
# manage simulations
#-------------------

#' supported_strict_methods
#' 
#' Show supported strict deconvolution function. That is, functions which take
#' some reference/signature matrix Z and some bulk/convoluted signals matrix Y
#' and return the vector of predictions or proportions, p.
#' 
#' @param varname Variable/column in lute_transfer_learning.rda dataset with 
#' strict deconvolution method names.
#' @param fname Name of the file containing supported deconvolution methods 
#' info.
#' @returns List of supported strict deconvolution function names and 
#' descriptions.
#' @export 
supported_strict_methods <- function(varname = "strict_deconvolution_method_used",
                                     fname = "lute_transfer_learning"){
  if(verbose){message("Loading data from ",fname,"...")}
  ltl <- get(load(data(fname)))
  if(!varname %in% colnames(ltl)){
    stop("Error, ",varname," is not a variable in dataset.")}
  return(ltl[,varname])
}

#' predtype
#' 
#' Get predictions by type. Takes deconvolution objects and returns vector of 
#' predictions corresponding to identities of signature/reference matrix types 
#' (e.g. returns object pres, where type Z[,1] is pres[1], etc.).
#' 
#' @param Z Signature/reference matrix, where rows correspond to markers and 
#' columns correspond to types.
#' @param Y Convoluted signals matrix, where rows correspond to markers and 
#' columns correspond to types.
#' @param strict_method Name of a supported strict deconvolution method to use.
#' See `supported_strict_methods()` for details about supported methods.
#' @param proportions Whether to return predicted proportions. If False, returns 
#' unmodified point prediction outputs.
#' @param verbose Whether to show verbose status updates.
#' @returns 
#' @examples
#' @seealso decon_results, supported_strict_methods
#' @export
predtype <- function(Z, Y, strict_method = "nnls", proportions = TRUE, 
                     verbose = FALSE){
  if(strict_method == "nnls"){
      p <- nnls::nnls(Z, Y)$x
  } else{
    stop("Error, method not supported. Choose one of either: ",
         paste0(names(supported_strict_methods()), collapse = ","))
  }
  p <- as.numeric(p)
  if(type.prop){
    if(verbose){message("Computing proportions from outputs.")}
    p <- p/sum(p)
    } else{
    if(verbose){message("Returning unmodified point prediction outputs.")}
  }
  return(p)
}


#' decon_results
#' 
#' Run a series of deconvolution simulations and return the results. For list
#' arguments `lgv`, `lpv`, and `lsv`, list index corresponds to the simulation 
#' rep.
#' 
#' @param lgv List of marker expression for reference/signature matrix Z.
#' @param lpv List of type proportions to make Y and compare p predictions.
#' @param lsv List of size factor values. 
#' @param strict_method Type of strict deconvolution method to use (see 
#' `supported_strict_methods()` for details)
#' @param verbose Whether to show verbose status updates.
#' @param 
#' @returns 
#' @examples
#' 
#' # make example data
#' lgv <- list(list(c(1,2),c(2,1),c(1,1)), list(c(2,2),c(2,1),c(1,2)))
#' lpv <- list(c(0.1, 0.8, 0.1),c(0.3, 0.6, 0.1),c(0.2, 0.2, 0.6))
#' lsv <- list(c(1, 10, 10), c(2, 3, 2), c(1, 1, 1))
#' 
#' # run simulations
#' lres <- decon_results(lgv, lpv, lsv)
#' 
#' @seealso
#' @export
decon_results <- function(lgv, lpv, lsv, strict_method = "nnls", 
                          proportions = TRUE, verbose = FALSE){
  if(verbose){message("found ",length(lgv)," expt to run...")}
  lres <- lapply(seq(length(lgv)), function(ii){
    if(verbose){message("Running expt ", ii, " of ",length(lgv),"...")}
    lgi <- lgv[[ii]]; G <- length(lgi[[1]])
    S <- lsv[[ii]]; P <- lpv[[ii]]
    Z <- do.call(cbind, lgi); ZS <- sweep(Z, 2, S, "*")
    Y <- t(t(P) %*% t(ZS))
    if(verbose){message("Getting type predictions...")}
    p1 <- try(predtype(Z = Z, Y = Y, strict_method = strict_method, 
                   proportions = proportions, verbose = verbose))
    p2 <- try(predtype(Z = ZS, Y = Y, strict_method = strict_method, 
                   proportions = proportions, verbose = verbose))
    if(is(p1, "try-error")){
      message("Warning, couldn't get predictions for unadjusted Z test.")}
    if(is(p2, "try-error")){
      message("Warning, couldn't get predictions for S-adjusted Z test.")}
    if(verbose){message("Making result data.frame...")}
    dfres <- do.call(rbind, lapply(list(p1, p2), pdiff, P))
    dfres <- as.data.frame(dfres)
    dfres$expt <- paste0("expt", ii)
    dfres$zs_transform <- c(FALSE, TRUE)
    if(verbose){message("Making results return list...")}
    lexpt <- list(Z = Z, ZS = ZS, Y = Y, method = strict.method)
    lpred <- list(p1 = p1, p2 = p2)
    lres <- list(lexpt = lexpt, lpred = lpred, dfres = dfres)
    return(lres)
  })
  names(lres) <- paste0("expt", seq(length(lres)))
  return(lres)
}

#' decon_analysis
#' 
#' Do simulations, get results df and plots.
#' 
#' @param lgv List of marker expression for reference/signature matrix Z.
#' @param lpv List of type proportions to make Y and compare p predictions.
#' @param lsv List of size factor values. 
#' @param verbose Whether to show verbose status updates.
#' @returns 
#' @examples
#' 
#' # example:
# lgv <- list(list(c(1,2),c(2,1),c(1,1)), list(c(2,2),c(2,1),c(1,2)))
# lpv <- list(c(0.1, 0.8, 0.1),c(0.3, 0.6, 0.1),c(0.2, 0.2, 0.6))
# lsv <- list(c(1, 10, 10), c(2, 3, 2), c(1, 1, 1))
# lres <- decon_results(lgv, lpv, lsv)
#'
#' @seealso decon_results, 
#' @export
decon_analysis <- function(lgv, lpv, lsv, verbose = FALSE){
  # decon_analysis
  #
  # Do simulations, get results df and plots.
  #
  if(verbose){
    message("Running deconvolution simulations...")
    lres <- decon_results(lgv = lgv, lpv = lpv, lsv = lsv, verbose = verbose)
  } else{
    lres <- suppressWarnings(
      decon_results(lgv = lgv, lpv = lpv, lsv = lsv, verbose = verbose))
  }
  if(verbose){message("Appending results data.frames together...")}
  dfres <- do.call(rbind, lapply(lres, function(ii){ii$dfres}))
  if(verbose){message("Appending experiment data to results data.frame...")}
  prop1 <- unlist(lapply(lpv, function(ii){ii[1]}))
  prop2 <- unlist(lapply(lpv, function(ii){ii[2]}))
  sfact1 <- unlist(lapply(lsv, function(ii){ii[1]}))
  sfact2 <- unlist(lapply(lsv, function(ii){ii[2]}))
  dfres$prop_k1 <- rep(prop1, each = 2)
  dfres$prop_k2 <- rep(prop2, each = 2)
  dfres$sfact_k1 <- rep(sfact1, each = 2)
  dfres$sfact_k2 <- rep(sfact2, each = 2)
  if(verbose){message("Making results ggplots...")}
  lgg <- results_plots(dfres = dfres)
  return(list(dfres = dfres, lgg = lgg))
}

#-------------------
# analysis functions
#-------------------
# functions supporting analysis of deconvolution simulation results

#' pdiff
#' 
#' Compare two sets of type predictions across types within an iteration. 
#' 
#' @param pi Prediction set vector. Indices correspond to those in `P`.
#' @param P Prediction set vector. Indices correspond to those in `pi`.
#' @param verbose Whether to show verbose status updates.
#' @returns Comparisons across types, including bias (prediction - true), the
#' root mean squared error (RMSE), and Pearson and Spearman correlation 
#' coefficients.
#' @examples
#' @seealso
#' @export
pdiff <- function(pi, P, verbose = FALSE){
  pi <- as.numeric(pi)
  P <- as.numeric(P)
  bias <- pi - P
  rmse <- sqrt(mean(bias^2))
  corr.p <- try(cor.test(pi, P, method = "pearson")$estimate, silent = verbose)
  corr.s <- try(cor.test(pi, P, method = "spearman")$estimate, silent = verbose)
  corr.p <- ifelse(is(corr.p, "try-error"), NA, corr.p)
  corr.s <- ifelse(is(corr.s, "try-error"), NA, corr.p)
  rv <- c(bias = bias, rmse = rmse, corr.p = corr.p, corr.s = corr.s)
  return(rv)
}

#' results_plots
#' 
#' Makes standard plots to analyze deconvolution simulation results.
#' 
#' @param dfres Data.frame of deconvolution simulation results.
#' @returns List of ggplot2 objects analyzing deconvolution simulation results, 
#' including scatter plots and a violin plot.
#' @param verbose Whether to show verbose status updates.
#' @examples
#' @seealso
#' @export
results_plots <- function(dfres, verbose = FALSE){
  if(verbose){
    message("Making scatter plots of RMSE by first type predictions...")}
  ggpt1 <- ggplot(dfres, aes(x = prop_k1, y = rmse, color = zs_transform)) + 
    geom_point() + ggtitle("RMSE by proportion type 1")
  ggpt1 <- ggpt1 + facet_wrap(~zs_transform)
  if(verbose){message("Making violin plots of RMSE by type...")}
  ggvp <- ggplot(dfres, aes(x = zs_transform , y = rmse, color = zs_transform)) +
    geom_violin(draw_quantiles = 0.5) + ggtitle("RMSE by type")
  if(verbose){
    message("Making scatter plots of RMSE, with vs. without S-transform...")}
  dfp <- data.frame(no_stransform = dfres[dfres$zs_transform==F,]$rmse,
                    with_stransform = dfres[dfres$zs_transform==T,]$rmse)
  ggpt2 <- ggplot(dfp, aes(x = no_stransform, y = with_stransform)) + geom_point() +
    geom_abline(intercept = 0, slope = 1) + ggtitle("RMSE, Z vs. ZS")
  return(list(ggpt1 = ggpt1, ggvp = ggvp, ggpt2 = ggpt2))
}