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
#' @param ltl.varname Variable/column in lute_transfer_learning.rda dataset 
#' containing the strict deconvolution method names.
#' @returns List of supported strict deconvolution function names and descriptions.
#' @export 
supported_strict_methods <- function(ltl.varname = "strict_deconvolution_method_used"){
  ltl <- get(load(data("lute_transfer_learning")))
  return(ltl[,ltl.varname])
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
#' @param proportions Whether to return predicted proportions. If False, returns 
#' unmodified point prediction outputs.
#' @param verbose Whether to show verbose status updates.
#' @returns 
#' @examples
#' @seealso decon_results, lute_methods
#' @export
predtype <- function(Z, Y, strict.method = "nnls", proportions = TRUE, verbose = FALSE){
  if(method == "nnls"){p <- nnls(Z, Y)$x} else{
    stop("Error, method not supported. Choose one of either: ".
         paste0(names(supported_strict_methods()), collapse = ","))
  }
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
#' @param strict.method Type of strict deconvolution method to use (see 
#' `supported_strict_methods()` for details)
#' @param 
#' @returns 
#' @examples
#' @seealso
#' @export
decon_results <- function(lgv, lpv, lsv, strict.method = "nnls", 
                          type.prop = proportions, verbose = F){
  # decon_results
  #
  # example:
  # lgv <- list(list(c(1,2),c(2,1),c(1,1)), list(c(2,2),c(2,1),c(1,2)))
  # lpv <- list(c(0.1, 0.8, 0.1),c(0.3, 0.6, 0.1),c(0.2, 0.2, 0.6))
  # lsv <- list(c(1, 10, 10), c(2, 3, 2), c(1, 1, 1))
  # lres <- decon_results(lgv, lpv, lsv)
  #
  if(verbose){message("found ",length(lgv)," expt to run...")}
  lres <- lapply(seq(length(lgv)), function(ii){
    if(verbose){message("running expt ", ii, " of ",length(lgv),"...")}
    lgi <- lgv[[ii]]; G <- length(lgi[[1]])
    S <- lsv[[ii]]; P <- lpv[[ii]]
    Z <- do.call(cbind, lgi); ZS <- sweep(Z, 2, S, "*")
    Y <- t(t(P) %*% t(ZS))
    p1 <- predtype(Z = Z, Y = Y, strict.method = strict.method, 
                   proportions = proportions)
    p2 <- predtype(Z = ZS, Y = Y, strict.method = strict.method, 
                   proportions = proportions)
    dfres <- do.call(rbind, 
                     lapply(list(p1, p2), 
                            pdiff, P)) # make results df
    dfres <- as.data.frame(dfres)
    dfres$expt <- paste0("expt", ii)
    dfres$zs_transform <- c(FALSE, TRUE)
    # make return list
    lexpt <- list(Z = Z, ZS = ZS, Y = Y, method = "nnls")
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
#' @param lgv
#' @param lpv 
#' @param lsv
#' @returns 
#' @examples
#' # example:
# lgv <- list(list(c(1,2),c(2,1),c(1,1)), list(c(2,2),c(2,1),c(1,2)))
# lpv <- list(c(0.1, 0.8, 0.1),c(0.3, 0.6, 0.1),c(0.2, 0.2, 0.6))
# lsv <- list(c(1, 10, 10), c(2, 3, 2), c(1, 1, 1))
# lres <- decon_results(lgv, lpv, lsv)
#' @seealso decon_results, 
#' @export
decon_analysis <- function(lgv, lpv, lsv){
  # decon_analysis
  #
  # Do simulations, get results df and plots.
  #
  lres <- suppressWarnings(decon_results(lgv, lpv, lsv, verbose = F))
  dfres <- do.call(rbind, lapply(lres, function(ii){ii$dfres}))
  # append proportions for k1
  prop1 <- unlist(lapply(lpv, function(ii){ii[1]}))
  prop2 <- unlist(lapply(lpv, function(ii){ii[2]}))
  sfact1 <- unlist(lapply(lsv, function(ii){ii[1]}))
  sfact2 <- unlist(lapply(lsv, function(ii){ii[2]}))
  dfres$prop_k1 <- rep(prop1, each = 2)
  dfres$prop_k2 <- rep(prop2, each = 2)
  dfres$sfact_k1 <- rep(sfact1, each = 2)
  dfres$sfact_k2 <- rep(sfact2, each = 2)
  lgg <- results_plots(dfres = dfres)
  return(list(dfres = dfres, lgg = lgg))
}

#---------------------
# comparator functions
#---------------------

#'
#' @param 
#' @returns 
#' @examples
#' @seealso
#' @export
pdiff <- function(pi, P){
  bias <- pi - P
  rmse <- sqrt(mean(bias^2))
  corr.p <- try(
    cor.test(pi, P, method = "pearson")$estimate,
    silent = T)
  corr.p <- ifelse(is(corr.p, "try-error"), NA, corr.p)
  corr.s <- try(
    cor.test(pi, P, method = "spearman")$estimate,
    silent = T)
  corr.s <- ifelse(is(corr.s, "try-error"), NA, corr.p)
  rv <- c(bias = bias, rmse = rmse, 
          corr.p = corr.p, 
          corr.s = corr.s)
  return(rv)
}

results_plots <- function(dfres){
  # plot rmse by proportion k1
  ggpt1 <- ggplot(dfres, aes(x = prop_k1, y = rmse, color = zs_transform)) + 
    geom_point() + ggtitle("RMSE by proportion type 1")
  ggpt1 <- ggpt1 + facet_wrap(~zs_transform)
  ggvp <- ggplot(dfres, aes(x = zs_transform , y = rmse, color = zs_transform)) +
    geom_violin(draw_quantiles = 0.5) + ggtitle("RMSE by type")
  dfp <- data.frame(no_stransform = dfres[dfres$zs_transform==F,]$rmse,
                    with_stransform = dfres[dfres$zs_transform==T,]$rmse)
  ggpt2 <- ggplot(dfp, aes(x = no_stransform, y = with_stransform)) + geom_point() +
    geom_abline(intercept = 0, slope = 1) + ggtitle("RMSE, Z vs. ZS")
  return(list(ggpt1 = ggpt1, ggvp = ggvp, ggpt2 = ggpt2))
}