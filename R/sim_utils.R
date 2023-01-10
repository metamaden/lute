#!/usr/bin/env R

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
#' @returns return 
#' @examples
#' # example
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
  if(proportions){
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
#' @param proportions Whether to predict proportions. If FALSE, returns the 
#' unadjusted point estimates.
#' @param verbose Whether to show verbose status updates.
#' @returns return
#' @examples
#' # make example data
#' lgv <- list(list(c(1,2),c(2,1),c(1,1)), list(c(2,2),c(2,1),c(1,2)))
#' lpv <- list(c(0.1, 0.8, 0.1),c(0.3, 0.6, 0.1),c(0.2, 0.2, 0.6))
#' lsv <- list(c(1, 10, 10), c(2, 3, 2), c(1, 1, 1))
#' # run simulations
#' lres <- decon_results(lgv, lpv, lsv)
#' @seealso decon_analysis
#' @export
decon_results <- function(lgv, lpv, lsv, strict_method = "nnls", 
                          proportions = TRUE, verbose = FALSE){
  if(verbose){message("found ",length(lgv)," expt to run...")}
  # check provided object lengths
  if(length(lpv) >= length(lgv)){
    lpv <- lpv[1:length(lgv)]
  } else{
    stop("Error, not enough lpv entries provided.")
  }
  if(!is(lsv, "NULL")){
    if(length(lsv) >= length(lgv)){
      lsv <- lsv[1:length(lgv)]
    } else{
      stop("Error, provided lsv length should equal lgv length.")
    }
  }
  # run simulations
  lres <- lapply(seq(length(lgv)), function(ii){
    if(verbose){message("Running expt ", ii, " of ",length(lgv),"...")}
    lgi <- lgv[[ii]]
    G <- length(lgi[[1]])
    P <- lpv[[ii]]
    Z <- do.call(cbind, lgi)
    if(verbose){message("Getting type predictions...")}
    if(is(lsv, "NULL")){
      if(verbose){message("Making Y without S transform.")}
      Y <- t(t(P) %*% t(Z))
    } else{
      S <- lsv[[ii]]
      ZS <- sweep(Z, 2, S, "*")
      Y <- t(t(P) %*% t(ZS))
    }
    lexpt <- list(Z = Z, Y = Y, method = strict_method) # metadata for return
    p1 <- try(predtype(Z = Z, Y = Y, strict_method = strict_method, 
                       proportions = proportions, verbose = verbose))
    if(is(p1, "try-error")){
      message("Warning, couldn't get predictions for unadjusted Z test.")}
    lp <- list(p1)
    if(!is(lsv, "NULL")){
      p2 <- try(predtype(Z = ZS, Y = Y, strict_method = strict_method, 
                         proportions = proportions, verbose = verbose))
      if(is(p2, "try-error")){
        message("Warning, couldn't get predictions for S-adjusted Z test.")}
      lp[[2]] <- p2
    }
    dfres <- do.call(rbind, lapply(lp, pdiff, P)) # compute results stats
    dfres <- as.data.frame(dfres)
    dfres$expt <- paste0("expt", ii)
    dfres$zs_transform <- "NA"; dfres$zs_transform[1] <- FALSE
    if(verbose){message("Making results return list...")}
    if(!is(lsv, "NULL")){lexpt[["ZS"]] <- ZS;dfres$zs_transform[2] <- TRUE}
    lpred <- lp; names(lpred) <- paste0("p", seq(length(lp)))
    lres <- list(lexpt = lexpt, lpred = lpred, dfres = dfres)
    return(lres)
  })
  names(lres) <- paste0("expt", seq(length(lres)))
  lres[["lpv"]] <- lpv
  lres[["lgv"]] <- lgv
  lres[["lsv"]] <- lsv
  return(lres)
}

#' decon_analysis
#' 
#' Do simulations, get results df and plots.
#' 
#' @param lpv List of type proportions to make Y and compare p predictions. The
#' length of this list is considered the target number of simulation iterations 
#' for the run.
#' @param lsv List of size factor values. If length(lsv) > length(lpv), only use
#' up to the number of iterations in lpv. If NULL, S-transformation experiments
#' aren't performed, and Y is calculated as $P*Z$ rather than $P*ZS$.
#' @param verbose Whether to show verbose status updates.
#' @param lgv List of marker expression for reference/signature matrix Z. If 
#' length(lgv) > length(lpv), only use up to the number of iterations in lpv.
#' @param sce SingleCellExperiment or SummarizedExperiment object.
#' @param ... Additional arguments passed to `kexpr_sce()`.
#' @returns return
#' @examples
#' # example:
# lgv <- list(list(c(1,2),c(2,1),c(1,1)), list(c(2,2),c(2,1),c(1,2)))
# lpv <- list(c(0.1, 0.8, 0.1),c(0.3, 0.6, 0.1),c(0.2, 0.2, 0.6))
# lsv <- list(c(1, 10, 10), c(2, 3, 2), c(1, 1, 1))
# lres <- decon_results(lgv, lpv, lsv)
#' @seealso decon_results, 
#' @export
decon_analysis <- function(lgv, lpv, lsv, verbose = FALSE, 
                           sce = NULL, ...){
  num.iter <- length(lgv)
  num.types <- length(lpv[[1]])
  if(verbose){message("Prepping ",num.iter," simulation iterations...")}
  # check iterations for each object
  if(verbose){
    message("Running deconvolution simulations...")
    lres <- decon_results(lgv = lgv, lpv = lpv, lsv = lsv, verbose = verbose)
  } else{
    lres <- suppressWarnings(
      decon_results(lgv = lgv, lpv = lpv, lsv = lsv, verbose = verbose))
  }
  if(verbose){message("Appending results data.frames together...")}
  lresf <- lres[grepl("expt", names(lres))]
  dfres <- do.call(rbind, lapply(lres, function(ii){ii$dfres}))
  if(verbose){message("Appending experiment data to results data.frame...")}
  # reassign param lists from results list
  lpv <- lres$lpv
  lsv <- lres$lsv
  lgv <- lres$lgv
  kv <- length(lpv[[1]])
  for(ki in seq(kv)){
    indexv <- as.numeric(gsub("expt", "", dfres$expt)) # index (e.g. for rep)
    new.propv <- unlist(lapply(lpv, function(ii){ii[1]}))[indexv]
    new.sv <- unlist(lapply(lsv, function(ii){ii[1]}))[indexv]
    dfres$newprop <- new.propv; dfres$news <- new.sv
    colnames(dfres)[(ncol(dfres)-1):ncol(dfres)] <- paste0(
      c("prop_k", "sfact_k"), ki)
  }
  if(verbose){message("Making results ggplots...")}
  lgg <- results_plots(dfres = dfres, lsv = lsv)
  lr <- list(dfres = dfres, lgg = lgg)
  if(nrow(dfres) > 3){
    if(verbose){message("Getting by type across simulations...")}
    lr[["dfres.k"]] <- dfres_k(dfres)
  }
  return(lr)
}

#' dfres_k
#'
#' Make k-wise results data.frame from simulation outcomes.
#'
#' @param dfres Results data.frame containing simulation outcomes.
#' @param verbose Whether to show verbose status updates.
#' @returns dfk, data.frame of k/type-wise results summaries.
#' @export
dfres_k <- function(dfres, verbose = FALSE){
  require(dplyr)
  cnv <- colnames(dfres)
  cnv.bias <- cnv[grepl("bias.*", cnv)]
  cnv.prop <- cnv[grepl("prop_.*", cnv)]
  if(!length(cnv.prop)==length(cnv.bias)){
    message("Warning, miss-matched biases and proportions provided.")}
  kv <- gsub("^prop_k", "", cnv.prop)
  ev <- unique(dfres$zs_transform)
  if(verbose){message("Found ",length(kv)," types.")}
  dfk <- do.call(rbind, lapply(ev, function(ei){
    do.call(rbind, lapply(kv, function(ki){
      cnf <- c(cnv.bias[grepl(paste0("bias", ki, "$"), cnv.bias)],
               cnv.prop[grepl(paste0("prop_k",ki,"$"), cnv.prop)])
      filt <- dfres$zs_transform == ei
      dff <- dfres[filt,cnf]
      if(nrow(dff)>3){
        dff$real <- dff[,1] + dff[,2]
        rmse <- sqrt(mean((dff[,2]-dff[,3])^2))
        corr.p <- cor.test(dff[,2], dff[,3], method = "pearson")
        corr.s <- cor.test(dff[,2], dff[,3], mehtod = "spearman")
        return(c(rmse, corr.p$estimate, corr.p$p.value, corr.s$estimate,
                 corr.s$p.value))
      }
    }))
  }))
  dfk <- as.data.frame(dfk)
  colnames(dfk) <- c("rmse", "cor.est.pearson", "cor.pval.pearson",
                    "cor.est.spearman", "cor.pval.spearman")
  dfk$expt <- rep(ev, each = length(kv))
  dfk$k <- rep(kv, length(ev))
  rownames(dfk) <- paste0("k:", dfk$k, ";expt:", dfk$expt)
  return(dfk)
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
#' # examples
#' @seealso decon_analysis, decon_results
#' @export
pdiff <- function(pi, P, verbose = FALSE){
  pi <- as.numeric(pi)
  P <- as.numeric(P)
  if(!length(pi) == length(P)){
    stop("Error, predictions objects must be of the same length.")
  }
  bias <- pi - P
  rmse <- sqrt(mean(bias^2))
  if(length(pi) > 2 & length(P) > 2){
    corr.p <- try(cor.test(pi, P, method = "pearson")$estimate, silent = verbose)
    corr.s <- try(cor.test(pi, P, method = "spearman")$estimate, silent = verbose)
    corr.p <- ifelse(is(corr.p, "try-error"), NA, corr.p)
    corr.s <- ifelse(is(corr.s, "try-error"), NA, corr.p)
  } else{
    if(verbose){message("Two predictions provided; skipping correlations.")}
    corr.p <- corr.s <- NA
  }
  rv <- c(bias = bias, rmse = rmse, corr.p = corr.p, corr.s = corr.s)
  return(rv)
}

#--------------
# results plots
#--------------
# standard plots of simulation results

#' results_plots
#' 
#' Makes standard plots to analyze deconvolution simulation results.
#' 
#' @param dfres Data.frame of deconvolution simulation results.
#' @param lsv List of size factor values. If length(lsv) > length(lpv), only use
#' up to the number of iterations in lpv. If NULL, S-transformation experiments
#' aren't performed, and Y is calculated as $P*Z$ rather than $P*ZS$.
#' @param refline.color Color of the reference line for the second scatterplot 
#' of RMSE by experiment type.
#' @param verbose Whether to show verbose status updates.
#' @returns List of ggplot2 objects analyzing deconvolution simulation results, 
#' including scatter plots and a violin plot.
#' @examples
#' # example
#' @seealso decon_analysis
#' @export
results_plots <- function(dfres, lsv = NULL, refline.color = "black", 
                          verbose = FALSE){
  require(ggplot2); lgg <- list()
  if(is(lsv, "NULL")){
    lgg[["ggpt1"]] <- plot_ggpt_rmsebyp1(dfres = dfres, facet = FALSE, 
                                         verbose = verbose)
    lgg[["ggvp"]] <- plot_ggvp_rmse(dfres = dfres, facet = FALSE, 
                                    verbose = verbose)
    lgg[["ggpt_bias"]] <- plot_ggpt_bias(dfres = dfres, facet = FALSE,
                                         verbose = verbose)
  } else{
    lgg[["ggpt1"]] <- plot_ggpt_rmsebyp1(dfres = dfres, facet = TRUE, 
                                         verbose = verbose)
    lgg[["ggvp"]] <- plot_ggvp_rmse(dfres = dfres, facet = TRUE, 
                                    verbose = verbose)
    lgg[["ggpt_bias"]] <- plot_ggpt_bias(dfres = dfres, facet = TRUE,
                                         verbose = verbose)
    lgg[["ggpt_rmse"]] <- plot_ggpt_bias(dfres = dfres, verbose = verbose)
  }
  return(lgg)
}

#' plot_ggpt_rmsebyp1
#'
#' Makes a ggplot scatterplot of the RMSE (y-axis) by the true proportion of 
#' type 1 (x-axis).
#'
#' @param dfres Simulation series results data.frame.
#' @param facet Whether to facet scatterplots for each S condition, if present.
#' @param verbose Whether to return verbose status messages.
#' @returns ggplot scatterplot object
#' @export
plot_ggpt_rmsebyp1 <- function(dfres, facet = TRUE, verbose = FALSE){
  if(verbose){
    message("Making scatter plots of RMSE by first type predictions...")}
  ggpt <- ggplot(dfres, aes(x = prop_k1, y = rmse)) + 
    geom_point() + ggtitle("RMSE by proportion type 1")
  if(facet){ggpt <- ggpt + facet_wrap(~zs_transform)}
  return(ggpt)
}

#' plot_ggpt_bias
#'
#' Makes a ggplot scatterplot of the true proportions (x-axis) by the predicted
#' proportions (y-axis), where point colors and shapes indicate the type.
#'
#' @param dfres Simulation series results data.frame.
#' @param facet Whether to facet scatterplots for each S condition, if present.
#' @param verbose Whether to return verbose status messages.
#' @returns ggplot scatterplot object
#' @export
plot_ggpt_bias <- function(dfres, facet = TRUE, verbose = FALSE){
  if(verbose){
    message("Making bias scatter plots, with vs. without S-transform...")}
  # get plot data
  cnv <- colnames(dfres)
  typev <- gsub("^bias", "", cnv[grepl("^bias.*", cnv)])
  dfp <- do.call(rbind, lapply(typev, function(typei){
    prop_truei <- dfres[,paste0("prop_k", typei)]
    prop_predi <- dfres[,paste0("bias", typei)] + prop_true
    data.frame(prop_true = prop_truei, 
               prop_pred = prop_predi,
               expt_type = dfres$zs_transform,
               celltype = rep(typei, nrow(dfres)))
  }))
  # get new expt lvl labels
  lvlstr.false <- "no_scaling"; lvlstr.true <- "with_scaling"
  # get rmse to print
  rmse.false <- sqrt(mean((dfp[dfp$expt_type==FALSE,]$prop_true-
                             dfp[dfp$expt_type==FALSE,]$prop_pred)^2))
  rmse.true <- sqrt(mean((dfp[dfp$expt_type==TRUE,]$prop_true-
                            dfp[dfp$expt_type==TRUE,]$prop_pred)^2))
  dfp$rmse <- ifelse(dfp$expt_type==TRUE, rmse.true, rmse.false)
  df.rmse <- data.frame(expt_type = c(lvlstr.false, lvlstr.true),
                        rmse = c(format(rmse.false, digits = 2), 
                                 format(rmse.true, digits = 2)))
  df.rmse$xpos <- min(dfp$prop_true)+0.1*min(dfp$prop_true)
  df.rmse$ypos <- max(dfp$prop_pred)-0.1*max(dfp$prop_pred)
  df.rmse$hjustpos <- df.rmse$vjustpos <- 0
  df.rmse$rmse <- paste0("RMSE: ", df.rmse$rmse)
  # format expt_type variable
  dfp$expt_type <- ifelse(dfp$expt_type == "TRUE", lvlstr.true, lvlstr.false)
  dfp$expt_type <- factor(dfp$expt_type, levels = c(lvlstr.false, lvlstr.true))
  # new plot object
  ggpt <- ggplot() + theme_bw() +
    geom_text(data = df.rmse, alpha = 0.8,
              mapping = aes(x = xpos, y = ypos, 
                            hjust = hjustpos, vjust = vjustpos,
                            label = rmse)) +
    geom_point(dfp, mapping = aes(x = prop_true, y = prop_pred, 
                                  shape = celltype, color = celltype),
               alpha = 0.5, size = 3) + 
    geom_abline(intercept = 0, slope = 1, col = "black") +
    xlim(0.38, 1) + ylim(0.43, 1) +
    xlab("True cell composition (cc)") +
    ylab("Estimated cc")
  if(facet){ggpt <- ggpt + facet_grid(cols=vars(expt_type))}
  return(ggpt)
}

#' plot_ggvp_rmse
#'
#' Makes a violin plot of the RMSE (y-axis) by experiment condition (x-axis).
#'
#' @param dfres Simulation series results data.frame.
#' @param facet Whether to facet violin plots for each S condition, if present.
#' @param verbose Whether to return verbose status messages.
#' @returns ggplot violin plot object
#' @export
plot_ggvp_rmse <- function(dfres, facet = TRUE, verbose = FALSE){
  if(verbose){message("Making violin plots of RMSE by type...")}
  ggvp <- ggplot(dfres, aes(x = zs_transform , y = rmse)) +
    geom_violin(draw_quantiles = 0.5) + ggtitle("RMSE by type")
  if(facet){ggvp <- ggvp + facet_wrap(~zs_transform)}
  return(ggvp)
}

#' plot_ggpt_rmse
#'
#' Plots RMSE for zs_transform == FALSE (x-axis) by zs_transform == TRUE 
#' (y-axis). Note, this requires that the S factor condition varied within
#' simulations.
#'
#' @param dfres Simulation series results data.frame.
#' @param verbose Whether to return verbose status messages.
#' @returns ggplot scatterplot object
#' @export
plot_ggpt_rmse <- function(dfres, verbose = FALSE){
  if(verbose){
    message("Making scatter plots of RMSE, with vs. without S-transform...")}
  dfp <- data.frame(no_stransform = dfres[dfres$zs_transform==F,]$rmse,
                    with_stransform = dfres[dfres$zs_transform==T,]$rmse)
  ggpt <- ggplot(dfp, aes(x = no_stransform, y = with_stransform)) + 
    geom_point() + geom_abline(intercept = 0, slope = 1, color = "black") + 
    ggtitle("RMSE, Z vs. ZS")
  return(ggpt)
}