#!/usr/bin/env R

# Author: Sean Maden
#
# Methods for expression assay matrices (rows = genes/markers, columns = 
# cells/samples/types)
#

#---------------
# anova analysis
#---------------

#' analyze_anova
#'
#' Main function to perform ANOVA analyses on a series of assays, such as from
#' a SingleCellExperiment or SummarizedExperiment object.
#' 
#' @param sce A SingleCellExperiment or SummarizedExperiment object.
#' @param pheno.df Table containing phenotype information. Should include all
#' variables in provided model except "expr". If NULL, use sce colData.
#' @param ngene.sample Number of genes to sample at random.
#' @param assayv Vector of names of assays to analyze.
#' @param model Character string of the model to use. By default, provides an
#' interaction between donor and cell type, and the response/dependent variable
#' is the individual gene.
#' @param return.var Vector of variables to return from ANOVA test results (see 
#' `?get_anova_df` for valid options). By default, returns percentage of sum of
#' squared variances.
#' @param return.sce Whether to return the SingleCellExperiment object with 
#' results stored as new object in metadata. If FALSE, returns results list.
#' @param plot.results Whether to append jitter plots and scatterplots of the
#' results.
#' @param md.name Name of new sce metadata. Ignored if return.sce is FALSE.
#' @param seed.num Random seed to set for reproducibility.
#' @param zoom.panel Whether to show zoomed panels for jitter plots. 
#' @param zoom.ylim The y-axis limmits for zoomed panels on jitter plots.
#' @param verbose Whether to show verbose status messages.
#' @param ... Additional arguments passed to `anova_scatter_plots()` (see 
#' `?anova_scatter_plots` for details).
#' @returns Either a SingleCellExperiment object with results as a new metadata
#' entry, or a list of ANOVA results objects.
#' @seealso get_anova_df, anova_jitter_plots, anova_scatter_plots
#' @examples 
#' # get example data
#' sce <- random_sce()
#' sce[["donor"]] <- c(rep("donor1", 2), rep("donor2", 10))
#' 
#' # simple anova
#' sce <- analyze_anova(sce, model = "expr ~ celltype", ngene.sample = 5)
#' 
#' # anova with batch term -- fails, can't overwrite results metadata
#' sce <- analyze_anova(sce, model = "expr ~ celltype + donor", 
#'                      ngene.sample = 5)
#' 
#' # anova with batch term -- success, writing results to new metadata
#' new.md.name <- "anova_results_batch_model"
#' sce <- analyze_anova(sce, model = "expr ~ celltype + donor", 
#'                      ngene.sample = 5, md.name = new.md.name)
#' 
#' # anova with batch interaction term
#' new.md.name <- "anova_results_batch_interaction_model"
#' sce <- analyze_anova(sce, model = "expr ~ celltype * donor", 
#'                      ngene.sample = 5, md.name = new.md.name)
#' 
#' @export
analyze_anova <- function(sce, pheno.df = NULL, ngene.sample = 1000, 
                          assayv = "counts", model = "expr ~ celltype * donor",
                          return.var = c("percvar"), return.sce = TRUE, 
                          plot.results = TRUE, md.name = "anova_results", 
                          seed.num = 0, zoom.panel = TRUE, zoom.ylim = c(0, 20),
                          verbose = FALSE, ...){
  if(verbose){message("Performing ANOVAs...")}; lr <- list()
  dfa <- get_anova_df(sce = sce, pheno.df = pheno.df, 
                      ngene.sample = ngene.sample, assayv = assayv,
                      model = model, return.var = return.var, 
                      seed.num = seed.num, verbose = verbose)
  lr[["df.anova"]] <- dfa
  if(verbose){message("Getting plots of results...")}
  if(plot.results){
    lggj <- anova_jitter_plots(dfa, zoom.panel = zoom.panel, 
                               zoom.ylim = zoom.ylim)
    lr[["ggplot.jitter"]] <- lggj
    lggpt <- list()
    if(length(assayv) > 1){ # make scatterplots from assay combinations
      lc <- combn(c(assayv), 2) # make combinations
      for(ii in seq(ncol(lc))){ # parse assay combinations
        new.list.name <- paste0("a1:",lc[1,ii],";a2:",lc[2,ii])
        lggpti <- anova_scatter_plots(dfa, a1 = assayv[1], a2 = assayv[2], ...)
        lggpt[[new.list.name]] <- lggpti
      }
      lr[["ggplot.scatterplots"]] <- lggpt
    }
  }
  # parse return option
  if(return.sce){
    if(md.name %in% names(metadata(sce))){
      stop("Error, sce already contains metadata of the same name")}
    metadata(sce)[[md.name]] <- lr
    return(sce)
  }
  return(lr)
}

#' get_anova_df
#'
#' Get a data.frame of ANOVA statistics.
#' 
#' @param sce SingleCellExperiment or SummarizedExperiment object containing
#' at least one expression assay data matrix.
#' @param pheno.df Data.frame containing cell/sample/type phenotype info for the
#' model matrix. Variables named in argument `model` should be listed here, 
#' with the exception of "expr", which is defined from the sce expression data 
#' at runtime.
#' @param ngene.sample Number of genes to sample at random. If NULL, fit every
#' gene in sce.
#' @param model Character string of the model to use. By default, provides an
#' interaction between donor and cell type, and the response/dependent variable
#' is the individual gene.
#' @param seed.num Random seed to set for reproducibility.
#' @param verbose Whether to return verbose status messages.
#' @returns Data.frame containing ANOVA results terms.
#' @examples 
#' sce <- random_sce()
#' sce[["donor"]] <- c(rep("donor1", 2), rep("donor2", 8), rep("donor1", 2))
#' pheno.df <- data.frame(donor = sce[["donor"]], celltype = sce[["celltype"]])
#' dfa <- get_anova_df(sce = sce, pheno.df = pheno.df)
#' 
#' @export
get_anova_df <- function(sce, pheno.df = NULL, ngene.sample = NULL, 
                         assayv = "counts",
                         model = "expr ~ celltype * donor",
                         return.var = c("percvar"), 
                         seed.num = 0, verbose = FALSE){
  set.seed(seed.num); lr <- lgg <- lggj <- lggpt <- list()
  # parse bg gene sampling options
  sampv <- seq(nrow(sce))
  samp.cond <- !is(ngene.sample, "NULL") & ngene.sample < nrow(sce)
  if(samp.cond){sampv <- sample(sampv, ngene.sample, replace = FALSE)}
  # parse pheno.df options
  if(is(pheno.df, "NULL")){pheno.df <- colData(sce)}
  dfa <- do.call(rbind, lapply(assayv, function(ai){
    if(verbose){message("working on assay: ", ai, "...")}
    # parse assay subset
    mi <- assays(sce)[[ai]];mi <- mi[sampv,,drop=F] 
    # filter all-zeros
    maxv <- rowMaxs(mi); filt.max <- maxv > 0
    mi <- mi[filt.max,] 
    # filter nas
    num.na <- apply(mi, 1, function(ri){length(which(is.na(ri)))}) 
    mi <- mi[which(num.na == 0),] 
    dfa.mi <- do.call(rbind, lapply(seq(nrow(mi)), function(ii){
      if(verbose){message("working on gene number ", ii, "...")}
      df.model <- pheno.df; df.model$expr <- mi[ii,]
      av.str <- paste0("aov(formula = ",model,", data = df.model)")
      avi <- eval(parse(text = av.str))
      ## old method:
      # lmi <- lm(expr ~ celltype * donor, data = dfi)
      # avi <- anova(lmi)
      # perc.var <- 100*avi$`Sum Sq`/sum(avi$`Sum Sq`)
      ## use limma:
      # dmat <- model.matrix(~ celltype * donor, data = dfi)
      # avi <- limma::lmFit(object = mi, design = dmat)
      dfi <- as.data.frame(matrix(nrow = 1, ncol = 0))
      namev <- gsub(" ", "", rownames(summary(avi)[[1]]))
      if("percvar" %in% return.var){
        ssqv <- summary(avi)[[1]][[2]] # get sum of squared variances
        perc.var <- 100*ssqv/sum(ssqv)
        for(ii in seq(length(perc.var))){
          dfi[,ncol(dfi) + 1] <- perc.var[ii]
          colnames(dfi)[ncol(dfi)] <- paste0("percvar.", namev[ii])
        }
      }
      dfi
    }))
    dfa.mi$marker <- rownames(mi); dfa.mi$assay <- ai
    dfa.mi <- dfa.mi[!is.na(dfa.mi[,1]),] # filter nas
    dfa.mi
  }))
  if(verbose){message("finished anovas")}
  return(dfa)
}


#' anova_jitter_plots
#'
#' Get jitter scatter plots, grouped by assay, from an anova results data.frame.
#'
#' @param dfa Data.frame containing ANOVA results data, from `get_anova_df()`.
#' @param zoom.panel Whether to include a zoomed panel. Made by calling 
#' `ggforce::facet_zoom()`.
#' @param zoom.ylim Vector of minimum and maximum y-axis limits for zoomed 
#' panel. Ignored if the `zoom.panel` argument is FALSE.
#' @examples 
#' sce <- random_sce()
#' sce[["donor"]] <- c(rep("donor1", 2), rep("donor2", 8), rep("donor1", 2))
#' pheno.df <- data.frame(donor = sce[["donor"]], celltype = sce[["celltype"]])
#' dfa <- get_anova_df(sce = sce, pheno.df = pheno.df)
#' lggj <- anova_jitter_plots(dfa)
#' @seealso anova_scatter_plots, get_anova_df, analyze_anova
#' @export
anova_jitter_plots <- function(dfa, zoom.panel = TRUE, zoom.ylim = c(0, 20)){
  require(ggplot2); require(ggforce)
  # get variables
  cnv <- colnames(dfa)
  cnv <- cnv[!cnv %in% c("marker", "assay")]
  varv <- unique(gsub(".*\\.", "", cnv)) # get variable names
  metricv <- unique(gsub("\\..*", "", cnv)) # get metric names
  if(length(metricv)==0|length(varv)==0){
    stop("Error, couldn't find metrics and variables from dfa colnames.")}
  # get plots list
  lggj <- lapply(metricv, function(mi){
    cnvf <- cnv[grepl(paste0(mi,"\\..*"), cnv)]
    varvf <- unique(gsub(".*\\.", "", cnvf))
    li <- lapply(varvf, function(vi){
      ci <- cnvf[grepl(paste0(".*\\.", vi, "$"), cnvf)]
      title.str <- paste0(ci)
      dfp <- data.frame(value = dfa[,ci], assay = dfa$assay)
      ggj <- ggplot(dfp, aes(x = assay, y = value)) + geom_jitter(alpha = 0.5) + 
        stat_summary(geom = "crossbar", fun = "mean", color = "red") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        ggtitle(title.str)
      if(zoom.panel){ggj <- ggj + facet_zoom(ylim = zoom.ylim)}
      return(ggj)
    })
    names(li) <- varvf
    return(li)
  })
  names(lggj) <- metricv
  return(lggj)
}

#' anova_scatter_plots
#'
#' Get scatter plots showing variances explained from ANOVAs performed on two
#' assay types, a1 (x-axis) and a2 (y-axis)
#' 
#' @param dfa.all Data.frame of ANOVA results, returned from get_anova_df().
#' @param a1 Name of first assay, for x-axis.
#' @param a2 Name of second assay, for y-axis.
#' @param regex.cnv Regex pattern to find column names to plot
#' @returns List of plot objects.
#' @seealso anova_jitter_plots, get_anova_df, analyze_anova
#' @export
anova_scatter_plots <- function(dfa, a1 = "counts", a2 = "counts_ds_combat", 
                                regex.cnv = "^percvar\\..*"){
  # filter dfa
  df1 <- dfa[dfa.all$assay==a1,]; df2 <- dfa[dfa.all$assay==a2,]
  # match dfs
  df1 <- df1[!duplicated(df1$marker),]; df2 <- df2[!duplicated(df2$marker),]
  markerv.int <- intersect(df1$marker, df2$marker)
  df1 <- df1[df1$marker %in% markerv.int,]
  df2 <- df2[df2$marker %in% markerv.int,]
  df2 <- df2[order(match(df2$marker, df1$marker)),]
  # check match success
  cond <- identical(as.character(df2$marker), as.character(df1$marker))
  if(!cond){stop("Error, couldn't match markers for assays.")}
  cnv <- colnames(df1); cnv <- cnv[grepl(regex.cnv, cnv)] # colnames to plot
  lggi <- lapply(cnv, function(ci){ # plot colname matches
    dfp <- data.frame(pv1 = df1[,ci], pv2 = df2[,ci])
    title.str <- paste0("Percent var. ", gsub(".*\\." , "", ci))
    ggplot(dfp, aes(x = pv1, y = pv2)) + geom_point(alpha = 0.5) + 
      geom_abline(slope = 1, intercept = 0, col = "black") +
      ggtitle(title.str) + xlab(a1) + ylab(a2)
  })
  names(lggi) <- cnv
  return(lggi)
}

#---------------------------
# dispersion point estimates
#---------------------------

#' analyze_dispersion_est
#'
#' Perform an analysis of dispersion coefficients by fitting negative binomial
#' models to genes
#'
#' @param sce A SingleCellExperiment object.
#' @param genes.markerv A vector of marker genes.
#' @param type.varname Name of the cell types variable in sce colData columns.
#' @param marker.name.str Character string for marker genes.
#' @param bg.name.str Character string for background genes.
#' @param new.md.name Name of new metadata object in sce. Ignored if 
#' `return.sce` argument is FALSE.
#' @param assay.name Name of assay in sce to analyze.
#' @param num.genes.bg Number of random genes to select for background 
#' estimates.
#' @param seed.num Number to set the seed for computational reproducibility.
#' @param make.plots Whether to make boxplots and jitter plots of the 
#' dispersion coefficients by groups from the specified type variable.
#' @param method.str Method to fit negative binomial models (see 
#' `?mexpr_nbcoef()` for details). Defaults to `glmGamPoi` to use `glm_gp()`.
#' @param return.sce Whether to return a SingleCellExperiment.
#' @param verbose Whether to show verbose status messages.
#' @param ... Additional arguments passed to `plot_dispersion_coeff()`
#' @details Analyze the dispersion coefficients from a SingleCellExperiment 
#' object.
#' 
#' Returned results format is specified by the `return.sce` argument. If this is
#' TRUE, the SingleCellExperiment object is returned with a new metadata item 
#' appended called "dispersion_coefficient_analysis". Otherwise, the results are
#' returned in a list.
#' 
#' @returns Either a SingleCellExperiment object with new metadata, or a list
#' containing analysis results.
#' @seealso plot_dispersion_coeff, get_dispersion_coeff_df
#' @examples 
#' sce <- random_sce()
#' sce <- analyze_dispersion_est(sce, rownames(sce)[seq(5)], "celltype",
#'                               assay.name = "counts", num.genes.bg = 10,
#'                               new.md.name = "dispersion_results")
#' metadata(sce)[["dispersion_results"]]$plots
#' 
#' @export
analyze_dispersion_est <- function(sce, genes.markerv, type.varname = "k2", 
                                   marker.name.str = "top-markers", 
                                   bg.name.str = "bg",
                                   new.md.name = "dispersion_coefficient_analysis",
                                   assay.name = "counts_adj", num.genes.bg = 100,
                                   seed.num = 0, make.plots = TRUE, 
                                   method.str = "glmGamPoi",
                                   return.sce = TRUE, verbose = FALSE, ...){
  ldist <- list()
  # get dispersions by type
  mexpr <- assays(sce)[[assay.name]]; typev <- sce[[type.varname]]
  dfp <- get_dispersion_coef_df(mexpr = mexpr, typev,
                                bg.name.str = bg.name.str,
                                marker.name.str = marker.name.str,
                                num.genes.bg = num.genes.bg, 
                                method.str = method.str,
                                seed.num = seed.num, verbose = verbose)
  ldisp <- list(dfp = dfp) # append results
  # parse plot options
  if(make.plots){
    if(verbose){message("Plotting dispersion point estimates...")}
    lplot <- plot_dispersion_coef(dfp, ...)
    ldisp[["plots"]] <- lplot
  }
  # parse return options
  if(return.sce){
    metadata(sce)[[new.md.name]] <- ldisp
    return(sce)
  } else{
    return(ldisp)
  }
  return(NULL)
}

#' get_dispersion_coef_df
#' 
#' @param mexpr An expression assay matrix (rows = genes/markers, cols = 
#' cells/samples/types).
#' @param typev Vector of length equal to number of columns in mexpr, containing
#' the cell type labels.
#' @param bg.name.str Character string label for background genes.
#' @param marker.name.str Character string label for marker genes.
#' @param num.genes.bg Number of background genes to select at random.
#' @param method.str Method to get dispersion coefficient estimates. This is 
#' passed to `mexpr_nbcoef` (see `?mexpr_nbcoef` for details).
#' @param seed.num Number for random seed.
#' @param verbose Whether to show verbose status messages.
#' @param ... Additional arguments passed to `mexpr_nbcoef()`.
#' @returns Tall data.frame containing dispersion estimates grouped by types. 
#' @seealso plot_dispersion_coef, analyze_dispersion_est
#' @export
get_dispersion_coef_df <- function(mexpr, typev, bg.name.str = "bg", 
                              marker.name.str = "top-markers", 
                              num.genes.bg = 1000, method.str = "glmGamPoi",
                              seed.num = 0, verbose = FALSE, ...){
  set.seed(seed.num)
  # get bg genes
  bg.name <- paste0(bg.name.str, num.genes.bg)
  genes.samplev <- sample(seq(nrow(mexpr)), num.genes.bg, replace = FALSE)
  # define categories
  catv <- c(unique(typev), "all") 
  # get plot data
  if(verbose){message("Calculating dispersion coefficient point estimates...")}
  dfp <- do.call(rbind, lapply(catv, function(typei){
    if(verbose){message("Working on type category: ",typei,"...")}
    # parse filter
    type.filt <- seq(ncol(mexpr))
    if(!typei == "all"){type.filt <- typev == typei}
    mef <- mexpr[,type.filt]
    
    # get dispersions
    if(method.str == "glmGamPoi"){
      lglm.bg <- mexpr_nbcoef(mexpr[genes.samplev,], 
                              method.str = method.str,
                              verbose = verbose, ...)
      lglm.top <- mexpr_nbcoef(mexpr[genes.markerv,], 
                               method.str = method.str,
                               verbose = verbose, ...)
      # get plot data
      dfp1 <- data.frame(disp = lglm.bg$fit$overdispersions)
      dfp1$marker.type <- bg.name
      dfp2 <- data.frame(disp = lglm.top$fit$overdispersions)
      dfp2$marker.type <- marker.name.str
      dfpi <- rbind(dfp1, dfp2); dfpi$celltype <- typei
    } else{
      stop("Error, didn't recognize method.str. ",
           "Must be valid method for mexpr_nbcoef().")
    }
    return(dfpi)
  }))
  return(dfp)
}

#' plot_dispersion_coef
#'
#' Plot dispersion coefficients by gene groups and cell types.
#' 
#' @param dfp Data.frame for plotting.
#' @param make.boxplot Whether to make ggplot boxplots.
#' @param make.jitter Whether to make ggplot jitter plots.
#' @param box.zoom.ymax Vector of y-axis maxima for two zoomed panels.
#' @param jitter.zoom.ymax Vector of y-axiz maxima for two zoomed panels.
#' @param verbose Whether to show verbose status messages.
#' @returns List of ggplot plot objects.
#' @details Makes summary plots of dispersion point estimates. Boxplots and 
#' jitter plots are generated at 3 zoom levels using the `ggforce::facet_zoom()` 
#' function. These are: 1. not zoomed; 2. zoom level 1; and 3. zoom level 2. The 
#' first and second zoom levels are the maximal y-axis values as specified by 
#' the arguments `box.zoom.ymax` and `jitter.zoom.ymax` for boxplots and jitter 
#' plots, respectively.
#' 
#' @seealso analyze_dispersion_est
#' @export
plot_dispersion_coef <- function(dfp, make.boxplot = TRUE, 
                                  make.jitter = TRUE, 
                                  box.zoom.ymax = c(350, 50), 
                                  jitter.zoom.ymax = c(350, 50),
                                  verbose = FALSE){
  require(ggplot2)
  lplot <- list()
  if(make.boxplot){
    if(verbose){message("Getting boxplots...")}
    # boxplots at 3 zoom levels
    lggbox <- list()
    ggbox.zoom1 <- ggplot(dfp, aes(x = marker.type, y = disp)) + 
      geom_boxplot() + facet_wrap(~celltype)
    ggbox.zoom2 <- ggplot(dfp, aes(x = marker.type, y = disp)) + 
      geom_boxplot() + facet_wrap(~celltype) + ylim(0, box.zoom.ymax[1])
    ggbox.zoom3 <- ggplot(dfp, aes(x = marker.type, y = disp)) + 
      geom_boxplot() + facet_wrap(~celltype) + ylim(0, box.zoom.ymax[2])
    lggbox <- list(zoom1 = ggbox.zoom1, zoom2 = ggbox.zoom2, zoom3 = ggbox.zoom3)
    lplot[["ggplot.boxplot"]] <- lggbox
  }
  if(make.jitter){
    # jitter plots at 3 zoom levels
    lggj <- list
    ggj.zoom1 <- ggplot(dfp, aes(x = marker.type, y = disp)) + 
      geom_jitter(alpha = 0.5) + 
      stat_summary(geom = "crossbar", fun = "median", color = "red") + 
      facet_wrap(~celltype)
    ggj.zoom2 <- ggplot(dfp, aes(x = marker.type, y = disp)) + 
      geom_jitter(alpha = 0.5) + 
      stat_summary(geom = "crossbar", fun = "median", color = "red") + 
      facet_wrap(~celltype) + ylim(0, jitter.zoom.ymax[1])
    ggj.zoom3 <- ggplot(dfp, aes(x = marker.type, y = disp)) + 
      geom_jitter(alpha = 0.5) + 
      stat_summary(geom = "crossbar", fun = "median", color = "red") + 
      facet_wrap(~celltype) + ylim(0, jitter.zoom.ymax[2])
    lggj <- list(zoom1 = ggj.zoom1, zoom2 = ggj.zoom2, zoom3 = ggj.zoom3)
    lplot[["ggplot.jitter"]] <- lggj
  }
  return(lplot)
}





