#!/usr/bin/env R

# Author: Sean Maden
#
# Methods for SingleCellExperiment objects.
#

#-------------------
# summary statistics
#-------------------
#' sce_groupstat
#' 
#' Get group-level summary statistics, either on rowData or collapsed colData.
#' 
#' @param sce Filtered `SingleCellExperiment` object. Should reflect the data
#' for a specific type.
#' @param group.variable Variable containing group information on scef colData.
#' @param summarytype Whether to summarize data on rowData (e.g. one entry per 
#' row/gene), or otherwise use the colData collapsed on means (e.g. take the 
#' mean expression across cells for each gene, returning one row per type).
#' @param groupstat Summary statistics to calculate. Can be either of "count", 
#' "mean", "median", "var", "sd", or "numzero" (i.e. number of zero-value 
#' entries).
#' @param type.variable Variable containing type labels.
#' @param return.tall Whether to return a tall table. If FALSE, returns the 
#' rowData/colData-compatible format with groups in column names.
#' @param verbose Whether to show verbose status messages.
#' @returns data.frame containing group-level summary statistics for all groups 
#' specified in ugroupv.
#' @details Computes summary statistics for either rowData (e.g. genes, markers,
#' etc.) or colData (e.g. samples, cells, etc.) for an object of class 
#' SingleCellExperiment.
#' 
#' Dimensions in the returned object are determined by the argument 
#' `return.tall`. When this is TRUE, rows in the returned table correspond to
#' donors or donor-types (when `type.variable` is provided) if 
#' `summarytype`=="colData", and donor-genes otherwise. When this is FALSE, a 
#' colData/rowData-compatible object is returned as specified by `summarytype`.
#'  
#' @examples 
#' 
#' # make random data
#' sce = random_sce(zero.include = TRUE, zero.fract = 0.5)
#' colData(sce)$donor <- c(rep("donor1", 7), rep("donor2", 3))
#'
#' # get group summaries across types
#' sce_groupstat(sce, summarytype = "colData", return.tall = TRUE)
#'
#' # get group-gene summaries across types
#' sce_groupstat(sce, summarytype = "rowData", return.tall = TRUE)
#'
#' # get tall table of group-type summaries
#' sce_groupstat(sce, type.variable = "celltype", 
#'              summarytype = "colData", return.tall = TRUE)
#'
#' # get coldata-compatible group-type summaries
#' sce_groupstat(sce, type.variable = "celltype", 
#'              summarytype = "colData", return.tall = FALSE)
#'
#' # get rowdata-compatible group-type summaries
#' sce_groupstat(sce, type.variable = "celltype", 
#'              summarytype = "rowData", return.tall = FALSE)
#' 
#' @export
sce_groupstat <- function(sce, group.variable = "donor",  
                          assayname = "counts", summarytype = "colData", 
                          groupstat = c("count", "numzero", "var"),
                          return.tall = FALSE, type.variable = NULL, 
                          verbose = FALSE){
  # run checks
  # check sce
  if(!is(sce, "SingleCellExperiment")){
    stop("Error, sce must be a SingleCellExperiment.")}
  # check assay signals
  expr <- assays(sce)[[assayname]]
  cond <- is(expr, "matrix")|is(expr, "DelayedArray")
  if(!cond){
    stop("Error, assay signals should be either a matrix or DelayedArray.")}
  # filter group stats
  groupstat.filt <- groupstat %in% c("count", "mean", "median", 
                                     "var", "sd", "numzero")
  groupstat <- groupstat[groupstat.filt]
  if(verbose){message("Checking colData variables...")}; cd <- colData(sce)
  lvar <- check_coldata(cd = cd, var = c(group.variable, type.variable))
  # get group stats df
  ugroupv <- lvar[[group.variable]]$uvec
  groupv <- lvar[[group.variable]]$vec
  if(!is(type.variable, "NULL")){
    utypev <- lvar[[type.variable]]$uvec
    typev <- lvar[[type.variable]]$vec
  }
  
  ldf <- lapply(ugroupv, function(groupi){
    if(is(type.variable, "NULL")){
      exprf <- assays(sce[,groupv==groupi])[[assayname]]
      dfi <- get_groupstat_df(exprf, groupstat = groupstat, 
                              summarytype = summarytype)
      
    } else{
      dfi <- do.call(rbind, lapply(utypev, function(typei){
        sce.filt <- groupv==groupi & typev==typei
        exprff <- assays(sce[,sce.filt])[[assayname]]
        dfii <- get_groupstat_df(exprff, groupstat = groupstat, 
                                 summarytype = summarytype)
        dfii$type <- typei; return(dfii)
      }))
    }
    dfi$group <- groupi
    dfi
  })
  if(return.tall){
    dfr <- do.call(rbind, lapply(ldf, function(dfi){dfi}))
  } else{
    dfr <- do.call(cbind, lapply(ldf, function(dfi){
      groupi <- unique(dfi$group)[1]
      if(summarytype == "rowData" & !is(type.variable, "NULL")){
        dfi <- do.call(cbind, lapply(utypev, function(typei){
          dfii <- dfi[dfi$type == typei,,drop = F]
          dfii <- dfii[,!colnames(dfii) %in% c("group", "marker", "type")]
          colnames(dfii) <- paste0(groupi, ";", typei, ";", colnames(dfii))
          dfii
        }))
      } else{
        dfi <- dfi[,!colnames(dfi) %in% c("group", "type", "marker")]
        colnames(dfi) <- paste0(groupi, ";", colnames(dfi))
      }
      dfi
    }))
    cond <- summarytype == 'colData' & !is(type.variable, "NULL")
    if(cond){dfr$type <- utypev}
  }
  return(dfr)
}

#' check_coldata
#'
#' Check coldata columns.
#' 
#' @param cd Matrix of colData.
#' @param varv Vector of variable names to check
#' @returns list, contains vector of variable values and vector of unique levels
#' for variables found in cd.
#' @export
check_coldata <- function(cd, varv){
  varv <- varv[!is.null(varv)]
  lr <- lapply(varv, function(vari){
    if(!vari %in% colnames(cd)){
      stop("Error, didn't find variable ", vari, "in colData colnames")}
    lri <- list()
    lri[["vec"]] <- cd[,vari]
    lri[["uvec"]] <- unique(cd[,vari])
    return(lri)
  })
  names(lr) <- varv
  return(lr)
}

#' get_groupstat_df
#'
#' Get grouped statistics from an expression matrix.
#' 
#' @param exprf Filtered expression matrix.
#' @param groupstat Vector of valid group statistics.
#' @param summarytype Type of summary. Either colData or rowData.
#' @param round.digits Digits to round.
#' @returns Table of type data.frame containing summary statistics.
#' @examples
#' sce = random_sce(zero.include = TRUE, zero.fract = 0.5)
#' colData(sce)$donor <- c(rep("donor1", 7), rep("donor2", 3))
#' exprf <- assays(sce)$counts
#' get_groupstat_df(exprf)
#' 
#' @export
get_groupstat_df <- function(exprf, groupstat = c("count", "var"), 
                             summarytype = "colData", round.digits = 3){
  # parse summarytype
  ngene <- nrow(exprf); ncell <- ncol(exprf); exprff <- exprf
  if(summarytype == "colData"){exprff <- t(as.matrix(colMeans(exprf)))}
  # make na matrix
  groupstatf <- groupstat[!grepl("^count.*|^numzero.*", groupstat)]
  mna <- matrix(NA, ncol = length(groupstatf), nrow = nrow(exprff))
  dfti <- as.data.frame(mna); colnames(dfti) <- groupstatf
  cond <- length(which(grepl("^count.*", groupstat)))
  if(cond > 0){
    dfti$count.cells <- ncell; dfti$count.genes <- ngene
  }
  for(ri in seq(nrow(exprff))){
    if(!is(exprff, "NULL")){
      # get group stats
      ei <- exprff[ri,,drop = F]
      if("mean" %in% groupstat){
        meani <- rowMeans(ei)
        dfti[ri,]$mean <- round(meani, digits = round.digits)
      }
      if("median" %in% groupstat){
        mediani <- rowMedians(ei)
        dfti[ri,]$median <- round(mediani, digits = round.digits)
      }
      if("var" %in% groupstat){
        vari <- rowVars(ei)
        dfti[ri,]$var <- round(vari, digits = round.digits)
      }
      if("sd" %in% groupstat){
        sdi <- rowSds(ei)
        dfti[ri,]$sd <- round(sdi, digits = round.digits)
      }
    }
  }
  # parse zero counts
  cond <- length(which(grepl("^numzero.*", groupstat)))
  if(cond > 0){
    if(summarytype == "colData"){
      # get summaries across cells
      zct.cell <- median(apply(exprf, 1, function(ri){length(which(ri==0))}))
      zfr.cell <- median(zct.cell/ncell)
      dfti$median.count.cells.zero <- round(zct.cell, digits = round.digits)
      dfti$median.fract.cells.zero <- round(zfr.cell, digits = round.digits)
      # get summaries across genes
      zct.gene <- median(apply(exprf, 2, function(ci){length(which(ci==0))}))
      zfr.gene <- median(zct.gene/ngene)
      dfti$median.count.genes.zero <- round(zct.gene, digits = round.digits)
      dfti$median.fract.genes.zero <- round(zfr.gene, digits = round.digits)
    } else{
      dfti$count.cells.zero <- apply(exprf, 1, function(ri){
        length(which(ri == 0))})
    }
  }
  if(summarytype == "rowData"){dfti$marker <- rownames(exprf)}
  return(dfti)
}

#------
# plots
#------

#' mexpr_downsample
#'
#' Perform downsampling on an expression assay matrix (rows = markers/genes,
#' columns = cells/samples).
#' 
#' @param mexpr Expression assay matrix.
#' @param ds.method Name of the downsampling method to use. The default option
#' "scuttle" uses `scuttle::downsampleMatrix()`.
#' @param as.matrix Whether to convert the downsampled expressiond data to
#' a matrix object.
#' @param verbose Whether to show verbose status messages.
#' @returns Matrix containing the downsampled expression assays matrix.
#' @details This function provides a wrapper to call downsampling functions for 
#' an expression assay matrix. The input matrix show have rows corresponding to
#' genes or markers, and columns corresponding to cells or samples, and it 
#' should be a matrix-type object.
#' 
#' Downsampling is a technique to select subsets of expression data based on the 
#' minimum batch-wise observed expression. This can help mitigate batch effects 
#' in transcriptomics datasets. The operation is a function of the probability
#' density function or proportions specified for each cell/sample/column in
#' mexpr. By default, this function uses the formula of min(colsums)/colsums to
#' define this vector of probabilty densities.
#' 
#' @seealso sce_dispersion, get_dispersion_data, plot_dispersion
#' @export
mexpr_downsample <- function(mexpr, ds.method = "scuttle", as.matrix = TRUE,
                             verbose = FALSE){
  if(!is(mexpr, "matrix")){stop("Error, mexpr should be a matrix.")}
  if(verbose){message("Getting the proportions vector...")}
  sumv <- colSums(mexpr); propv <- min(sumv)/sumv
  if(ds.method == "scuttle"){
    require(scuttle)
    if(verbose){message("Using downsampling method ",ds.method,"...")}
    mexpr <- scuttle::downsampleMatrix(mexpr, prop = propv, bycol = T)
  } else{
    stop("Error, didn't recognize ds.method.")
  }
  if(as.matrix){mexpr <- as.matrix(mexpr)}
  if(verbose){message("Completed downsampling. Returning...")}
  return(mexpr)
}

#' mexpr_nbcoef
#'
#' Computes the negative binomial distribution coeffcients from an expression
#' assay matrix
#'
#' @param mexpr An expression assay matrix (rows = genes/markers, columns = 
#' cells/samples).
#' @param method.str Method used to fit the negative binomial model and compute 
#' distribution coefficients. If set to "gmlGamPoi" (the default), calls
#' `glmGamPoi::glm_gp()` for each gene/row in mexpr.
#' @param verbose Whether to show verbose status messages.
#' @returns Results object containing the model coefficients for each gene and
#' group in mexpr.
#' @details Fits a negative binomial model for each gene in an expression assay
#' matrix, returning the distribution coefficients. These are useful for running
#' simulations using empirical estimates of the distribution coefficients, and
#' for comparing gene-wise expression distributions. This function is called by 
#' `sce_dispersion()` (see `?sce_dispersion` for details).
#' 
#' @seealso sce_dispersion
#' @export
mexpr_nbcoef <- function(mexpr, method.str = "glmGamPoi", verbose = FALSE){
  if(verbose){
    message("Getting negative binomial model coefficients from mexpr.")}
  lr <- list()
  if(method.str == "glmGamPoi"){
    require(glmGamPoi)
    if(verbose){message("Using method ", method.str, 
                        " to get model coefficients...")}
    lr[["fit"]] <- glmGamPoi::glm_gp(mexpr)
  } else{
    stop("Error, didn't recognize method.")
  }
  return(lr)
}

#' sce_dispersion
#'
#' Perform dispersion analysis for a SingleCellExperiment object.
#' 
#' @param expr.data Either a matrix, a SingleCellExperiment, or 
#' SummarizedExperiment object.
#' @param group.data Either a colData variable name (if expr.data is of type 
#' SingleCellExperiment or SummarizedExperiment), or a vector of labels of 
#' length equal to the number of columns in expr.data (if expr.data is a 
#' matrix).
#' @param plot.dispersion Whether to make plots of dispersion using 
#' `plot_dispersion()` (see `?plot_dispersion` for details).
#' @param highlight.markers Vector of marker identifiers to highlight in 
#' dispersion plots. Ignored if NULL (default).
#' @param assayname Name of the expression assay matrix. Only used if expr.data
#' is a SingleCellExperiment or SummarizedExperiment.
#' @param downsample Whether to perform downsampling by calling 
#' `mexpr_downsample()` (see `?mexpr_downsample` for details).
#' @param ds.method Name of the downsampling method to use (only used if 
#' `downsample` is `TRUE`).
#' @param get.nbstat Compute the negative binomial distribution coefficients 
#' for each gene/marker. Calls the function `mexpr_nbcoef()` (see 
#' `?mexpr_nbcoef` for details).
#' @param method.nbstat Method to compute the negative binomial distribution
#' coefficients. Should be a valid method recognized by the function 
#' `mexpr_nbcoef()` (see `?mexpr_nbcoef` for details).
#' @param verbose Whether to show verbose status messages.
#' @param ... Additional arguments passed to `plot_dispersion()`.
#' @returns List of dispersion results, including the summary statistics table
#' and the dispersion plot.
#' @details This is the main function to manage a dispersion analysis for a 
#' SingleCellExperiment, SummarizedExperiment, or matrix of expression assay
#' data. This generates and checks the group-level summary statistics, then 
#' produces plots of gene/marker-wise points or model smooths showing the 
#' mean and variance relationship. The coefficients for the negative binomial 
#' distribution can also be generated for each gene in each provided group.
#' 
#' @seealso get_dispersion_data, plot_dispersion, mexpr_downsample, mexpr_nbcoef
#' @examples 
#' # get some random example sce data
#' sce.exe <- random_sce()
#' 
#' # analyze dispersion
#' #
#' # get dispersion results with defaults
#' ld1 <- sce_dispersion(sce.exe, group.data = "celltype", verbose = TRUE)
#' #
#' # show marker highlights in plots
#' ld2 <- sce_dispersion(sce.exe, group.data = "celltype", verbose = TRUE,
#' highlight.markers = rownames(sce.exe)[seq(10)])
#' # 
#' # get negative binomial coefficients
#' ld3 <- sce_dispersion(sce.exe, get.nbstat = TRUE)
#'
#' @export
sce_dispersion <- function(expr.data, group.data = NULL, plot.dispersion = TRUE,
                           highlight.markers = NULL, assayname = "counts", 
                           downsample = TRUE, ds.method = "scuttle",
                           get.nbstat = FALSE, method.nbstat = "glmGamPoi",
                           verbose = FALSE, ...){
  lr <- list() # begin return list
  lr[["metadata"]] <- list(expr.data.class = class(expr.data)[1], 
                           assayname = assayname, 
                           downsample = list(downsample, ds.method),
                           nbstat = list(get.nbstat, method.nbstat),
                           markers = highlight.markers)
  cond.se <- is(expr.data, "SingleCellExperiment")|
    is(expr.data, "SummarizedExperiment")
  group.vector <- NULL
  if(cond.se){
    sce <- expr.data
    mexpr <- assays(sce)[[assayname]]
    mexpr <- as.matrix(mexpr)
    if(!is(group.data, "NULL")){
      group.vector <- as.character(sce[[group.data]])
    }
  } else{
    mexpr <- expr.data
    if(!is(group.data, "NULL")){
      group.vector <- group.data
    }
  }
  cond.matrix <- is(mexpr, "matrix")
  if(!cond.matrix){
    stop("Error, couldn't get matrix of expression data from expr.data.")}
  ds.str <- "" # set plot title character string variable
  if(downsample){
    if(verbose){message("Performing downsampling...")}
    mexpr <- mexpr_downsample(mexpr = mexpr, ds.method = ds.method, 
                              verbose = verbose)
    ds.str <- "d.s." # update plot title character string variable
  }
  if(get.nbstat){
    if(verbose){message("Computing the negative binomial coefficients...")}
    nbstat <- mexpr_nbcoef(mexpr, method.str = method.nbstat, 
                           verbose = verbose)
    lr[["neg.binom.statistic"]] <- nbstat
  }
  dfstat <- get_dispersion_data(mexpr = mexpr, group.vector = group.vector)
  lr[["dfstat"]] <- dfstat
  title.str <- paste0("Dispersion, ",ds.str, " ", assayname)
  plot.obj <- plot_dispersion(dfstat = dfstat, title.str = title.str,
                              highlight.markers = highlight.markers, ...)
  lr[["ggplot.dispersion"]] <- plot.obj
  return(lr)
}

#' get_dispersion_data
#'
#' Get data for dispersion analyses and plots.
#' 
#' @param dfstat Table containing statistics for dispersion analysis. Should 
#' include columns for "mean" and "var".
#' @param mexpr Matrix of expression assay data (rows = genes/markers, 
#' columns = cells/samples).
#' @param sce SingleCellExperiment object containing a assays expression matrix
#' specified by the `assayname` argument.
#' @param assayname Name of expression assay data contained in provided sce
#' object and/or used in plot title.
#' @param group.variable Name of sce colData variable containing group labels.
#' @param group.vector Vector of group labels of length equal to number of 
#' columns in mexpr or sce. This is only used if `mexpr` was provided.
#' @param verbose Whether to show verbose status messages.
#' @returns Table dfstat containing statistics (e.g. "mean" and "var" columns)
#' to use for dispersion analyses.
#' @seealso plot_dispersion, sce_dispersion
#' @examples 
#' sce.exe <- random_sce() # get a random sce
#' dfstat <- get_dispersion_data(sce = sce.exe, group.variable = "celltype")
#' 
#' @export
get_dispersion_data <- function(dfstat = NULL, mexpr = NULL, sce = NULL,
                                group.variable = NULL, group.vector = NULL,
                                assayname = "counts", verbose = FALSE){
  # parse provided expression data
  if(is(dfstat, "NULL")){
    if(is(mexpr, "NULL")){
      if(is(sce, "NULL")){
        stop("Error, provide either sce or mexpr.")
      } else{
        if(assayname %in% names(assays(sce))){
          assays(sce)[[assayname]] <- as.matrix(assays(sce)[[assayname]])
        } else{
          stop("Error, assayname not in sce assays.")
        }
        if(verbose){message("Getting dfstat from sce object...")}
        dfstat <- sce_groupstat(sce = sce, group.variable = group.variable, 
                                assayname = assayname, 
                                groupstat = c("mean", "var"), 
                                summarytype = "rowData", return.tall = TRUE,
                                verbose = verbose)
      }
    } else{
      if(verbose){message("Getting dfstat from mexpr object...")}
      if(is(group.vector, "NULL")){
        dfstat <- data.frame(var = rowVars(mexpr), mean = rowMeans(mexpr),
                             marker = rownames(mexpr))
      } else{
        if(!length(group.vector)==ncol(mexpr)){
          stop("Error, length of group.vector should equal ",
               "total columns in mexpr.")}
        ugroupv <- unique(group.vector)
        dfstat <- do.call(rbind, lapply(ugroupv, function(gi){
          mfilt <- group.vector == gi; mef <- mexpr[,mfilt]
          dfsi <- data.frame(var = rowVars(mef), mean = rowMeans(mef))
          dfsi$group = gi; dfsi
        }))
        dfstat$marker <- rep(rownames(mexpr), length(ugroupv))
      }
    }
  }
  if(verbose){message("Checking dfstat...")}
  cond <- "mean" %in% colnames(dfstat) & "var" %in% colnames(dfstat)
  if(!cond){stop("Error, dfstat requires columns for mean and variance.")}
  return(dfstat)
}

#' plot_dispersion
#'
#' @param dfstat Table containing summary statistics for mean ("mean") and 
#' variance ("var").
#' @param highlight.markers Vector of markers to highlight in plots. Ignored if
#' NULL (default).
#' @param hl.color Color of points or smooth for highlighted markers. See details.
#' @param hl.alpha Transparency of highlighted marker points. See details.
#' @param show.marker.labels Whether to show marker labels using 
#' `ggrepel::geom_text_repel()` (see `?ggrepel::geom_text_repel` for details).
#' @param hl.padding Box outline padding amount for when marker labels are 
#' shown with `show.marker.labels` (see `ggrepel::geom_text_repel` for details).
#' @param nrow.facet Number of rows for facet plots. This is only evaluated if 
#' dfstat contains a column called "group" containing the group labels.
#' @param title.str Character string of main plot title.
#' @param xlab.str Character string of the x-axis title.
#' @param ylab.str Character string of the y-axais title.
#' @param ref.linecol Color to use for the reference line.
#' @param show.smooth Whether to show the smooth lines using `geom_smooth()` 
#' with defaults.
#' @param smooth.linecol Color to use for smoothed model line.
#' @param show.points Whether to show the scatter plot points using 
#' `geom_point()`.
#' @param point.alpha Transparency level for scatter plot points.
#' @param verbose Whether to show verbose status messages.
#' @return A ggplot2 dispersion plot object.
#' @details Plot the dispersion (means versus variances) for some summary 
#' statistics as specified in the `dfstat` object. The `dfstat` object should 
#' include columns for "mean", "var", and optionally for "marker" and "group". 
#' it should be tall, meaning that each row corresponds to a marker-group 
#' summary statistic and there will be >= 1 row per marker if there are >= 1 
#' unique group labels.
#' 
#' The dispersion plot is highly customizable. It is shown on either the 
#' identity or the log10 scale (see `axis.scale` argument). One can plot either
#' the model smooth of the data (e.g. when `show.smooth==TRUE`, the default),
#' or the scatter plot points (e.g. when `show.points==TRUE`, the default), or 
#' both. If more than one group is detected, groups are plotted in separate 
#' facets using `facet_wrap()` with the number of facet rows as specified by
#' `nrow.facet` (the default is 1 row).
#' 
#' Markers can be highlighted in the dispersion plots. By default, when points
#' are plotted, separate points are plotted for the markers provided with
#' `highlight.markers` which are also discovered in dfstat. In this case, the 
#' highlighted points are shown with the color specified by `hl.color` and 
#' transparency level from `hl.alpha`. These parameters are also inherited to
#' show the marker labels with `show.marker.labels` is `TRUE`. When 
#' `show.points` is `FALSE` but `show.smooth` is `TRUE`, this will instead show 
#' the smooth of the discovered markers, where the model color is as specified 
#' in `hl.color`.
#'
#' @seealso get_dispersion_data
#' @examples 
#' # get dispersion data for plotting
#' sce.exe <- random_sce() # get a random sce
#' dfstat <- get_dispersion_data(sce = sce.exe, group.variable = "celltype")
#' 
#' # plot dispersion
#' #
#' # 1. plot without highlighting markers
#' ggds1 <- plot_dispersion(dfstat)
#' markerv <- dfstat$marker[seq(10)] # get markers for highlighting
#' #
#' # 2. highlight some markers without labels
#' markerv <- dfstat$marker[seq(10)]
#' ggds2 <- plot_dispersion(dfstat, highlight.markers = markerv, 
#'                         show.marker.labels = FALSE)
#' #
#' # 3. highlight some markers with labels
#' ggds3 <- plot_dispersion(dfstat, highlight.markers = markerv)
#' #
#' # 4. show smooths of all genes and markers, without points
#' ggds4 <- plot_dispersion(dfstat, highlight.markers = markerv, 
#'                          show.points = FALSE)
#' 
#' @export
plot_dispersion <- function(dfstat, highlight.markers = NULL, 
                            hl.color = "black", hl.alpha = 0.5,
                            show.marker.labels = TRUE, hl.padding = 0.5,
                            axis.scale = "log10", nrow.facet = 1,
                            title.str = "Dispersion plot", ref.linecol = "red", 
                            xlab.str = "Mean (log10-scaled)",
                            ylab.str = "Variance (log10-scaled)",
                            show.smooth = TRUE, smooth.linecol = "blue", 
                            show.points = TRUE, point.alpha = 0.5, 
                            verbose = FALSE){
  if(verbose){message("Plotting dispersion...")}
  # parse marker highlight settings
  hlm.cond <- !is(highlight.markers, "NULL") # get marker plot cond
  if(hlm.cond){
    # get markerv for new cond
    if("marker" %in% colnames(dfstat)){
      marker.filt <- dfstat$marker %in% highlight.markers
      markerv <- dfstat[marker.filt,]$marker
    } else{
      marker.filt <- rownames(dfstat) %in% highlight.markers
      markerv <- rownames(dfstat[marker.filt,])
    }
    hlm.cond <- length(markerv) > 0 # overwrite marker plot cond
    if(hlm.cond){
      dfm <- dfstat[marker.filt,]
    } else{
      message("Warning, no provided highligh markers were found in dfstat.")
    }
  }
  # get basic plot
  ggds <- ggplot(dfstat, aes(x = mean, y = var))
  ggds <- ggds + theme_bw() + ggtitle(title.str) + 
    geom_abline(slope = 1, intercept = 0, color = ref.linecol) +
    xlab(xlab.str) + ylab(ylab.str)
  # parse plot content arguments
  if(axis.scale == "log10"){ggds <- ggds + scale_x_log10() + scale_y_log10()}
  if(show.points){
    ggds <- ggds + geom_point(alpha = point.alpha)
    if(hlm.cond){
      ggds <- ggds + geom_point(data = dfm, aes(x = mean, y = var), 
                                color = hl.color, alpha = hl.alpha)
      if(show.marker.labels){
        require(ggrepel)
        ggds <- ggds + geom_text_repel(data = dfm, 
                                       aes(x = mean, y = var, label = marker), 
                                       color = hl.color, alpha = 1,
                                       box.padding = hl.padding)
      }
    }
  }
  if(show.smooth){
    ggds <- ggds + geom_smooth(color = smooth.linecol)
    if(hlm.cond & !show.points){ # show marker smooth if points not plotted
      ggds <- ggds + geom_smooth(data = dfm, aes(x = mean, y = var), 
                                 color = hl.color, alpha = hl.alpha)
    }
  }
  if("group" %in% colnames(dfstat)){
    ggds <- ggds + facet_wrap(~group, nrow = nrow.facet)
  }
  if(verbose){message("Finished dispersion plot(s). Returning...")}
  return(ggds)
}
