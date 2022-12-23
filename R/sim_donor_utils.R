#!/usr/bin/env R

# Utilities to perform simulations and randomizations with one or more donors or
# cell type data sources.
#

#--------------------
# experiment function
#--------------------

#' donor_marker_experiment
#'
#' @param ... Arguments passed to functions rand_donor_marker_table() and
#' run_decon_expr().
#' @returns Results of a donor marker experiment, including randomized marker
#' signal table, results of PCA on marker table, and results and plots of 
#' deconvolution predictions.
#'
donor_marker_experiment <- function(...){
  
  rand_donor_marker_table()
  pcaplots_donor()
}

#---------------------
# experiment utilities
#---------------------

#' rand_donor_marker_table
#' 
#' Get a flat table of random donor marker signals by types.
#' 
#' @param ndonor Number of donors to simulate.
#' @param gindexv Vector of marker indices. Index values correspond to the k types,
#' and each index position represents a marker (e.g. c(1,2,2) means two markers 
#' for the second type, etc.).
#' @param ktotal Total types to simulate.
#' @param lambda.pos Value of lambda (Poisson dist. mean) for "positive" marker 
#' status (e.g. mean of dist. for k when marker is positive for k, negative for 
#' not-k).
#' @param lambda.neg Value of lambda (Poisson dist. mean) for "negative" marker 
#' status (e.g. mean of dist. for k when marker is positive for not-k, negative 
#' for k).
#' @param mean.offset.pos Poisson dist mean for randomization of offsets for
#' positive marker signals.
#' @param mean.offset.neg Poisson dist mean for randomization of offsets for
#' negative marker signals.
#' @param seed.num Token to set the random seed.
#' @param vebose Whether to return verbose status messages.
#' @param ... Additional parameters passed to `random_lgv()` to get marker 
#' signals.
#' @returns return 
#' @examples
#' 
#' get_donor_marker_flattable(ndonor = 2, gindexv = c(1,2))
#' 
#' get_donor_marker_flattable(ndonor = 10, gindexv = c(1,1,2))
#' 
#' get_donor_marker_flattable(ndonor = 10, gindexv = c(rep(1, 10), rep(2, 20)))
#' 
#' @seealso random_lgv
#' @export
rand_donor_marker_table <- function(ndonor = 2, gindexv = c(1, 2), ktotal = 2, 
                                    lambda.pos = 20, lambda.neg = 2,
                                    mean.offset.pos = 10, mean.offset.neg = 2, 
                                    seed.num = 0, verbose = FALSE, ...){
  set.seed(seed.num)
  nmarkers <- length(gindexv)
  # draw random offsets from normal dist
  offposv <- rnorm(n = ndonor, mean = mean.offset.pos)
  offnegv <- rnorm(n = ndonor, mean = mean.offset.neg)
  # get value vectors
  meanv.pos <- offposv + lambda.pos
  meanv.neg <- offnegv + lambda.neg
  # convert negative means
  meanv.pos[meanv.pos < 0] <- -1*meanv.pos
  meanv.neg[meanv.neg < 0] <- -1*meanv.neg
  # get matrix of markers (rows) by donors (cols)
  md <- do.call(cbind, lapply(seq(ndonor), function(ii){
    unlist(random_lgv(gindexv, num.iter = 1, ktotal = ktotal,
                      lambda.pos = meanv.pos[ii],
                      lambda.neg = meanv.neg[ii], 
                      ...))
  }))
  md <- as.data.frame(md)
  colnames(md) <- paste0("donor", seq(ndonor))
  if(ndonor > 1){
    if(verbose){message("Getting donor summary columns...")}
    which.cnv.donor <- which(grepl("donor", colnames(md)))
    md$donor.combn.all.mean <- apply(md[,which.cnv.donor], 1, mean)
    md$donor.combn.all.median <- apply(md[,which.cnv.donor], 1, median)
  }
  md$type <- paste0("type", rep(seq(ktotal), each = nmarkers))
  md$marker <- paste0("marker", rep(seq(nmarkers), times = ktotal))
  md$marker.type <- paste0("type", gindexv)
  return(md)
}

#' pcaplots_donor
#'
#' Make plots of two PCAs: (1) by donor, across markers and types; (2) by donor
#' and type, across markers.
#' 
#' @param dt Donor marker signals table.
#' @param verbose Whether to show verbose status messages.
#' @param ... Additional arguments passed to PCA functions.
#' @returns list of PCA results, plots, and metadata
#' @export
pcaplots_donor <- function(dt, verbose = FALSE, ...){
  list(pca.bydonor = pca_bydonor(dt, ...), 
       pca.bydonortype = pca_bydonortype(dt, ...),
       metadata = test.md)
}

#------
# plots
#------

#' pca_bydonor
#'
#' Get analysis results and plots for PCA by donor, across markers and types.
#'
#' @param dt Data table containing the donors (rows) by marker signals (columns)
#' for each marker and type.
#' @param test.md Metadata info to be returned with test results
#' @param verbose Whether to show verbose status messages.
#' @returns lr, results list containing PCA results data, ggplots, and metadata.
#' @export
pca_bydonor <- function(dt, test.md = list(test = "pca", test.type = "by donor"), 
                        verbose = FALSE){
  require(ggplot2)
  # run pca
  df.pca <- t(dt[,grepl("donor", colnames(dt))]); rpca <- prcomp(df.pca)
  # assign pc labels
  percv <- round(100*rpca$sdev/sum(rpca$sdev),0)
  colnames(rpca$x) <- paste0(colnames(rpca$x), " (",percv,"%)")
  # get scatterplot data
  num.markers <- length(unique(dt$marker))
  num.types <- length(unique(dt$type))
  dfp <- as.data.frame(rpca$x)
  dfp$x <- dfp[,1]; dfp$y <- dfp[,2]
  dfp$donor <- rownames(df.pca)
  dfp$donor.summary <- ifelse(grepl("mean|median", dfp$donor), TRUE, FALSE)
  title.str <- paste0("PCA across markers; Num. markers = ", num.markers, 
                      "; Num. types = ", num.types)
  # make scatterplot
  gg.pt <- ggplot(dfp, aes(x = x, y = y, color = donor, 
                           shape = donor.summary)) + 
    geom_point(size = 4, alpha = 0.5) + xlab(colnames(dfp)[1]) +
    ylab(colnames(dfp)[2]) + ggtitle(title.str)
  # get screeplot data
  dfp <- data.frame(pc = colnames(rpca$x), sd = rpca$sdev)
  title.str.scree <- paste0("Screeplot; Num. markers = ", ncol(df.pca)) 
  # make screeplot
  gg.bp <- ggplot(dfp, aes(x = pc, y = sd)) + geom_bar(stat="identity") + 
    theme_bw() + ggtitle(title.str.scree) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  # make return list
  lr <- list(pca.results = rpca, scatterplot.pc1.pc2 = gg.pt, 
             screeplot = gg.bp, metadata = test.md)
  return(lr)
}

#' pca_bydonortype
#'
#' Get analysis results and plots for PCA by donor and type, across markers.
#' 
#' @param dt Data table containing the donors (rows) by marker signals (columns)
#' for each marker and type.
#' @param test.md Metadata info to be returned with test results
#' @param verbose Whether to show verbose status messages.
#' @returns lr, results list containing PCA results data, ggplots, and metadata.
#' @export
pca_bydonortype <- function(dt, 
                            test.md = list(test = "pca", 
                                           test.type = "by donor;type"), 
                            verbose = FALSE){
  require(ggplot2)
  # run pca
  ntype <- length(unique(dt$type))
  cndv <- colnames(dt)[grepl("^donor.*", colnames(dt))]
  ndonorcat <- length(cndv)
  catv <- paste0(rep(cndv, each = ntype), ";", 
                 rep(paste0("type", seq(ntype)), times = ndonorcat))
  df.pca <- do.call(rbind, lapply(catv, function(cati){
    colfilt <- grepl(gsub(";.*", "", cati), colnames(dt))
    rowfilt <- grepl(gsub(".*;", "", cati), dt$type)
    dtf <- dt[rowfilt, colfilt]
  }))
  rownames(df.pca) <- catv; rpca <- prcomp(df.pca)
  # assign pc labels
  percv <- round(100*rpca$sdev/sum(rpca$sdev),0)
  colnames(rpca$x) <- paste0(colnames(rpca$x), " (",percv,"%)")
  # make scatterplot
  dfp <- as.data.frame(rpca$x); dfp$x <- dfp[,1]; dfp$y <- dfp[,2]
  dfp$donor <- gsub(";.*", "", rownames(df.pca))
  dfp$type <- gsub(".*;", "", rownames(df.pca))
  title.str <- paste0("PCA by donor, marker; Num. markers = ", ncol(df.pca))
  # plot scatterplot, first 2 pc's
  gg.pt <- ggplot(dfp, aes(x = x, y = y, color = donor, shape = type)) + 
    geom_point(size = 4, alpha = 0.5) + ggtitle(title.str) +
    xlab(colnames(dfp)[1]) + ylab(colnames(dfp)[2])
  # get screeplot data
  dfp2 <- data.frame(pc = colnames(rpca$x), sd = rpca$sdev)
  title.str.scree <- paste0("Num. markers = ", ncol(df.pca))
  # make screeplot
  gg.bp <- ggplot(dfp2, aes(x = pc, y = sd)) + geom_bar(stat="identity") + 
    theme_bw() + ggtitle(title.str.scree) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  # get return object
  lr <- list(pca.results = rpca, scatterplot.pc1.pc2 = gg.pt, 
             screeplot = gg.bp, metadata = test.md)
  if(ncol(rpca$x) > 2){ # plot ggpairs
    if(verbose){message("Making ggpairs scatterplots for >2 PCs...")}
    title.str.pairs <- paste0("Num. markers = ", ncol(df.pca))
    gg.pairs <- ggpairs(dfp, top = list(continuous="na"), 
                        columns = seq(ncol(df.pca)), 
                        map = ggplot2::aes(color=donor, shape=type),
                        title = title.str.pairs)
    lr$pairs = gg.pairs
  }
  return(lr)
}