#!/usr/bin/env R

# Utilities to perform simulations and randomizations with one or more donors or
# cell type data sources.
#

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
rand_donor_marker_table <- function(ndonor, gindexv = c(1, 2), ktotal = 2, 
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

#------
# plots
#------

#' pcaplots_donor
#'
#' Make plots of two PCAs: (1) by donor, across markers and types; (2) by donor
#' and type, across markers.
#' 
#' @param dt Donor marker signals table.
#'
#'
pcaplots_donor <- function(dt, verbose = FALSE, ...){
  require(ggplot2)
  
}

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
  lr <- list(pca.results = rpca, scatterplot = gg.pt, screeplot = gg.bp, 
             metadata = test.md)
  return(lr)
}

#' pca_bydonortype
#'
#' Get analysis results and plots for PCA by donor and type, across markers.
#' 
#' 
#' 
#'
pca_bydonortype <- function(dt, verbose = FALSE){
  
}