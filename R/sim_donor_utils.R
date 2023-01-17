#!/usr/bin/env R

# Utilities to perform simulations and randomizations with one or more donors or
# cell type data sources.
#

#--------------------
# experiment function
#--------------------

#' donor_marker_sfactorsim
#'
#' @param gindexv Vector of type indices for the G markers. See `?random_lgv` 
#' for details.
#' @param ndonor Total number of donors to simulate.
#' @param ktotal Total K types to simulate.
#' @param num.sim Number of simulations for deconvolution experiments.
#' @param sd.offset.pos Poisson dist mean for randomization of offsets for
#' positive marker signals.
#' @param sd.offset.neg Poisson dist mean for randomization of offsets for
#' negative marker signals.
#' @param lpv List of length num.sim containing true proportions for each 
#' simulated type. Automatically generated if not provided.
#' @param lsv List of length num.sim containing the cell sizes for each type. 
#' Automatically generated with `seed.num` argument if not provided.
#' @param run.decon Whether to run deconvolution experiments, returning 
#' predictions. Automatically generates values for lpv, lsv if none passed.
#' @param seed.num Seed value for random sizes in lsv, in case lsv is NULL.
#' @param ... Arguments passed to functions `rand_donor_marker_table()`.
#' @returns Results of a donor marker experiment, including randomized marker
#' signal table, results of PCA on marker table, and results and plots of 
#' deconvolution predictions.
#' @examples 
#' donor_marker_experiment()
#' @export
donor_marker_sfactorsim <- function(gindexv = c(1, 2), ndonor = 2, ktotal = 2, 
                                    num.sim = 1, sd.offset.pos = 5, 
                                    sd.offset.neg = 5, lpv = NULL, lsv = NULL, 
                                    run.decon = TRUE, seed.num = 0, 
                                    plot.title.append = NULL,
                                    verbose = FALSE, ...){
  if(verbose){message("Getting random marker table...")}
  dt <- rand_donor_marker_table(ndonor = ndonor, gindexv = gindexv, 
                                ktotal = ktotal, 
                                sd.offset.pos = sd.offset.pos,
                                sd.offset.neg = sd.offset.neg, 
                                ...)
  lpca.markers <- pcaplots_donor(dt = dt, title.append = plot.title.append)
  lr <- list(marker.table = dt, lpca.markers = lpca.markers)
  # manage deconvolution experiments
  if(run.decon){
    if(verbose){message("Running deconvolution experiment...")}
    if(is(lpv, "NULL")){lpv <- make_lpv(ktotal)[seq(num.sim)]}
    if(is(lsv, "NULL")){
      set.seed(seed.num)
      sizev <- sample(100, ktotal)
      lsv <- lapply(seq(num.sim), function(ii){sizev})
    }
    # run decon
    cndv <- colnames(dt)[grepl("donor", colnames(dt))]
    ld <- lapply(cndv, function(donori){
      cnvf <- c(donori, "type"); dtf <- dt[,cnvf]
      lgvi <- lapply(seq(ktotal), function(jj){
        dtf[dtf[,2]==paste0("type", jj),1]
      })
      # rep up to num sim
      lgv.in <- lapply(seq(num.sim), function(ii){lgvi})
      decon_analysis(lgv = lgv.in, lpv = lpv, lsv = lsv)
    })
    names(ld) <- cndv; lr$decon.results<- ld
  }
  return(lr)
}

#' donor_marker_biasexpt
#'
#' Compare predictions with and without donor bias corrections.
#' 
#' @param offsetv Vector of donor offset values. This is used to define the
#' test groups for donor bias experiments. Equal offsets are applied to positive
#' and negative marker signal distributions.
#' @param P Vector of true proportions.
#' @param donor.adj.method Method to adjust for donor bias.
#' @param gindexv Vector of type indices for the G markers. See `?random_lgv` 
#' for details.
#' @param ndonor Total number of donors to simulate.
#' @param ktotal Total K types to simulate.
#' @param sd.offset.pos Poisson dist mean for randomization of offsets for
#' positive marker signals.
#' @param sd.offset.neg Poisson dist mean for randomization of offsets for
#' negative marker signals.
#' @param lpv List of length num.sim containing true proportions for each 
#' simulated type. Automatically generated if not provided.
#' @param run.decon Whether to run deconvolution experiments, returning 
#' predictions. Automatically generates values for lpv, lsv if none passed.
#' @param seed.num Seed value for random sizes in lsv, in case lsv is NULL.
#' @param ... Arguments passed to functions `rand_donor_marker_table()`.
#' @returns List of experiment results and experiment objects.
#' @export
donor_marker_biasexpt <- function(offsetv = c(1, 1e3), P = c(0.25, 0.75),
                                  donor.adj.method = 'var_denom',
                                  gindexv = c(1, 2), ndonor = 10, ktotal = 2,
                                  seed.num = 0, verbose = FALSE, ...){
  if(verbose){message("Making new pseudobulk sample...")}
  P <- c(0.25, 0.75)
  Ypb <- t(t(P) %*% t(Z)) 
  if(verbose){message("Getting randomized donor marker data...")}
  # offsetv <- c(1, 1e3)
  ktotal <- length(P)
  ldonordf <- lapply(offsetv, function(offi){
    rand_donor_marker_table(ndonor = ndonor, gindexv = gindexv, ktotal = ktotal,
                            sd.offset.pos = offi, sd.offset.neg = offi)
  })
  names(ldonordf) <- paste0("offset:", offsetv)
  if(verbose){message("Getting type predictions...")}
  cname.data <- "donor.combn.all.mean"
  dfr <- do.call(rbind, lapply(seq(length(ldonordf)), function(ii){
    # simulate donor data
    namei <- names(ldonordf)[ii]; df <- ldonordf[[namei]]
    donor.unadj <- df[,cname.data]
    # donor.adjv <- donoradj(donor.unadj, df, donor.adj.method = "var_denom", ...)
    donor.adjv <- donoradj(donor.unadj, df, donor.adj.method = donor.adj.method)
    # get marker tables
    Zunadj <- matrix(donor.unadj, ncol = ktotal)
    Zadj <- matrix(donor.adjv, ncol = ktotal)
    # get predictions
    punadj <- predtype(Z = Zunadj, Y = Ypb, strict_method = "nnls",
                       proportions = TRUE, verbose = TRUE)
    padj <- predtype(Z = Zadj, Y = Ypb, strict_method = "nnls",
                     proportions = TRUE, verbose = TRUE)
    # append results
    data.frame(prop.type = c(rep("punadj", 2), rep("padj", 2)),
               prop.pred = c(punadj, padj), prop.true = c(rep(P, 2)),
               type.index = c(rep(seq(ktotal), 2)),
               offset = rep(gsub(".*:", "", namei), 4))
  }))
  # get return object
  lmd.adj <- list(donor.adj.method = donor.adj.method, ...)
  lmd <- list(offsetv = offsetv, P = P, donor.adj.info = lmd.adj)
  lr <- list(dfres = dfr, ldonordf, Ypb = Ypb, metadata = lmd)
  return(lr)
}

#' donoradj
#'
#' Apply some specified bias adjustment to a vector of marker data.
#'
#' @param donorv Vector of markrer signals (e.g. for a donor, for some summaries 
#' across donors, etc.).
#' @param donordf A data.frame containing the donor information used for bias
#' corrections. Should contain donor-specific marker info identifiable by 
#' column names of the format: "donor"+"[0-9]".
#' @param method Method to correct for bias. Supports "var_denom", "sd_denom"
#' @param denom_offset Denominator offset, usually very small, for methods which 
#' divide by some quantity.
#' @param bounds_thresh Threshold for denominator offset. Absolute denominator
#' values above this threshold are set to equal this value.
#' @param verbose Whether to show verbose status messages.
#' @returns Vector of same length as donorv, containing the bias-adjusted 
#' values.
#' @export
donoradj <- function(donorv, donordf, method = "var_denom", denom_offset = 1e-3,
                     bounds_thresh = NULL, verbose = FALSE, ...){
  donor.adj <- NA
  if(verbose){message("Getting donor marker data...")}
  cnv <- colnames(donordf)
  donorcol <- cnv[grepl("^donor\\d+$", cnv)]
  if(verbose){message("Found ",length(donorcol),
                      " columns of donor marker data in donordf")}
  dff <- donordf[,donorcol]
  if(verbose){message("Getting adjusted data...")}
  if(grepl(".*_denom$", method)){
    if(verbose){message("Parsing denominator adjustment...")}
    if(method == "var_denom"){
      denomv <- rowVars(as.matrix(dff))
    } else if(method == "sd_denom"){
      denomv <- rowSds(as.matrix(dff))
    } else{
      stop("Error, invalid adjustment method provided.")
    }
    denomv <- denomv + denom_offset
    if(!is(bounds_thresh, "NULL")){
      denomv[which(abs(denomv) >= bounds_thresh)] <- bounds_thresh
    }
    donor.adj <- donorv/denomv
  } else{
    stop("Error, invalid adjustment method provided.")
  }
  return(donor.adj)
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
#' @param sd.offset.pos Poisson dist mean for randomization of offsets for
#' positive marker signals.
#' @param sd.offset.neg Poisson dist mean for randomization of offsets for
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
                                    sd.offset.pos = 10, sd.offset.neg = 2, 
                                    seed.num = 0, verbose = FALSE, ...){
  set.seed(seed.num)
  nmarkers <- length(gindexv)
  # draw random offsets from normal dist
  offposv <- rnorm(n = ndonor, mean = 0, sd = sd.offset.pos)
  offnegv <- rnorm(n = ndonor, mean = 0, sd = sd.offset.neg)
  # get value vectors
  meanv.pos <- offposv + lambda.pos
  meanv.neg <- offnegv + lambda.neg
  # convert negative means
  meanv.pos[meanv.pos < 0] <- -1*meanv.pos[meanv.pos < 0]
  meanv.neg[meanv.neg < 0] <- -1*meanv.neg[meanv.neg < 0]
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
#' @param title.append Optional string to append to plot titles.
#' @param verbose Whether to show verbose status messages.
#' @param ... Additional arguments passed to PCA functions.
#' @returns list of PCA results, plots, and metadata
#' @export
pcaplots_donor <- function(dt, title.append = NULL, verbose = FALSE, ...){
  list(pca.bydonor = pca_bydonor(dt= dt, title.append = title.append, ...), 
       pca.bydonortype = pca_bydonortype(dt = dt, title.append = title.append, ...))
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
#' @param title.append Optional string to append to plot titles.
#' @param verbose Whether to show verbose status messages.
#' @returns lr, results list containing PCA results data, ggplots, and metadata.
#' @export
pca_bydonor <- function(dt, test.md = list(test = "pca", test.type = "by donor"), 
                        title.append = NULL, verbose = FALSE){
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
  title.str.pt <- paste0("PCA across markers; Num. markers = ", num.markers, 
                      "; Num. types = ", num.types)
  if(!is(title.append, "NULL")){
    title.str.pt <- paste0(title.append, title.str.pt)}
  # make scatterplot
  gg.pt <- ggplot(dfp, aes(x = x, y = y, color = donor, 
                           shape = donor.summary)) + 
    geom_point(size = 4, alpha = 0.5) + xlab(colnames(dfp)[1]) +
    ylab(colnames(dfp)[2]) + ggtitle(title.str.pt)
  # get screeplot data
  dfp <- data.frame(pc = colnames(rpca$x), sd = rpca$sdev)
  title.str.scree <- paste0("Screeplot; Num. markers = ", ncol(df.pca)) 
  if(!is(title.append, "NULL")){
    title.str.scree <- paste0(title.append, title.str.scree)}
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
#' @param title.append Optional string to append to plot titles.
#' @param verbose Whether to show verbose status messages.
#' @returns lr, results list containing PCA results data, ggplots, and metadata.
#' @export
pca_bydonortype <- function(dt, 
                            test.md = list(test = "pca", 
                                           test.type = "by donor;type"),
                            title.append = NULL, verbose = FALSE){
  require(ggplot2)
  # run pca
  ntype <- length(unique(dt$type))
  cndv <- colnames(dt)[grepl("^donor.*", colnames(dt))]
  ndonorcat <- length(cndv)
  catv <- paste0(rep(cndv, each = ntype), ";", 
                 rep(paste0("type", seq(ntype)), times = ndonorcat))
  df.pca <- do.call(rbind, lapply(catv, function(cati){
    donor.str.filt <- paste0(gsub(";.*", "", cati), "$")
    type.str.filt <- paste0(gsub(".*;", "", cati), "$")
    colfilt <- grepl(donor.str.filt, colnames(dt))
    rowfilt <- grepl(type.str.filt, dt$type)
    dt[rowfilt, colfilt]
  }))
  rownames(df.pca) <- catv; rpca <- prcomp(df.pca)
  # assign pc labels
  percv <- round(100*rpca$sdev/sum(rpca$sdev),0)
  colnames(rpca$x) <- paste0(colnames(rpca$x), " (",percv,"%)")
  # make scatterplot
  dfp <- as.data.frame(rpca$x); dfp$x <- dfp[,1]; dfp$y <- dfp[,2]
  dfp$donor <- gsub(";.*", "", rownames(df.pca))
  dfp$type <- gsub(".*;", "", rownames(df.pca))
  title.str.pt <- paste0("PCA by donor, marker; Num. markers = ", ncol(df.pca))
  if(!is(title.append, "NULL")){
    title.str.pt <- paste0(title.append, title.str.pt)}
  # plot scatterplot, first 2 pc's
  gg.pt <- ggplot(dfp, aes(x = x, y = y, color = donor, shape = type)) + 
    geom_point(size = 4, alpha = 0.5) + ggtitle(title.str.pt) +
    xlab(colnames(dfp)[1]) + ylab(colnames(dfp)[2])
  # get screeplot data
  dfp2 <- data.frame(pc = colnames(rpca$x), sd = rpca$sdev)
  title.str.scree <- paste0("Num. markers = ", ncol(df.pca))
  if(!is(title.append, "NULL")){
    title.str.scree <- paste0(title.append, title.str.scree)}
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
    if(!is(title.append, "NULL")){
      title.str.pairs <- paste0(title.append, title.str.pairs)}
    gg.pairs <- ggpairs(dfp, top = list(continuous="na"), 
                        columns = seq(ncol(df.pca)), 
                        map = ggplot2::aes(color=donor, shape=type),
                        title = title.str.pairs)
    lr$pairs = gg.pairs
  }
  return(lr)
}