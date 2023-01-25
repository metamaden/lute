#!/usr/bin/env R

# Utilities to perform simulations and randomizations with one or more donors or
# cell type data sources.
#
# Contents:
#
# 1. main experiment functions: main experiment wrapper scripts
#
# 2. experiment utility methods: methods to manage experiments
# 
# 3. donor bias adjustment methods: methods for performing bias adjustments
#
# 4. methods for plotting results: plot results of bias adjustments and pca, 
# returning the plot objects
#


#-----------------------------
# 1. main experiment functions
#-----------------------------
# main experiment wrapper scripts

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
#' @param verbose Whether to show verbose status messages.
#' @param ... Arguments passed to functions `random_donordf()`.
#' @returns Results of a donor marker experiment, including randomized marker
#' signal table, results of PCA on marker table, and results and plots of 
#' deconvolution predictions.
#' @examples 
#' donor_marker_sfactorsim()
#' @export
donor_marker_sfactorsim <- function(gindexv = c(1, 2), ndonor = 2, ktotal = 2, 
                                    num.sim = 1, sd.offset.pos = 5, 
                                    sd.offset.neg = 5, lpv = NULL, lsv = NULL, 
                                    run.decon = TRUE, seed.num = 0, 
                                    plot.title.append = NULL,
                                    verbose = FALSE, ...){
  if(verbose){message("Getting random marker table...")}
  dt <- random_donordf(ndonor = ndonor, gindexv = gindexv, 
                       sd.offset.pos = sd.offset.pos, 
                       sd.offset.neg = sd.offset.neg, ...)
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

#' run_donor_bias_expt
#'
#' Run a deconvolution experiment evaluating the impact of donor bias on 
#' predictions.
#' 
#' @param donordf Table of type `donor.data.frame`. If NULL, attempt to makes a 
#' new table from the provided arguments.
#' @param method Randomization method passed to random_lgv(). Supports either 
#' "nbinom" for negative binomial distribution (a.k.a. gamma poisson 
#' distribution) or "poisson" for poisson distribution.
#' @param lambda.pos The mean or mu value when marker status is positive.
#' @param lambda.neg The mean or mu value when marker status is negative.
#' @param lambda.sdoff.pos Offset SD for lambda when marker status positive.
#' @param lambda.sdoff.neg Offset SD for lambda when marker status negative.
#' @param gamma.pos Magnitude of gamma dispersion value when marker status is 
#' positive.
#' @param gamma.neg Magnitude of gamma dispersion value when marker status is
#' negative.
#' @param P Vector of true proportions, where each entry corresponds to a type.
#' @param donor.adj.method Method to adjust for donor bias. Can be either 
#' "limma", "var_denom", "sd_denom", "combat", or NULL. If NULL, skip this step.
#' @param plot.biasadj Whether to make scatterplot of donor summary signals
#' before and after bias adjustment.
#' @param plot.pca Whether to include PCA results plots using simulated donor
#' signals data.frame.
#' @param cname.donorsummary Name of column containing the donor summary data
#' with which to perform experiment (default "donor.combn.all.mean" for cross-
#' donor means).
#' @param gindexv Vector of type indices for the G markers. See `?random_lgv` 
#' for details.
#' @param ndonor Total number of donors to simulate.
#' @param seed.num Seed value for random sizes in lsv, in case lsv is NULL.
#' @param verbose Whether to show verbose status messages.
#' @param ... Arguments passed to function `donoradj()`.
#' @details Parses various methods for performing the donor bias adjustment, 
#' returning predicted type proportions for adjusted and unadjusted signature 
#' matrices. Optionally, can also run PCA and return scatterplots of PCA 
#' results, deconvolution predictions, and bias.
#' 
#' Each individual donor bias experiment consists of the following steps:
#' 
#' * 1. Make the reference pseudobulked sample `Ypb` from a single simulated 
#'      donor.
#' * 2. Simulate donor signals, use the cross-donor means to make the first
#'      signature matrix `Z1`, and use this to get the first predictions `P1`.
#' * 3. Adjust the donor signals with a covariate for type, use the means 
#'      to make the second signature matrix `Z2`, and use this to get the second
#'      set of predictions `P2`.
#' * 4. Return table of type `donor.data.frame`, experiment results df, and 
#'      plots.
#' 
#' @returns List of experiment results and experiment objects.
#' @examples 
#' lb <- donor_marker_biasexpt()
#' @seealso biasexpt
#' @export
run_donor_bias_expt <- function(donordf = NULL, method = "nbinom", 
                                  lambda.pos = 20, lambda.neg = 2,
                                  lambda.sdoff.pos = 5, lambda.sdoff.neg = 2,
                                  gamma.pos = 20, gamma.neg = 2, P = c(0.25, 0.75),
                                  donor.adj.method = "combat", plot.biasadj = TRUE,
                                  plot.pca = TRUE,
                                  cname.donorsummary = "donor.combn.all.mean",
                                  gindexv = c(1, 2), ndonor = 10, seed.num = 0,
                                  verbose = TRUE, ...){
  set.seed(seed.num); lr <- list()
  if(verbose){message("Making pseudobulk sample from types matrix...")}
  df <- random_donordf(ndonor = 1, gindexv = gindexv,
                       lambda.sdoff.pos = 0, lambda.sdoff.neg = 0)
  ktotal <- length(P)
  Z <- matrix(df[,"donor1"], ncol = ktotal)
  Ypb <- ypb_fromtypes(Z = Z, P = P)
  
  if(is(donordf, "NULL")){
    if(verbose){message("Getting randomized donor marker data...")}  
    donordf <- random_donordf(ndonor = ndonor, method = method,
                              gindexv = gindexv, 
                              lambda.sdoff.pos = lambda.sdoff.pos,
                              lambda.sdoff.neg = lambda.sdoff.neg,
                              gamma.pos = gamma.pos, gamma.neg = gamma.neg,
                              seed.num = seed.num)
  } else{
    if(verbose){message("Checking if provided donordf passes checks...")}
    if(!check_donordf(donordf)){
      stop("Error, provided donordf is invalid. ",
           "Check that it contains all required columns.")}
  }
  
  if(verbose){message("Getting type predictions...")}
  type.indexv <- seq(ktotal)
  donor.unadj <- donordf[,cname.donorsummary] # get donor summary datas
  lbias <- biasexpt(df = donordf, Ypb = Ypb, P = P, donor.unadj = donor.unadj,
                    donor.adj.method = donor.adj.method,
                    plot.biasadj = plot.biasadj, verbose = verbose, ...)
  
  # get return object
  lr[["dfres"]] <- lbias$dfres # get results df
  lr[["donordf"]] <- donordf
  lr[["Ypb"]] <- Ypb
  lr[["adj.method"]] <- method
  # make new plots
  if(plot.pca){lr[["lpca"]] <- pcaplots_donor(dt = donordf)}
  if(plot.biasadj){lr[["ggpt.biasadj"]] <- lbias$ggpt.biasadj}
  
  return(lr)
}

#' biasexpt
#'
#' Main function to run a donor bias adjustment experiment.
#' 
#' @param donor.unadj Vector of unadjusted donor signals.
#' @param df Donor data.frame.
#' @param Ypb Pseudobulk sample data.
#' @param P Vector of true proportions, where each entry corresponds to a type.
#' @param donor.adj.method Method to adjust for donor bias. Can be either 
#' "limma", "var_denom", "sd_denom", or NULL. If NULL, skip this step.
#' @param plot.biasadj Whether to make scatterplot of donor summary signals
#' before and after bias adjustment.
#' @param verbose Whether to show verbose status messages.
#' @param ... Additional arguments passed to `donoradj()`.
#' @returns lexpt, experiment results list.
#' @examples 
#' # simulate donor signals data
#' donordf <- random_donordf()
#' 
#' # get ypb
#' P <- c(0.75, 0.25)
#' ktotal <- length(P)
#' Z <- matrix(df[,"donor1"], ncol = ktotal)
#' Ypb <- ypb_fromtypes(Z = Z, P = P)
#' 
#' # get bias expt results
#' li <- biasexpt(df = donordf, Ypb = Ypb)
#' 
#' @seealso run_donor_bias_expt
#' @export
biasexpt <- function(df, Ypb, P, donor.unadj = NULL,donor.adj.method = "combat", 
                     plot.biasadj = TRUE, verbose = FALSE, ...){
  lr <- list() # begin return list
  if(!check_donordf(df)){
    stop("Data.frame is not a valid simulated donor signals data.frame.")}
  if(is(donor.unadj, "NULL")){donor.unadj <- df[,"donor.combn.all.mean"]}
  ktotal <- length(unique(df$type))
  Zunadj <- matrix(donor.unadj, ncol = ktotal)
  punadj <- predtype(Z = Zunadj, Y = Ypb, strict_method = "nnls",
                     proportions = TRUE, verbose = verbose)
  # initial variable defs
  prop.typev <- rep("punadj", ktotal)
  ppredv <- punadj; ptruev <- P
  type.indexv <- seq(ktotal)
  lr[["donor.unadj"]] <- donor.unadj
  if(!is(donor.adj.method, "NULL")){
    donor.adjv <- donoradj(df = df, donor.unadj = donor.unadj,
                           donor.adj.method = donor.adj.method, ...)
    lr[["donor.adj"]] <- donor.adjv
    Zadj <- matrix(donor.adjv, ncol = ktotal)
    padj <- predtype(Z = Zadj, Y = Ypb, strict_method = "nnls",
                     proportions = TRUE, verbose = verbose)
    # append to variable defs
    ptruev <- c(ptruev, P); ppredv <- c(ppredv, padj)
    prop.typev <- c(prop.typev, rep("padj", ktotal))
    type.indexv <- rep(type.indexv, 2)
    # parse plot arg
    if(plot.biasadj){
      dfp <- data.frame(unadj = lr[["donor.unadj"]], adj = lr[["donor.adj"]],
                        marker = df$marker, type = df$type)
      lr[["ggpt.biasadj"]] <- ggpt_donorbias(dfp, method.str = donor.adj.method)
    }
  }
  # append results
  biasv <- ptruev - ppredv
  lr[["dfres"]] <- data.frame(prop.type = prop.typev, prop.pred = ppredv,
                            prop.true = ptruev, bias = biasv,
                            type.index = type.indexv)
  return(lr)
}

#------------------------------
# 2. experiment utility methods
#------------------------------
# methods to manage experiments

#' ypb_fromtypes
#'
#' Makes a pseuobulked convoluted signals sample from types data. Uses the 
#' formula:
#' 
#' Y = P * Z
#' 
#' @param Z Signature matrix containing pure type signals (rows = markers, 
#' columns = types).
#' @param P Vector of true proportions.
#' @returns Ypb, a new pseudobulked sample of convoluted signals.
#' @examples 
#' P <- c(0.75, 0.25)
#' Z <- matrix(c(0, 1, 0, 1), nrow = 2)
#' Ypb <- ypb_fromtypes(Z, P)
#' @export
ypb_fromtypes <- function(Z, P){
  t(t(P) %*% t(Z))
}

#' pcaplots_donor
#'
#' Make plots of two PCAs: (1) by donor, across markers and types; (2) by donor
#' and type, across markers.
#' 
#' @param donordf Table of class `donor.data.frame`.
#' @param title.append Optional string to append to plot titles.
#' @param verbose Whether to show verbose status messages.
#' @param ... Additional arguments passed to PCA functions.
#' @returns list of PCA results, plots, and metadata
#' @export
pcaplots_donor <- function(donordf, title.append = NULL, verbose = FALSE, ...){
  list(pca.bydonor = pca_bydonor(dt = donordf, title.append = title.append, 
                                 verbose = verbose, ...), 
       pca.bydonortype = pca_bydonortype(dt = donordf, 
                                         title.append = title.append, 
                                         verbose = verbose, ...))
}

#---------------------------------
# 3. donor bias adjustment methods
#---------------------------------
# methods for performing bias adjustments

#' donoradj
#'
#' Apply some specified bias adjustment to a vector of marker data.
#'
#' @param df A data.frame containing the donor information used for bias
#' corrections. Should contain donor-specific marker info identifiable by 
#' @param donorv Vector of markrer signals (e.g. for a donor, for some summaries 
#' across donors, etc.).
#' column names of the format: "donor"+"[0-9]".
#' @param method Method to correct for bias. Supports "var_denom", "sd_denom",
#' "combat".
#' @param denom_offset Denominator offset, usually very small, for methods which 
#' divide by some quantity.
#' @param bounds_thresh Threshold for denominator offset. Absolute denominator
#' values above this threshold are set to equal this value.
#' @param verbose Whether to show verbose status messages.
#' @returns Vector of same length as donorv, containing the bias-adjusted 
#' values.
#' @examples 
#' donordf <- random_donordf()
#' donorv <- donordf$donor.combn.all.mean
#' donoradj(donorv, donordf, method = "combat")
#' 
#' @seealso biasexpt, 
#' @export
donoradj <- function(df, donorv = NULL, method = "combat", denom_offset = 1e-3,
                     bounds_thresh = NULL, verbose = FALSE, ...){
  donor.adj <- NA
  if(is(donorv, "NULL")){donorv <- df[,"donor.combn.all.mean"]}
  if(verbose){message("Getting donor marker data...")}
  cnv <- colnames(df); donorcol <- cnv[grepl("^donor\\d+$", cnv)]
  if(verbose){message("Found ",length(donorcol),
                      " columns of donor marker data in donordf")}
  dff <- df[,donorcol]
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
  } else if(method == "combat"){
    # return the means of adjusted expression
    donor.adj <- donoradj_combat(df = df, return.type = "donor.adj")
  } else{stop("Error, invalid adjustment method provided.")}
  return(donor.adj)
}

#' donoradj_combat
#'
#' Use method sva::ComBat() to adjust for donor bias.
#' 
#' @param df Data.frame containing marker signals and pheno data.
#' @param return.type Type of data to return. Accepts either "donor.adj" for a 
#' vector of adjusted donor signals, "mexpr" for a matrix of adjusted 
#' signals, or "donordf" for a donordf-type data.frame.
#' @param verbose Whether to show verbose status messages.
#' @returns madj, matrix of adjusted marker signals.
#' @examples 
#' df <- random_donordf()
#' donoradj_combat(df)
#' @export
donoradj_combat <- function(df, return.type = "donor.adj", verbose = FALSE){
  require(sva)
  if(verbose){message("Getting combat variables...")}
  mexpr <- mexpr_from_donordf(df) # make expr matrix
  cnv <- colnames(mexpr)
  pheno <- data.frame(donor = gsub(";.*", "", cnv),
                      type = gsub(".*;", "", cnv))
  mod <- model.matrix(~type, data = pheno)
  batch <- pheno$donor
  if(verbose){message("Running combat...")}
  if(verbose){
    madj <- ComBat(dat = mexpr, batch = batch, mod = mod,
                   par.prior = TRUE, prior.plots = FALSE)
  } else{
    madj <- suppressMessages(
      ComBat(dat = mexpr, batch = batch, mod = mod,
             par.prior = TRUE, prior.plots = FALSE)
    )
  }
  if(verbose){message("Parsing return.type...")}
  if(return.type == "donor.adj"){
    ltype <- list(type = gsub(".*;", "", colnames(madj)))
    dfa <- aggregate(t(madj), by = ltype, FUN = "mean")
    dfa <- t(dfa[,2:ncol(dfa)])
    donor.adj <- unlist(lapply(seq(nrow(dfa)), function(ri){dfa[ri,]}))
    return(donor.adj)
  } else if(return.type == "mexpr"){
    return(madj)
  } else if(return.type == "donordf"){
    return(donordf_from_mexpr(mexpr = madj, verbose = verbose))
  } else{
    stop("Error, invalid return.type specified. ",
         "Should be either 'donor.adj' or 'mexpr'.")
  }
  return(NULL)
}

#--------------------------------
# 4. methods for plotting results
#--------------------------------
# plot results of bias adjustments and pca, returning the plot objects

#' ggpt_donorbias
#'
#' Make a scatterplot of the unadjusted and adjusted donor signals.
#' 
#' @param dfp Plot data.frame, containing columns titled: "unadj" (unadjusted 
#' donor signals), "adj" (adjusted donor signals), "marker" (marker labels), and
#' "type" (type labels).
#' @param method.str Character string describing adjustment method used.
#' @returns A ggplot scatterplot object.
#' @export
ggpt_donorbias <- function(dfp, method.str = "combat"){
  require(ggplot2)
  plot.titlestr <- paste0("Donor signals\n")
  plot.titlestr <- paste0(plot.titlestr, "Adj. method: ", method.str)
  ggplot(dfp, aes(x = unadj, y = adj, color = marker, shape = type)) + 
    theme_bw() + geom_point(alpha = 0.5, size = 4) + 
    geom_abline(slope = 1, intercept = 0, color = "black") +
    ggtitle(plot.titlestr)
}

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
  df.pca <- t(dt[,grepl("^donor.*", colnames(dt))])
  rpca <- prcomp(df.pca)
  # assign pc labels
  percv <- round(100*rpca$sdev/sum(rpca$sdev),0)
  colnames(rpca$x) <- paste0(colnames(rpca$x), " (",percv,"%)")
  # get scatterplot data
  num.markers <- length(unique(dt$marker))
  num.types <- length(unique(dt$type))
  dfp <- as.data.frame(rpca$x)
  dfp$x <- dfp[,1]; dfp$y <- dfp[,2]
  dfp$donor <- rownames(df.pca)
  dfp$donor.summary <- ifelse(grepl(".*\\.mean$|.*\\.median$", dfp$donor), 
                              TRUE, FALSE)
  title.str.pt <- paste0("PCA across markers; Num. markers = ", num.markers, 
                      "; Num. types = ", num.types)
  if(!is(title.append, "NULL")){
    title.str.pt <- paste0(title.append, title.str.pt)}
  # make scatterplot
  gg.pt <- ggplot(dfp, aes(x = x, y = y, color = donor, 
                           shape = donor.summary)) + theme_bw() +
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
    theme_bw() + geom_point(size = 4, alpha = 0.5) + ggtitle(title.str.pt) +
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