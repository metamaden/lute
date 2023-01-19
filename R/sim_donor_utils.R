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
#' @param ... Arguments passed to functions `rand_donor_marker_table()`.
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
  dt <- rand_donor_marker_table(ndonor = ndonor, gindexv = gindexv, 
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
#' Compare predictions with and without donor bias corrections. Parses various
#' methods for performing the donor bias adjustment, returning predicted 
#' type proportions for adjusted and unadjusted signature matrices.
#' 
#' @param offsetv Vector of donor offset values. This is used to define the
#' test groups for donor bias experiments. Equal offsets are applied to positive
#' and negative marker signal distributions.
#' @param P Vector of true proportions, where each entry corresponds to a type.
#' @param donor.adj.method Method to adjust for donor bias. Can be either 
#' "limma", "var_denom", "sd_denom", "combat", or NULL. If NULL, skip this step.
#' @param plot.biasadj Whether to make scatterplot of donor summary signals
#' before and after bias adjustment.
#' @param plot.pca Whether to include PCA results plots using simulated donor
#' signals data.frame.
#' @param cname.donorsummary Name of column containing the donor summary data
#' with which to perform experiment.
#' @param gindexv Vector of type indices for the G markers. See `?random_lgv` 
#' for details.
#' @param ndonor Total number of donors to simulate.
#' @param seed.num Seed value for random sizes in lsv, in case lsv is NULL.
#' @param verbose Whether to show verbose status messages.
#' @param ... Arguments passed to function `donoradj()`.
#' @returns List of experiment results and experiment objects.
#' @examples 
#' lb <- donor_marker_biasexpt()
#' @export
donor_marker_biasexpt <- function(offsetv = c(1, 10), P = c(0.25, 0.75),
                                  donor.adj.method = NULL, 
                                  plot.biasadj = TRUE, plot.pca = TRUE,
                                  cname.donorsummary = "donor.combn.all.mean",
                                  gindexv = c(1, 2), ndonor = 10,
                                  seed.num = 0, verbose = FALSE, ...){
  set.seed(seed.num); lr <- list()
  if(verbose){message("Making pseudobulk sample from types matrix...")}
  df <- rand_donor_marker_table(ndonor = 1, gindexv = gindexv,
                                sd.offset.pos = 0, sd.offset.neg = 0)
  ktotal <- length(P)
  Z <- matrix(df[,"donor1"], ncol = ktotal)
  Ypb <- ypb_fromtypes(Z = Z, P = P)
  if(verbose){message("Getting randomized donor marker data...")}
  ldonordf <- lapply(offsetv, function(offi){
    rand_donor_marker_table(ndonor = ndonor, gindexv = gindexv,
                            sd.offset.pos = offi, sd.offset.neg = offi,
                            seed.num = seed.num)
  })
  names(ldonordf) <- paste0("offset:", offsetv)
  if(verbose){message("Getting type predictions...")}
  type.indexv <- seq(ktotal)
  lexpt <- lapply(seq(length(ldonordf)), function(ii){
    namei <- names(ldonordf)[ii]
    df <- ldonordf[[namei]] # get full donordf
    offsetv <- rep(gsub(".*:", "", namei), ktotal)
    donor.unadj <- df[,cname.donorsummary] # get donor summary datas
    li <- biasexpt(df = df, Ypb = Ypb, P = P, donor.unadj = donor.unadj,
                    donor.adj.method = donor.adj.method,
                    plot.biasadj = plot.biasadj,
                    verbose = verbose, ...)
    # append offset values
    li$dfi$offset <- rep(offsetv, nrow(li$dfi)/length(offsetv))
    return(li)
  })
  names(lexpt) <- names(ldonordf)
  # get results df
  dfres <- do.call(rbind, lapply(lexpt, function(ii){ii$dfi}))
  ldonorv <- lapply(lexpt, function(ii){ii[c("donor.unadj", "donor.adj")]})
  names(ldonorv) <- names(ldonordf)
  # get return object
  lmd.adj <- list(donor.adj.method = donor.adj.method, ...)
  lmd <- list(offsetv = offsetv, P = P, donor.adj.info = lmd.adj)
  lr[["dfres"]] <- dfres
  lr[["ldonorv"]] <- ldonorv
  lr[["ldonordf"]] <- ldonordf
  lr[["Ypb"]] <- Ypb
  lr[["metadata"]] <- lmd
  # get plot objects
  if(plot.pca){
    lpca <- lapply(ldonordf, function(dfi){
      pcaplots_donor(dt = dfi, title.append = NULL)
    })
    names(lpca) <- names(ldonordf); lr[["lpca"]] <- lpca
  }
  if(plot.biasadj){
    lpt <- lapply(lexpt, function(ii){ii$ggpt.biasadj})
    names(lpt) <- names(lexpt)
    lr[["ggpt.biasadj"]] <- lpt
  }
  return(lr)
}

#' biasexpt
#'
#' Run a single donor bias experiment.
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
#' donordf <- rand_donor_marker_table()
#' # get ypb
#' P <- c(0.75, 0.25)
#' ktotal <- length(P)
#' Z <- matrix(df[,"donor1"], ncol = ktotal)
#' Ypb <- ypb_fromtypes(Z = Z, P = P)
#' # get bias expt results
#' li <- biasexpt(df = donordf, Ypb = Ypb)
#' @export
biasexpt <- function(df, Ypb, P, donor.unadj = NULL, 
                     donor.adj.method = "combat", 
                     plot.biasadj = TRUE, verbose = FALSE, ...){
  lr <- list() # begin return list
  if(!check_donordf(df)){
    stop("Data.frame is not a valid simulated donor signals data.frame.")}
  if(is(donor.unadj, "NULL")){donor.unadj <- df[,"donor.combn.all.mean"]}
  ktotal <- length(unique(df$type))
  Zunadj <- matrix(donor.unadj, ncol = ktotal)
  punadj <- predtype(Z = Zunadj, Y = Ypb, strict_method = "nnls",
                     proportions = TRUE, verbose = TRUE)
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
                     proportions = TRUE, verbose = TRUE)
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
  lr[["dfi"]] <- data.frame(prop.type = prop.typev, prop.pred = ppredv,
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

#' check_donordf
#'
#' Checks whether data.frame is a valid simulated donor signals data.frame by
#' evaluating the column names.
#' 
#' @param df Data.frame to check.
#' @returns boolean, TRUE if df is valid, FALSE otherwise
#' @export
check_donordf <- function(df){
  cnv <- colnames(df); lcond <- list()
  lcond[["cond.donorcol"]] <- grepl("^donor\\d", cnv)
  lcond[["cond.typecol"]] <- grepl("^type$", cnv)
  lcond[["cond.markercol"]] <- grepl("^marker$", cnv)
  lcond[["cond.markertypecol"]] <- grepl("^marker\\.type$", cnv)
  # evaluate regex
  cond.allcol <- lapply(lcond, function(ii){length(which(ii)) > 0})
  cond.allcol <- unlist(cond.allcol)
  # get final eval
  cond.output <- length(cond.allcol[cond.allcol])==4
  if(cond.output){
    return(TRUE)
  } else{
    return(FALSE)
  }
  return(NULL)
}

#' rand_donor_marker_table
#' 
#' Get a flat table of random donor marker signals by types.
#' 
#' @param ndonor Number of donors to simulate.
#' @param gindexv Vector of marker indices. Index values correspond to the k types,
#' and each index position represents a marker (e.g. c(1,2,2) means two markers 
#' for the second type, etc.).
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
#' rand_donor_marker_table(ndonor = 2, gindexv = c(1,2))
#' 
#' rand_donor_marker_table(ndonor = 10, gindexv = c(1,1,2))
#' 
#' rand_donor_marker_table(ndonor = 10, gindexv = c(rep(1, 10), rep(2, 20)))
#' 
#' @seealso random_lgv
#' @export
rand_donor_marker_table <- function(ndonor = 2, gindexv = c(1, 2), 
                                    lambda.pos = 20, lambda.neg = 2,
                                    sd.offset.pos = 10, sd.offset.neg = 2, 
                                    seed.num = 0, verbose = FALSE, ...){
  set.seed(seed.num)
  nmarkers <- length(gindexv); ktotal <- length(unique(gindexv))
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
    unlist(random_lgv(gindexv, num.iter = 1, lambda.pos = meanv.pos[ii],
                      lambda.neg = meanv.neg[ii], ...))
  }))
  md <- as.data.frame(md); colnames(md) <- paste0("donor", seq(ndonor))
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
  list(pca.bydonor = pca_bydonor(dt = dt, title.append = title.append, 
                                 verbose = verbose, ...), 
       pca.bydonortype = pca_bydonortype(dt = dt, title.append = title.append, 
                                         verbose = verbose, ...))
}

#' donordf_from_mexpr
#'
#' Makes a valid donordf object from an expression matrix (e.g. rows = 
#' markers/genes, columns = samples/type data)
#'
#' @param mexpr An expression matrix (rows = markers/genes, columns = 
#' samples/type data).
#' @param verbose Whether to show verbose status messages.
#' @returns 
#' @examples 
#' df <- rand_donor_marker_table()
#' madj <- donoradj_combat(df, return.type = "mexpr")
#' df.adj <- donordf_from_mexpr(mexpr = madj)
#' @export
donordf_from_mexpr <- function(mexpr){
  typev <- unique(gsub(".*;", "", colnames(madj)))
  df <- do.call(rbind, lapply(seq(nrow(madj)), function(markeri){
    mi <- madj[markeri,,drop=F]
    dfi <- do.call(rbind, lapply(typev, function(typei){
      mi[,grepl(paste0(".*;", typei), colnames(mi))]
    }))
    dfi <- as.data.frame(dfi); colnames(dfi) <- gsub(";.*", "", colnames(dfi))
    dfi$type <- typev; dfi$marker <- paste0("marker",markeri)
    dfi
  }))
  which.donor.cnv <- grepl("^donor\\d", colnames(df))
  df$donor.combn.all.mean <- rowMeans(as.matrix(df[,which.donor.cnv]))
  df$donor.combn.all.median <- rowMedians(as.matrix(df[,which.donor.cnv]))
  # final check
  if(check_donordf(df)){
    if(verbose){message("donordf conversion success. Returning...")}
    return(df)
  } else{
    stop("Error, couldn't convert mexpr to donordf. ",
         "Do the mexpr colnames have format donor;type?")
  }
  return(NULL)
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
#' donordf <- rand_donor_marker_table()
#' donorv <- donordf$donor.combn.all.mean
#' donoradj(donorv, donordf, method = "combat")
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
#' vector of adjusted donor signals, or "mexpr" for a matrix of adjusted 
#' signals.
#' @param verbose Whether to show verbose status messages.
#' @returns madj, matrix of adjusted marker signals.
#' @examples 
#' df <- rand_donor_marker_table()
#' donoradj_combat(df)
#' @export
donoradj_combat <- function(df, return.type = "donor.adj", verbose = FALSE){
  require(sva)
  # make expr matrix
  filt.donor <- grepl("donor\\d", colnames(df))
  mexpr <- do.call(rbind, lapply(unique(df$marker), function(mi){
    dff <- df[df$marker==mi, ]
    unlist(lapply(unique(dff[dff$marker==mi,]$type), function(ti){
      datv <- dff[dff$type==ti, filt.donor]
      names(datv) <- paste0(colnames(dff[,filt.donor]), ";", ti)
      return(datv)
    }))
  }))
  rownames(mexpr) <- unique(df$marker)
  # make pheno 
  cnv <- colnames(mexpr)
  pheno <- data.frame(donor = gsub(";.*", "", cnv),
                      type = gsub(".*;", "", cnv))
  # get combat vars
  mod <- model.matrix(~1, data = pheno)
  batch <- pheno$donor
  if(verbose){
    madj <- ComBat(dat = mexpr, batch = batch, mod = mod,
                   par.prior = TRUE, prior.plots = FALSE)
  } else{
    madj <- suppressWarnings(
      ComBat(dat = mexpr, batch = batch, mod = mod,
             par.prior = TRUE, prior.plots = FALSE)
    )
  }
  if(return.type == "donor.adj"){
    ltype <- list(type = gsub(".*;", "", colnames(madj)))
    dfa <- aggregate(t(madj), by = ltype, FUN = "mean")
    dfa <- t(dfa[,2:ncol(dfa)])
    donor.adj <- unlist(lapply(seq(nrow(dfa)), function(ri){dfa[ri,]}))
    return(donor.adj)
  } else if(return.type == "mexpr"){
    return(madj)
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