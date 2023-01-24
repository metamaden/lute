#!/usr/bin/env R

# Author: Sean Maden
#
# Methods for table objects of type `donordf`.
#
# Sections:
#
# 1. checks: functions that check the object classes
#
# 2. conversions: functions facilitating conversion between object classes
#

#----------
# 1. checks
#----------
# functions that check the object classes.

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

#---------------
# 2. conversions
#---------------
# functions facilitating conversion between object classes

#' donordf_from_mexpr
#'
#' Makes a valid donordf object from an expression matrix (e.g. rows = 
#' markers/genes, columns = samples/type data)
#'
#' @param mexpr An expression matrix (rows = markers/genes, columns = 
#' samples/type data).
#' @param verbose Whether to show verbose status messages.
#' @returns mexpr, a new expression matrix
#' @examples 
#' df <- rand_donor_marker_table()
#' madj <- donoradj_combat(df, return.type = "mexpr")
#' df.adj <- donordf_from_mexpr(mexpr = madj)
#' @export
donordf_from_mexpr <- function(mexpr, verbose = FALSE){
  typev <- unique(gsub(".*;", "", colnames(mexpr)))
  df <- do.call(rbind, lapply(seq(nrow(mexpr)), function(markeri){
    mi <- mexpr[markeri,,drop=F]
    dfi <- do.call(rbind, lapply(typev, function(typei){
      mi[,grepl(paste0(".*;", typei), colnames(mi))]
    }))
    dfi <- as.data.frame(dfi); colnames(dfi) <- gsub(";.*", "", colnames(dfi))
    dfi$type <- typev; dfi$marker <- paste0("marker",markeri)
    dfi
  }))
  # get donor summaries
  which.donor.cnv <- grepl("^donor\\d", colnames(df))
  df$donor.combn.all.mean <- rowMeans(as.matrix(df[,which.donor.cnv]))
  df$donor.combn.all.median <- rowMedians(as.matrix(df[,which.donor.cnv]))
  # get marker.type from means
  markerv <- unique(df$marker)
  mapv <- unlist(lapply(markerv, function(mi){ # get marker.type mappings
    dff <- df[df$marker == mi,]
    max.val <- max(dff[,"donor.combn.all.mean"])
    max.filt <- dff$donor.combn.all.mean==max.val
    dff[max.filt,]$type[1]
  })); names(mapv) <- markerv
  df$marker.type <- "NA"
  for(mi in markerv){df[df$marker == mi,]$marker.type <- as.character(mapv[mi])}
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

#' mexpr_from_donordf
#'
#' Make and expression matrix from a donordf data.frame. The donordf contains 
#' one row per marker;type, while the expression matrix contains one row per
#' marker and one column per donor;type.
#'
#' @param df A donordf type data.frame.
#' @returns mexpr, a new expression matrix (rows = genes/markers, columns = 
#' samples/types)
#' @export
mexpr_from_donordf <- function(df){
  filt.donor <- grepl("donor\\d", colnames(df)) # get donor signals
  mexpr <- do.call(rbind, lapply(unique(df$marker), function(mi){
    dff <- df[df$marker==mi, ]
    ld <- lapply(unique(dff[dff$marker==mi,]$type), function(ti){
      datv <- dff[dff$type==ti, filt.donor]
      if(nrow(datv) > 1){
        stop("Error, more than one unique marker value for each type. ",
             "Does donordf contain multiple experiment groups?")}
      names(datv) <- paste0(colnames(dff[,filt.donor]), ";", ti)
      datv
    })
    unlist(ld)
  }))
  rownames(mexpr) <- unique(df$marker)
  return(mexpr)
}