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
#' @param sce
#' @param ngene.sample Number of genes to sample at random.
#' @param model Character string of the model to use. By default, provides an
#' interaction between donor and cell type, and the response/dependent variable
#' is the individual gene.
#' @param seed.num Random seed to set for reproducibility.
#' @param 
#' @returns 
#' @examples 
#' sce <- random_sce()
#' @export
analyze_anova <- function(sce, ngene.sample = 2000, 
                          model = "expr ~ celltype * donor",
                          return.var = c("sumsq", "perc.var"), 
                          type = "complex"){
  
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
#' @returns 
#' @examples 
#' sce <- random_sce()
#' sce[["donor"]] <- c(rep("donor1", 2), rep("donor2", 8), rep("donor1", 2))
#' pheno.df <- data.frame(donor = sce[["donor"]], celltype = sce[["celltype"]])
#' dfa <- get_anova_df(sce = sce, pheno.df = pheno.df)
#' @export
get_anova_df <- function(sce, pheno.df, ngene.sample = NULL, assayv = "counts",
                         model = "expr ~ celltype * donor",
                         return.var = c("perc.var"), 
                         seed.num = 0, verbose = FALSE){
  set.seed(seed.num)
  lr <- lgg <- lggj <- lggpt <- list()
  if(is(ngene.sample, "NULL")){
    sampv <- seq(nrow(sce))
  } else{
    sampv <- sample(seq(nrow(sce)), ngene.sample, replace = FALSE)
  }
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
      if("perc.var" %in% return.var){
        ssqv <- summary(avi)[[1]][[2]] # get sum of squared variances
        perc.var <- 100*ssqv/sum(ssqv)
        for(ii in seq(length(perc.var))){
          dfi[,ncol(dfi) + 1] <- perc.var[ii]
          colnames(dfi)[ncol(dfi)] <- paste0("perc.var.", namev[ii])
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


#'
#'
#'
#'
anova_jitter_plots <- function(){
  lggj <- list()
  # plot celltype perc var
  lggj[["celltype"]] <- 
    ggplot(dfa.all, aes(x = assay, y = perc.var.celltype)) +
    geom_jitter(alpha = 0.5) + 
    stat_summary(geom = "crossbar", fun = "mean", color = "red") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_zoom(ylim = c(0, 10))
  # plot donor perc var
  lggj[["donor"]] <- ggplot(dfa.all, aes(x = assay, y = perc.var.donor)) +
    geom_jitter(alpha = 0.5) + 
    stat_summary(geom = "crossbar", fun = "mean", color = "red") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_zoom(ylim = c(0, 2))
  # plot celltype.donor perc var
  lggj[["celltype.donor"]] <- 
    ggplot(dfa.all, aes(x = assay, y = perc.var.celltype.donor)) +
    geom_jitter(alpha = 0.5) + 
    stat_summary(geom = "crossbar", fun = "mean", color = "red") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_zoom(ylim = c(0, 2))
  # plot residuals perc var  
  lggj[["residuals"]] <- 
    ggplot(dfa.all, aes(x = assay, y = perc.var.residuals)) +
    geom_jitter(alpha = 0.5) + 
    stat_summary(geom = "crossbar", fun = "mean", color = "red") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_zoom(ylim = c(90, 100))
}

#' anova_scatter_plots
#'
#' @param dfa.all Data.frame of ANOVA results, returned from get_anova_df().
#' @param a1 Name of first assay, for x-axis.
#' @param a2 Name of second assay, for y-axis.
#' @param regex.cnv Regex pattern to find column names to plot
#'
anova_scatter_plots <- function(dfa.all, a1 = "counts", a2 = "counts_ds_combat", 
                                regex.cnv = "^perc\\.var\\..*"){
  # get dfs
  df1 <- dfa.all[dfa.all$assay==a1,]
  df2 <- dfa.all[dfa.all$assay==a2,]
  # match dfs
  df1 <- df1[!duplicated(df1$marker),]
  df2 <- df2[!duplicated(df2$marker),]
  markerv.int <- intersect(df1$marker, df2$marker)
  df1 <- df1[df1$marker %in% markerv.int,]
  df2 <- df2[df2$marker %in% markerv.int,]
  df2 <- df2[order(match(df2$marker, df1$marker)),]
  # get plots if match successful
  cond <- identical(as.character(df2$marker), as.character(df1$marker))
  if(cond){
    # get colnames to plot
    cnv <- colnames(df1); cnv <- cnv[grepl(regex.cnv, cnv)]
    # plot each colname
    lggi <- lapply(cnv, function(ci){
      dfp <- data.frame(pv1 = df1[,ci], pv2 = df2[,ci])
      title.str <- paste0("Percent var. ", gsub(".*\\." , "", ci))
      ggplot(dfp, aes(x = pv1, y = pv2)) + geom_point(alpha = 0.5) + 
        geom_abline(slope = 1, intercept = 0, col = "black") +
        ggtitle(title.str) + xlab(a1) + ylab(a2)
    })
    names(lggi) <- cnv
    return(lggi)
  } else{
    stop("Error, couldn't match markers across provided assays.")
  }
  return(NULL)
}



