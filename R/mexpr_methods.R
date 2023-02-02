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
#'
#'
get_anova_df <- function(sce, ngene.sample = 2000, 
                         model = "expr ~ celltype * donor",
                         return.var = c("sumsq", "perc.var"), 
                         type = "complex"){
  lr <- lgg <- list()
  sampv <- sample(seq(nrow(sce)), ngene.sample, replace = T)
  dfa.all <- do.call(rbind, lapply(names(assays(sce)), function(ai){
    message("working on assay: ", ai, "...")
    mi <- assays(sce)[[ai]];
    mi <- mi[sampv,] # get random subset
    maxv <- rowMaxs(mi); mi <- mi[maxv > 0,] # filter all-zeros
    dfi <- data.frame(expr = mi[1,], 
                      celltype = sce[["k2"]],
                      donor = sce[["donor"]])
    dfa.mi <- do.call(rbind, lapply(seq(nrow(mi)), function(ii){
      dfi$expr <- mi[ii,]
      lmi <- paste0("lm(",model,", data = dfi)")
      lmi <- eval(parse(text(lmi)))
      avi <- anova(lmi)
      perc.var <- 100*avi$`Sum Sq`/sum(avi$`Sum Sq`)
      data.frame(perc.var.celltype = perc.var[1],
                 perc.var.donor = perc.var[2],
                 perc.var.celltype.donor = perc.var[3],
                 perc.var.residuals = perc.var[4],
                 marker = rownames(mi)[ii])
    }))
    message("finished anovas")
    dfa.mi <- dfa.mi[!is.na(dfa.mi[,1]),] # filter nas
    dfa.mi$assay <- ai; dfa.mi
  }))
  lr[["dfa"]] <- dfa.all
  
  lgg <- list(gg.jitter = lggj, gg.pt = lggpt)
  lr <- list(dfa = dfa.all, lgg = lgg)
  return(lr)
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



