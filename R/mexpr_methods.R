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
sce <- random_sce()
sce[["donor"]] <- c(rep("donor1", 2), rep("donor2", 8), rep("donor1", 2))
pheno.df <- data.frame(donor = sce[["donor"]], celltype = sce[["celltype"]])
dfa <- get_anova_df(sce = sce, pheno.df = pheno.df)
#' @export
get_anova_df <- function(sce, pheno.df, ngene.sample = NULL, assayv = "counts",
                         model = "expr ~ celltype * donor",
                         return.var = c("percvar"), 
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

#---------------------------
# dispersion point estimates
#---------------------------

#'
#'
#'
analyze_dispersion_est <- function(){
  
}

#'
#'
#'
#'
est_dispersions <- function(){
  
}

#'
#'
#'
#'
plot_dispersion_coeff <- function(){
  
}

# get dispersions by type
# compare dispersion coefficients from fitted neg binom 
# parse params
typevar <- "k2"; marker.name <- "k2_top20"
assay.name <- "counts_adj"
# get bg genes
num.genes.bg <- 1000
bg.name <- paste0("bg_", num.genes.bg)
genes.samplev <- sample(seq(nrow(sce)), num.genes.bg)
# get marker genes
genes.markerv <- metadata(sce)[["k2.markers"]][["top20"]]$gene
# define categories
catv <- c(unique(typev), "all") 
# get plot data
dfp <- do.call(rbind, lapply(catv, function(typei){
  # parse filter
  type.filt <- seq(ncol(sce))
  if(!typei == "all"){type.filt <- sce[[typevar]] == typei}
  scef <- sce[,type.filt]
  
  # get dispersions
  mexpr <- assays(scef)[[assay.name]]
  lglm.bg <- glm_gp(mexpr[genes.samplev,], on_disk = F)
  lglm.top20 <- glm_gp(mexpr[genes.markerv,], on_disk = F)
  
  # get plot data
  dfp1 <- data.frame(disp = lglm.bg$overdispersions)
  dfp1$marker.type <- bg.name
  dfp2 <- data.frame(disp = lglm.top20$overdispersions)
  dfp2$marker.type <- marker.name
  dfp <- rbind(dfp1, dfp2)
  dfp$celltype <- typei
  return(dfp)
}))
# set return list
ldisp <- list(dfp = dfp)

# boxplots at 3 zoom levels
ldisp[["ggbox"]] <- list()
ldisp[["ggbox"]][["zoom1"]] <- ggplot(dfp, aes(x = marker.type, y = disp)) + 
  geom_boxplot() + facet_wrap(~celltype)
ldisp[["ggbox"]][["zoom2"]] <- ggplot(dfp, aes(x = marker.type, y = disp)) + 
  geom_boxplot() + facet_wrap(~celltype) + ylim(0, 350)
ldisp[["ggbox"]][["zoom3"]] <- ggplot(dfp, aes(x = marker.type, y = disp)) + 
  geom_boxplot() + facet_wrap(~celltype) + ylim(0, 50)

# jitter plots at 3 zoom levels
ldisp[["ggjitter"]] <- list()
ldisp[["ggjitter"]][["zoom1"]] <- ggplot(dfp, aes(x = marker.type, y = disp)) + 
  geom_jitter(alpha = 0.5) + 
  stat_summary(geom = "crossbar", fun = "median", color = "red") + 
  facet_wrap(~celltype)
ldisp[["ggjitter"]][["zoom2"]] <- ggplot(dfp, aes(x = marker.type, y = disp)) + 
  geom_jitter(alpha = 0.5) + 
  stat_summary(geom = "crossbar", fun = "median", color = "red") + 
  facet_wrap(~celltype) + ylim(0, 350)
ldisp[["ggjitter"]][["zoom3"]] <- ggplot(dfp, aes(x = marker.type, y = disp)) + 
  geom_jitter(alpha = 0.5) + 
  stat_summary(geom = "crossbar", fun = "median", color = "red") + 
  facet_wrap(~celltype) + ylim(0, 50)



