#!/usr/bin/env R

# Author: Sean Maden
#
# Methods for expression assay matrices (rows = genes/markers, columns = 
# cells/samples/types)
#

#---------------
# anova analysis
#---------------

#' get_anova_df
#'
#' Get a data.frame of ANOVA statistics.
#' 
#' @param mexpr
#' @param ncell.sample
#' @param modelv Character vector containing the labels to use in the model.
#'
#'
get_anova_df <- function(mi, ncell.sample = 2000, type = "complex"){
  lr <- lgg <- list()
  sampv <- sample(seq(nrow(sce)), ncell.sample, replace = T)
  dfa.all <- do.call(rbind, lapply(names(assays(sce)), function(ai){
    message("working on assay: ", ai, "...")
    mi <- assays(sce)[[ai]];
    mi <- mi[sampv,] # get random subset
    maxv <- rowMaxs(mi); mi <- mi[maxv > 0,] # filter all-zeros
    dfi <- data.frame(expr = mi[1,], celltype = sce[["k2"]],
                      donor = sce[["donor"]])
    dfa.mi <- do.call(rbind, lapply(seq(nrow(mi)), function(ii){
      dfi$expr <- mi[ii,]
      lmi <- lm(expr ~ celltype * donor, data = dfi)
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
  # plots
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
  
  # scatter plot of perc var, within assay
  assayv <- c("counts", "logcounts")
  lggpt <- lapply(assayv, function(ai){
    df1 <- dfa.all[dfa.all$assay==ai,]
    df2 <- dfa.all[dfa.all$assay==
                     paste0(assayname1, "_ds-donor_combat"),]
    df1 <- df1[!duplicated(df1$marker),]
    df2 <- df2[!duplicated(df2$marker),]
    markerv.int <- intersect(df1$marker, df2$marker)
    df1 <- df1[df1$marker %in% markerv.int,]
    df2 <- df2[df2$marker %in% markerv.int,]
    df2 <- df2[order(match(df2$marker, df1$marker)),]
    cond <- identical(as.character(df2$marker), 
                      as.character(df1$marker))
    if(cond){
      cnv <- colnames(df1)[seq(4)]
      lggi <- lapply(cnv, function(ci){
        dfp <- data.frame(pv1 = df1[,ci],
                          pv2 = df2[,ci])
        ggplot(dfp, aes(x = pv1, y = pv2)) + geom_point(alpha = 0.5) + 
          geom_abline(slope = 1, intercept = 0, col = "black") +
          ggtitle("Perc. var. ", ci) + xlab(paste0(ai)) +
          ylab(paste0(ai, " ds-donor, combat"))
      })
      names(lggi) <- cnv
      return(lggi)
    }
  })
  names(lggpt) <- assayv
  lgg <- list(gg.jitter = lggj, gg.pt = lggpt)
  lr <- list(dfa = dfa.all, lgg = lgg)
  return(lr)
}