#!/usr/bin/env R

# Author: Sean Maden
#
# Test methods for SingleCellExperiments
#

require(lute)

#---------------------
# test sce_groupstat()
#---------------------
# example

# make random data
sce = random_sce(zero.include = TRUE, zero.fract = 0.5)
colData(sce)$donor <- c(rep("donor1", 7), rep("donor2", 3))

# get group summaries across types
sce_groupstat(sce, summarytype = "colData", return.tall = TRUE)

# get group-gene summaries across types
sce_groupstat(sce, summarytype = "rowData", return.tall = TRUE)

# get tall table of group-type summaries
sce_groupstat(sce, type.variable = "celltype", 
              summarytype = "colData", return.tall = TRUE)

# get coldata-compatible group-type summaries
sce_groupstat(sce, type.variable = "celltype", 
              summarytype = "colData", return.tall = FALSE)

# get rowdata-compatible group-type summaries
sce_groupstat(sce, type.variable = "celltype", 
              summarytype = "rowData", return.tall = FALSE)

sce_groupstat(sce, summarytype = "rowData")
sce_groupstat(sce, type.variable = "celltype", return.type = "tall")
sce_groupstat(sce, summarytype = "rowData", return.type = "tall")



#sce
group.variable = "donor"
ugroupv = NULL
assayname = "counts"
summarytype = "rowData"
groupstat = c("count", "numzero", "var")
type.variable = NULL
verbose = FALSE

if(!is(sce, "SingleCellExperiment")){
  stop("Error, sce must be a SingleCellExperiment.")}
cd <- colData(sce)
# function
groupstat.filt <- groupstat %in% c("count", "mean", "median", 
                                   "var", "sd", "numzero")
groupstat <- groupstat[groupstat.filt]

if(verbose){message("Checking colData variables...")}
lvar <- check_coldata(cd = cd, var = c(group.variable, type.variable))
# get group stats df
ugroupv <- ufiltv <- lvar[[group.variable]]$uvec
groupv <- filtv <- lvar[[group.variable]]$vec
if(!is(type.variable, "NULL")){
  utypev <- lvar[[type.variable]]$uvec
  typev <- lvar[[type.variable]]$vec
  ufiltv <- paste0(rep(ugroupv, each = length(utypev)))
  ufiltv <- paste0(ufiltv, ";", utypev)
  filtv <- paste0(groupv, ";", typev)
}

ldf <- lapply(ugroupv, function(groupi){
  exprf <- assays(sce[,groupv==groupi])$counts
  if(is(type.variable, "NULL")){
    dfi <- get_groupstat_df(exprf, groupstat = groupstat, 
                            summarytype = summarytype)
    
  } else{
    dfi <- do.call(rbind, lapply(utypev, function(typei){
      exprff <- assays(sce[,groupv==groupi & typev==typei])$counts
      dfii <- get_groupstat_df(exprf, groupstat = groupstat, 
                               summarytype = summarytype)
      dfii$type <- typei; return(dfii)
    }))
  }
  dfi$group <- groupi
  dfi
})


if(return.type == "summarytype"){
  dfr <- do.call(cbind, lapply(ldf, function(dfi){
    groupi <- unique(dfi$group)[1]
    dfi <- dfi[,!colnames(dfi) %in% c("group", "marker")]
    colnames(dfi) <- paste0(groupi, ";", colnames(dfi))
    dfi
  }))
} else{
  dfr <- do.call(rbind, lapply(ldf, function(dfi)))
}

if(summarytype == "colData"){
  dfr <- do.call(cbind, lapply(ldf, function(dfi){
    colnames(dfi) <- paste0(unique(dfi$group), ";", colnames(dfi))
    dfi[,!colnames(dfi) %in% c("group", "marker")]
  }))
} else{
  dfr <- do.call(cbind, lapply(ldf, function(dfi){
    dfi[,!colnames(dfi)=="marker"]
  }))
  rownames(dfr) <- ldf[[1]]$marker
}


sce = random_sce()
colData(sce)$donor <- c(rep("donor1", 7), rep("donor2", 3))

sce_groupstat(sce, summarytype = "colData")
sce_groupstat(sce, summarytype = "rowData")
sce_groupstat(sce, type.variable = "celltype", summarytype = "colData")
sce_groupstat(sce, type.variable = "celltype", summarytype = "rowData")

# test parse_ldf()

parse_ldf <- function(ldf, summarytype)


#------------------------
# test get_groupstat_df()
#------------------------

# make random data
sce = random_sce(zero.include = TRUE, zero.fract = 0.5)
colData(sce)$donor <- c(rep("donor1", 7), rep("donor2", 3))

# params
exprf <- assays(sce)[["counts"]]
groupstat <- c("numzero", "count", "var")
summarytype <- "rowData"
round.digits = 3

# parse summarytype
ngene <- nrow(exprf); ncell <- ncol(exprf); exprff <- exprf
if(summarytype == "colData"){exprff <- t(as.matrix(colMeans(exprf)))}
# make na matrix
groupstatf <- groupstat[!grepl("^count.*|numzero", groupstat)]
mna <- matrix(NA, ncol = length(groupstatf), nrow = nrow(exprff))
dfti <- as.data.frame(mna); colnames(dfti) <- groupstatf
if(length(which(grepl("^count.*", groupstat))) > 0){
  dfti$count.cells <- ncell; dfti$count.genes <- ngene
}
for(ri in seq(nrow(exprff))){
  if(!is(exprff, "NULL")){
    # get group stats
    ei <- exprff[ri,,drop = F]
    if("mean" %in% groupstat){
      meani <- rowMeans(ei)
      dfti[ri,]$mean <- round(meani, digits = round.digits)
    }
    if("median" %in% groupstat){
      mediani <- rowMedians(ei)
      dfti[ri,]$median <- round(mediani, digits = round.digits)
    }
    if("var" %in% groupstat){
      vari <- rowVars(ei)
      dfti[ri,]$var <- round(vari, digits = round.digits)
    }
    if("sd" %in% groupstat){
      sdi <- rowSds(ei)
      dfti[ri,]$sd <- round(sdi, digits = round.digits)
    }
    if(length(which(grepl("^numzero.*", groupstat)))>0){
      if(summarytype == "colData"){
        # get summaries across cells
        num.cells.zero <- unlist(lapply(seq(nrow(exprf)), function(ri){
          datv <- exprf[ri,]; length(which(datv==0))
        }))
        zct.cell <- mean(num.cells.zero)
        zfr.cell <- mean(num.cells.zero/ncell)
        dfti$mean.count.cells.zero <- round(zct.cell, digits = round.digits)
        dfti$mean.fract.cells.zero <- round(zfr.cell, digits = round.digits)
        # get summaries across genes
        num.genes.zero <- unlist(lapply(seq(ncol(exprf)), function(ci){
          datv <- exprf[ci,]; length(which(datv==0))
        }))
        zct.gene <- mean(num.cells.zero)
        zfr.gene <- mean(num.cells.zero/ncell)
        dfti$mean.count.genes.zero <- round(zct.gene, digits = round.digits)
        dfti$mean.fract.genes.zero <- round(zfr.gene, digits = round.digits)
      } else{
        dfti$ncell.zero[ri] <- length(which(ei[ri,]==0))
      }
    }
  }
}
if(summarytype == "rowData"){dfti$marker <- rownames(exprf)}