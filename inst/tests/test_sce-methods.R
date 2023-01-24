#!/usr/bin/env R

# Author: Sean Maden
#
# Test methods for SingleCellExperiments
#

#---------------------
# test sce_groupstat()
#---------------------
sce = random_sce()
group.variable = "donor"
type.variable = "celltype"
ugroupv = NULL
utypev = NULL
assayname = "counts"
summarytype = "colData"
groupstat = c("count", "numzero", "var")
verbose = TRUE

# example
sce = random_sce()
colData(sce)$donor <- c(rep("donor1", 7), rep("donor2", 3))

# function
groupstat.filt <- groupstat %in% c("count", "mean", "median", "var", "sd", "numzero")
groupstat <- groupstat[groupstat.filt]

cd <- colData(sce)

# parse group variable
if(group.variable %in% colnames(cd)){
  groupv <- cd[, group.variable]
  if(is(ugroupv, "NULL")){
    ugroupv <- unique(groupv)
  } else{
    ugroupv <- ugroupv[ugroupv %in% groupv]
    if(length(ugroupv)==0){
      stop("Error, no labels in ugroupv found in sce colData")}
  }
} else{
  stop("Error, didn't find group.variable in sce data.")
}

# parse type variable
if(type.variable %in% colnames(cd)){
  typev <- cd[, type.variable]
  if(is(utypev, "NULL")){
    utypev <- unique(typev)
  } else{
    utypev <- utypev[utypev %in% typev]
    if(length(utypev)==0){
      stop("Error, no labels in ugroupv found in sce colData")}
  }
} else{
  stop("Error, didn't find group.variable in sce data.")
}

get_groupstat_df <- function(exprf, groupstat, summarytype){
  # make na matrix
  mna <- matrix(NA, ncol = length(groupstat), 
                nrow = ifelse(summarytype=="colData", 1, nrow(sce)))
  dfti <- as.data.frame(mna); colnames(dfti) <- groupstat
  if(!is(exprf, "NULL")){
    # get group stats
    if("count" %in% groupstat){
      dfti$count <- ifelse(summarytype == "colData", ncol(exprf), nrow(exprf))
    }
    if("mean" %in% groupstat){dfti$mean <- rowMeans(exprf)}
    if("median" %in% groupstat){dfti$median <- rowMedians(exprf)}
    if("var" %in% groupstat){dfti$var <- rowVars(exprf)}
    if("sd" %in% groupstat){dfti$sd <- rowSds(exprf)}
    if("numzero" %in% groupstat){
      dfti$numzero <- unlist(
        apply(exprf, 1, function(ri){length(which(ri==0))}))
    }
  }
  return(dfti)
}

# get group stats df
filtvar <- ugroupv; filtv <- groupv
if(!is(type.variable, "NULL")){
  ufiltvar <- paste0(rep(ugroupv, each = length(utypev))
                    , ";", utypev)
  filtv <- paste0(groupv, ";", typev)
}

dfg <- do.call(cbind, lapply(ugroupv, function(groupi){
  ufiltvari <- ufiltv[grepl(paste0(groupi, ";"), ufiltv)]
  dfgi <- do.call(rbind, lapply(ufiltvari, function(filti){
    scef <- sce[, filtv==filti]
    if(ncol(scef) > 0){
      exprf <- assays(scef)[[assayname]]
      if(summarytype == "colData"){
        exprf <- matrix(colMeans(exprf), nrow = 1)
      } else{exprf <- t(exprf)}
    } else{
      exprf <- NULL
    }
    get_groupstat_df(exprf, groupstat = groupstat, summarytype = summarytype)
  }))
  colnames(dfgi) <- paste0(groupi, ";", colnames(dfgi)); dfgi
}))
if(!is(type.variable, "NULL")){dfg$type <- utypev}
return(dfg)


if(group.variable %in% colnames(colData(sce))){
  
} else{
  message("Warning: variable ",group.variable," not found in sce coldata.")
}