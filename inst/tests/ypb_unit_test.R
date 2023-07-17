#!/usr/bin/env R

#
# Tests y pseudobulk function
#

# dependencies
libv <- c("here", "dplyr", "ggplot2", "gridExtra", "SingleCellExperiment", 
          "SummarizedExperiment", "scran")
# library(lute)
sapply(libv, library, character.only = TRUE)

# mrb sce path
sce.mrb.path <- here("deconvo_method-paper", "outputs", "01_prepare-datasets", "sce-mrb_dlpfc.rda")
sce <- get(load(sce.mrb.path))
# donor filters
donor1.filter <- sce[["donor"]]=="donor1"

#---------------------------------------------------
# section 1 --- THIS SHOULD PASS WITH CONDITION FALSE
#---------------------------------------------------
# get ypb objects
ypb3 <- lute::ypb_from_sce(sce, "logcounts", "k2", "donor")
ypb1 <- lute::ypb_from_sce(sce[,donor1.filter], "logcounts", "k2", "donor")
# compare
identical(ypb1[,1], ypb3[,"donor1"]) 

#---------------------------------------------------
# section 2 --- try filter on 2 donors --- 
#---------------------------------------------------
# get ypb objects
donor1and2.filter <- sce[["donor"]] %in% c("donor1", "donor2")
ypb3 <- lute::ypb_from_sce(sce, "logcounts", "k2", "donor")
ypb2 <- lute::ypb_from_sce(sce[,donor1and2.filter], "logcounts", "k2", "donor")
# compare
identical(ypb2[,1], ypb3[,c("donor1")]) 

#---------------------------------------------------
# section 1b --- THIS SHOULD PASS WITH CONDITION TRUE
#---------------------------------------------------
# convert assays to matrices
logcounts(sce) <- as.matrix(logcounts(sce))
# get ypb objects
ypb3 <- lute::ypb_from_sce(sce, "logcounts", "k2", "donor")
ypb1 <- lute::ypb_from_sce(sce[,donor1.filter], "logcounts", "k2", "donor")
# compare
identical(ypb1[,1], ypb3[,"donor1"]) 

#------------------------------------------------------
# section 2 --- THIS SHOULD REFLECT RESULT OF SECTION 1
#------------------------------------------------------

library(lute)

# sce.all <- sce
# sce <- sce[,donor1.filter]

assay.name = "logcounts" 
celltype.variable = "k2" 
sample.id.variable = "donor"
S = NULL
sce.iter <- sce

# COMPARE ytable
sce <- sce.iter[,sce.iter[[sample.id.variable]]=="donor1"]
num.groups <- 1; unique.group.id.vector <- ""
if(!is(sample.id.variable, "NULL")){
  group.id.vector <- sce[[sample.id.variable]]
  unique.group.id.vector <- group.id.vector %>% unique() %>% as.character()
  num.groups <- unique.group.id.vector %>% length()
}
list.cell.types <- lute:::.get_celltypes_from_sce(
  sce = sce, celltype.variable = celltype.variable)
num.types <- list.cell.types[["unique.types"]] %>% length()
ypb.list <- lapply(unique.group.id.vector, function(group.id){
  sce.filter <- sce[,sce[[sample.id.variable]]==group.id]
  if(is(S, "NULL")){
    S <- rep(1, num.types)
    names(S) <- list.cell.types[["unique.types"]]
  }
  Znew <- lute:::.get_z_from_sce(sce.filter, assay.name, celltype.variable)
  P <- list.cell.types[["character"]] %>% table() %>% prop.table()
  order.p <- match(names(P), list.cell.types[["unique.types"]]) %>% order()
  P <- P[order.p]
  ZSnew <- lute:::.zstransform(Znew, S)
  ypb <- t(t(P) %*% t(ZSnew))
  return(ypb)
})
ypb.table1 <- do.call(cbind, ypb.list) %>% as.data.frame()
colnames(ypb.table1) <- unique.group.id.vector
sce <- sce.iter
num.groups <- 1; unique.group.id.vector <- ""
if(!is(sample.id.variable, "NULL")){
  group.id.vector <- sce[[sample.id.variable]]
  unique.group.id.vector <- group.id.vector %>% unique() %>% as.character()
  num.groups <- unique.group.id.vector %>% length()
}
list.cell.types <- lute:::.get_celltypes_from_sce(
  sce = sce, celltype.variable = celltype.variable)
num.types <- list.cell.types[["unique.types"]] %>% length()
ypb.list <- lapply(unique.group.id.vector, function(group.id){
  sce.filter <- sce[,sce[[sample.id.variable]]==group.id]
  
  if(is(S, "NULL")){
    S <- rep(1, num.types)
    names(S) <- list.cell.types[["unique.types"]]
  }
  Znew <- lute:::.get_z_from_sce(sce.filter, assay.name, celltype.variable)
  P <- list.cell.types[["character"]] %>% table() %>% prop.table()
  order.p <- match(names(P), list.cell.types[["unique.types"]]) %>% order()
  P <- P[order.p]
  ZSnew <- lute:::.zstransform(Znew, S)
  ypb <- t(t(P) %*% t(ZSnew))
  return(ypb)
})
ypb.table2 <- do.call(cbind, ypb.list) %>% as.data.frame()
colnames(ypb.table2) <- unique.group.id.vector

head(ypb.table1)
#donor1
#MIR1302-2HG 8.474155e-05
#FAM138A     0.000000e+00
#OR4F5       0.000000e+00
#AL627309.1  8.017225e-02
#AL627309.3  0.000000e+00
#AL627309.2  9.163554e-04

head(ypb.table2[,"donor1",drop=F])
# donor1
# MIR1302-2HG 0.0001722307
# FAM138A     0.0000000000
# OR4F5       0.0000000000
# AL627309.1  0.0916942252
# AL627309.3  0.0000000000
# AL627309.2  0.0007158646

plot(ypb.table1[,1], ypb.table2[,"donor1"])

# COMPARE z objects --- SHOULD BE IDENTICAL
sce <- sce.iter[,sce.iter[[sample.id.variable]]=="donor1"]
num.groups <- 1; unique.group.id.vector <- ""
if(!is(sample.id.variable, "NULL")){
  group.id.vector <- sce[[sample.id.variable]]
  unique.group.id.vector <- group.id.vector %>% unique() %>% as.character()
  num.groups <- unique.group.id.vector %>% length()
}
list.cell.types <- lute:::.get_celltypes_from_sce(
  sce = sce, celltype.variable = celltype.variable)
num.types <- list.cell.types[["unique.types"]] %>% length()
zlist1 <- lapply(unique.group.id.vector, function(group.id){
  sce.filter <- sce[,sce[[sample.id.variable]]==group.id]
  
  if(is(S, "NULL")){
    S <- rep(1, num.types)
    names(S) <- list.cell.types[["unique.types"]]
  }
  Znew <- lute:::.get_z_from_sce(sce.filter, assay.name, celltype.variable)
  return(Znew)
})
names(zlist1) <- unique.group.id.vector
sce <- sce.iter
num.groups <- 1; unique.group.id.vector <- ""
if(!is(sample.id.variable, "NULL")){
  group.id.vector <- sce[[sample.id.variable]]
  unique.group.id.vector <- group.id.vector %>% unique() %>% as.character()
  num.groups <- unique.group.id.vector %>% length()
}
list.cell.types <- lute:::.get_celltypes_from_sce(
  sce = sce, celltype.variable = celltype.variable)
num.types <- list.cell.types[["unique.types"]] %>% length()
zlist2 <- lapply(unique.group.id.vector, function(group.id){
  sce.filter <- sce[,sce[[sample.id.variable]]==group.id]
  
  if(is(S, "NULL")){
    S <- rep(1, num.types)
    names(S) <- list.cell.types[["unique.types"]]
  }
  Znew <- lute:::.get_z_from_sce(sce.filter, assay.name, celltype.variable)
  return(Znew)
})
names(zlist2) <- unique.group.id.vector
head(zlist1[["donor1"]])
head(zlist2[["donor1"]])
identical(zlist1[["donor1"]], zlist2[["donor1"]])

# Compare Znew objects --- should be identical
sce <- sce.iter[,sce.iter[[sample.id.variable]]=="donor1"]
num.groups <- 1; unique.group.id.vector <- ""
if(!is(sample.id.variable, "NULL")){
  group.id.vector <- sce[[sample.id.variable]]
  unique.group.id.vector <- group.id.vector %>% unique() %>% as.character()
  num.groups <- unique.group.id.vector %>% length()
}
list.cell.types <- lute:::.get_celltypes_from_sce(
  sce = sce, celltype.variable = celltype.variable)
num.types <- list.cell.types[["unique.types"]] %>% length()
zlist1 <- lapply(unique.group.id.vector, function(group.id){
  sce.filter <- sce[,sce[[sample.id.variable]]==group.id]
  
  if(is(S, "NULL")){
    S <- rep(1, num.types)
    names(S) <- list.cell.types[["unique.types"]]
  }
  Znew <- lute:::.get_z_from_sce(sce.filter, assay.name, celltype.variable)
  P <- list.cell.types[["character"]] %>% table() %>% prop.table()
  order.p <- match(names(P), list.cell.types[["unique.types"]]) %>% order()
  P <- P[order.p]
  ZSnew <- lute:::.zstransform(Znew, S)
  return(ZSnew)
})
names(zlist1) <- unique.group.id.vector
sce <- sce.iter
num.groups <- 1; unique.group.id.vector <- ""
if(!is(sample.id.variable, "NULL")){
  group.id.vector <- sce[[sample.id.variable]]
  unique.group.id.vector <- group.id.vector %>% unique() %>% as.character()
  num.groups <- unique.group.id.vector %>% length()
}
list.cell.types <- lute:::.get_celltypes_from_sce(
  sce = sce, celltype.variable = celltype.variable)
num.types <- list.cell.types[["unique.types"]] %>% length()
zlist2 <- lapply(unique.group.id.vector, function(group.id){
  sce.filter <- sce[,sce[[sample.id.variable]]==group.id]
  
  if(is(S, "NULL")){
    S <- rep(1, num.types)
    names(S) <- list.cell.types[["unique.types"]]
  }
  Znew <- lute:::.get_z_from_sce(sce.filter, assay.name, celltype.variable)
  P <- list.cell.types[["character"]] %>% table() %>% prop.table()
  order.p <- match(names(P), list.cell.types[["unique.types"]]) %>% order()
  P <- P[order.p]
  ZSnew <- lute:::.zstransform(Znew, S)
  return(ZSnew)
})
names(zlist2) <- unique.group.id.vector
head(zlist1[["donor1"]])
head(zlist2[["donor1"]])
identical(zlist1[["donor1"]], zlist2[["donor1"]])


# FUNCTION

ypb_from_sce <- function(sce, assay.name = "counts", 
                         celltype.variable = "celltype", 
                         sample.id.variable = NULL, S = NULL){
  require(dplyr)
  num.groups <- 1; unique.group.id.vector <- ""
  if(!is(sample.id.variable, "NULL")){
    group.id.vector <- sce[[sample.id.variable]]
    unique.group.id.vector <- group.id.vector %>% unique() %>% as.character()
    num.groups <- unique.group.id.vector %>% length()
  }
  list.cell.types <- lute:::.get_celltypes_from_sce(
    sce = sce, celltype.variable = celltype.variable)
  num.types <- list.cell.types[["unique.types"]] %>% length()
  ypb.list <- lapply(unique.group.id.vector, function(group.id){
    sce.filter <- sce[,sce[[sample.id.variable]]==group.id]
    
    if(is(S, "NULL")){
      S <- rep(1, num.types)
      names(S) <- list.cell.types[["unique.types"]]
    }
    Znew <- lute:::.get_z_from_sce(sce.filter, assay.name, celltype.variable)
    P <- list.cell.types[["character"]] %>% table() %>% prop.table()
    order.p <- match(names(P), list.cell.types[["unique.types"]]) %>% order()
    P <- P[order.p]
    ZSnew <- lute:::.zstransform(Znew, S)
    ypb <- t(t(P) %*% t(ZSnew))
    return(ypb)
  })
  
  
  ypb.table <- do.call(cbind, ypb.list) %>% as.data.frame()
  colnames(ypb.table) <- unique.group.id.vector
  
  
  return(ypb.table)
}

# inspect ypb table