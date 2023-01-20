#!/usr/bin/env R

# Author: Sean Maden
#
# Test functions for NextFlow workflow.
#

require(lute)

# test filter_value_type()
set.seed(0)
sce <- random_sce(zero.include = T, zero.fract = 0.3)
scef <- filter_value_type(sce, filter.term = "zerocount", verbose = T,
                          max.gene.value.freq = 0.15)
metadata(scef)$filter.zerocount.by.type$df.type

# test filter_value_cells()
sce <- random_sce(zero.include = T, zero.fract = 0.3)
filter.term = "zerocount"
remove.cells = TRUE
max.value.freq = 0.25
assayname = "counts"
append.metadata = TRUE
new.metadata.name = NULL
verbose = T

if(!is(sce, "SingleCellExperiment")){
  stop("Error, sce must be a SingleCellExperiment.")}
cd <- colData(sce); mexpr <- assays(sce)[[assayname]]
value.countv <- apply(mexpr, 2, function(ci){
  if(filter.term == "NA"){
    length(which(is.na(ci)))
  } else{
    length(which(ci == 0))
  }
})
value.freqv <- value.countv/nrow(mexpr)
filt.cellv <- value.freqv > max.value.freq
if(remove.cells){
  scef <- sce[,!filt.cellv];filt.type <- "removed"
} else{
  scef <- sce; filt.type <- "flagged"
}
if(verbose){message("Filter on cell ",filter.term," values ",
                    filt.type, " ", length(which(filt.cellv)),
                    " cells.")}
if(append.metadata){
  if(verbose)(message("Appending new metadata."))
  lparam <- list(filter.term = filter.term,
                 max.value.freq = max.value.freq, 
                 assayname = assayname)
  df.cell <- data.frame(
    cell.uid = colnames(mexpr), 
    value.count = value.countv,
    value.freq = value.freqv, 
    above.max.value.freq = filt.cellv
  )
  lmd <- list(parameters = lparam, df.cell = df.cell)
  if(is(new.metadata.name, "NULL")){
    new.metadata.name <- paste0("filter.", filter.term, ".by.cell")
  }
  metadata(scef)[[new.metadata.name]] <- lmd
}

