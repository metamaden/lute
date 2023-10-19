#!/usr/bin/env R

### 
### Functions supporting r-nf_deconvolution workflow runs.
###

#' new_workflow_table
#' 
#' Makes a new experiment table for r-nf_deconvolution runs. 
#' 
#' @param singleCellExperimentNames Names of SingleCellExperiment files to load.
#' @param dataDirectory Directory containing datasets to load.
#' @param trueProportionsFilenameStem File name stem of true proportions values.
#' @param cellTypeVariable Name of variable containing cell type labels.
#' @param tableDirectory Directory to write table.
#' @param tableFileName The file name of the new table to write.
#' @param save Whether to save the new table.
#' @param overwrite Whether to overwrite old table files.
#' @param verbose Whether to show verbose messages (T/F).
#' @details Makes and returns/saves a r-nf_deconvolution experiment table. 
#' Checks for existence of provided files.
#' @returns New r-nf_deconvolution compatible table of experiment/run metadata.
#' @importFrom utils write.csv
#' 
#' @examples
#' new_workflow_table(save=FALSE)
#'
#' @export
new_workflow_table <- function(singleCellExperimentNames=NULL, dataDirectory="data",
                               trueProportionsFilenameStem="true_proportions_",
                               cellTypeVariable="celltype", tableDirectory=".", 
                               tableFileName="workflow_table.csv", save=TRUE, 
                               overwrite=TRUE, verbose=FALSE){
  rnf.colnames <- c("sce_filepath", "true_proportions_path", "decon_method", 
                    "decon_args", "run_info", "assay_name", "celltype_variable")
  tableFilePath <- file.path(tableDirectory, tableFileName)
  if(file.exists(tableFilePath) & !overwrite){
    stop("Found existing table ", tableFilePath, ". Stopping.")}
  newDataFrame <- matrix(nrow=0, ncol=length(rnf.colnames))
  newline <- c(file.path("$launchDir", dataDirectory),
               file.path("$launchDir", dataDirectory), "nnls", "NA", 
               "lung_adeno_first_benchmark", "counts", cellTypeVariable)
  checkFiles <- TRUE
  if(is(singleCellExperimentNames, "NULL")){
    if(verbose){message("Making example table...")}
    singleCellExperimentNames <- "[SCE_FILENAME_HERE]"
    checkFiles <- FALSE
  }
  for(iterateSingleCellExperiment in singleCellExperimentNames){
    if(verbose){
      message("Working on sce object ", iterateSingleCellExperiment, "...")}
    lineIterate <- newline
    singleCellExperimentFilepath <- 
      file.path(dataDirectory, paste0(iterateSingleCellExperiment, ".rda"))
    singleCellExperimentExists <- file.exists(singleCellExperimentFilepath)
    trueProportionsFilepath <- 
      paste0(trueProportionsFilenameStem, iterateSingleCellExperiment, ".rda")
    trueProportionsExists <- 
      file.exists(file.path(dataDirectory, trueProportionsFilepath))
    if(checkFiles){
      if(!singleCellExperimentExists){
        if(verbose){
          message("Didn't find file ",
                  singleCellExperimentFilepath,". Skipping data write.")}
      } else if(!trueProportionsExists){
        if(verbose){
          message("Didn't find file ",
                  trueProportionsFilepath,". Skipping data write.")}
      } else{if(verbose){message("Continuing.")}}
    }
    lineIterate[1] <- file.path(
      lineIterate[1], paste0(iterateSingleCellExperiment, ".rda"))
    lineIterate[2] <- file.path(
      lineIterate[2], 
      paste0(trueProportionsFilenameStem, iterateSingleCellExperiment, ".rda"))
    lineIterate[1] <- paste0('"', lineIterate[1], '"')
    lineIterate[2] <-paste0('"', lineIterate[2], '"')
    newDataFrame <- rbind(newDataFrame, lineIterate)
  }
  colnames(newDataFrame) <- rnf.colnames
  if(save){write.csv(newDataFrame, file=tableFilePath, row.names=FALSE)}
  return(newDataFrame)
}
