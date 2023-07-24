#!/usr/bin/env R

# 
# Functions supporting r-nf_deconvolution workflow runs.
#

#' new_workflow_table
#' 
#' Makes a new experiment table for r-nf_deconvolution runs. 
#' 
#' @param sce.names Names of SingleCellExperiment files to load.
#' @param data.directory Directory containing datasets to load.
#' @param true.proportions.filename.stem File name stem of true proportions values.
#' @param celltype.variable Name of variable containing cell type labels.
#' @param table.directory Directory to write new table.
#' @param table.filename Filename of new table to write.
#' @param overwrite Whether to overwrite old table files.
#' @param verbose Whether to show verbose messages.
#' @details Makes and returns/saves a r-nf_deconvolution experiment table. 
#' Checks for existence of provided files.
#' @returns New r-nf_deconvolution compatible table of experiment/run metadata.
#' 
#' @examples
#' new_workflow_table(save = F)
#'
#' @export
new_workflow_table <- function(sce.names = NULL, data.directory = "data",
                               true.proportions.filename.stem = "true_proportions_",
                               celltype.variable = "celltype", table.directory = ".", 
                               table.filename = "workflow_table.csv", save = TRUE, 
                               overwrite = TRUE, verbose = FALSE){
  rnf.colnames <- c("sce_filepath", "true_proportions_path", "decon_method", 
                    "decon_args", "run_info", "assay_name", "celltype_variable")
  table.filepath <- file.path(table.directory, table.filename)
  if(file.exists(table.filepath) & !overwrite){
    stop("Found existing workflow table at path ", table.fpath, ".")}
  dfnew <- matrix(nrow = 0, ncol = length(rnf.colnames))
  newline <- c(file.path("$launchDir", data.directory),
               file.path("$launchDir", data.directory), "nnls", "NA", 
               "lung_adeno_first_benchmark", "counts", celltype.variable)
  check.files <- TRUE
  if(is(sce.names, "NULL")){
    if(verbose){message("Making example table...")}
    sce.names <- "[SCE_FILENAME_HERE]"
    check.files <- FALSE
  }
  for(scei in sce.names){
    if(verbose){message("Working on sce object ", scei, "...")};linei <- newline
    sce.filepath <- file.path(data.directory, paste0(scei, ".rda"))
    sce.exists <- file.exists(sce.filepath)
    tp.filepath <- paste0(true.proportions.filename.stem, scei, ".rda")
    tp.exists <- file.exists(file.path(data.directory, tp.filepath))
    if(check.files){
      if(!sce.exists){
        if(verbose){
          message("Didn't find file ",sce.filepath,". Skipping data write.")}
      } else if(!tp.exists){
        if(verbose){
          message("Didn't find file ",tp.filepath,". Skipping data write.")}
      } else{if(verbose){message("Continuing.")}}
    }
    linei[1] <- file.path(linei[1], paste0(scei, ".rda"))
    linei[2] <- file.path(linei[2], 
                          paste0(true.proportions.filename.stem, scei, ".rda"))
    linei[1] <- paste0('"', linei[1], '"');linei[2] <-paste0('"', linei[2], '"')
    dfnew <- rbind(dfnew, linei)
  }
  colnames(dfnew) <- rnf.colnames
  if(save){write.csv(dfnew, file = table.filepath, row.names = FALSE)}
  return(dfnew)
}
