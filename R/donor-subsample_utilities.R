#!/usr/bin/env R

# Author: Sean Maden
#
# Utilities to perform subsampling experiments on donors/samples/batches.
#


#' prepare_subsample_experiment
#'
#' @param sce SingleCellExperiment
#' @param iterations Total iterations to perform, per valid method specified in 
#' the methods argument.
#' @param methods Vector of valid deconvolution methods to test.
#' @param celltype.variable Variable containing cell type labels.
#' @param group.variable Variable containing group labels for donor/sample/batch.
#' @param fraction.cells Fraction of cells, per group and type, to subsample in each iteration.
#' @param number.of.groups Total groups to sample for each iteration.
#' @param scale.factor Vector of scale factors corresponding to alphebetical order of
#' unique type labels in celltype.variable.
#' @param seed.num Random seed for computational reproducibility.
#' @param which.save Which items to save when running this function. Either "sce" for 
#' SingleCellExperiment, "tp" for true proportions, "ypb" for pseudobulk, or "wt"
#' for the `r-nf` compatible workflow table. These objects correspond to filepaths
#' in save.paths (see details).
#' @param save.paths List of save path details for new files.
#' @returns List of iterations data, saving new files as side-effect.
#' @details Prepares files for a subsampling experiment to evaluate the impact 
#' of a group-wise bias effect on deconvolution outcomes. Either real or simulated
#' data can be provided to the sce argument. This is then randomly subsampled 
#' where remaining arguments specify details for random iterations. 
#' 
#' The fraction of cells per type and group is specified by the fraction.cells 
#' argument. Since the total cells per type may vary across groups, the total 
#' cells per iteration and type can also vary unless the sce data has been 
#' pre-filtered to control for this.
#' 
#' Since it is assumed you will use the same bulk data and true proportions to
#' evaluate iterations, these data are calculated once using the full sce object.
#' The pseudobulk dataset is calculated as the produce of the mean signature 
#' gene expression and the scale.factor cell size scale factor object. The true
#' cell type proportions correspond to this pseudobulk sample.
#'
#' The final output files are saved at the base.path specified in save.paths list.
#' These are suitable for use in an r-nf workflow, or using some other means of
#' parallelization of the iteration runs.
#'
#' @examples 
#' sce <- random_sce()
#' sce[["donor"]] <- "donor1"
#' iterations.list <- prepare_subsample_experiment(sce, 
#' celltype.variable = "celltype", group.variable = "donor")
#' names(iterations.list)
#' 
#' @export
prepare_subsample_experiment <- function(sce, iterations = 10, 
                                         methods = c("nnls", "music"),
                                         celltype.variable = "k2", 
                                         group.variable = "Sample",
                                         fraction.cells = 25,
                                         number.of.groups = 3,
                                         scale.factor = c(10, 2),
                                         seed.num = 0,
                                         which.save = c("sce", "tp", "ypb", "wt"),
                                         save.paths = list(base.path = "./data/", 
                                                           sce.name = "sce.rda",
                                                           wt.name = "workflow-table.rda",
                                                           tp.name = "true-proportions.rda")){
  # randomly take cells from 3 donors at a time
  # get the exact numbers of cells of each type to sample
  dft <- as.data.frame(table(sce[[celltype.variable]], 
                             sce[[sample.variable]]))
  num.cells.glial <- round(min(dft[dft[,1]=="glial",3])*0.3, 0)
  num.cells.neuron <- round(min(dft[dft[,1]=="neuron",3])*0.3, 0)
  # get random cell indices as a list
  cd <- colData(sce)
  unique.samples <- unique(cd[,sample.variable])
  unique.types <- unique(cd[,celltype.variable])
  unique.types <- unique.types[order(unique.types)]
  # order s cell size factors
  S <- S[order(match(names(S), unique.types))]
  num.cellsv <- c("glial" = num.cells.glial, "neuron" = num.cells.neuron)
  lindex <- lapply(seq(num.iter), function(ii){
    set.seed(ii)
    # get filtered sce data as scef
    random.samples <-  sample(unique.samples, num.sample.iter)
    filt <- cd[,sample.variable] %in% random.samples
    cdf <- cd[filt,]
    vindex <- which(colnames(sce) %in% unlist(lapply(unique.types, function(typei){
      sample(rownames(cdf[cdf[,celltype.variable]==typei,]), num.cellsv[typei])
    })))
    # return results
    list(vindex = vindex, samples = random.samples)
  })
  
  # get pseudobulk data
  tp <- as.data.frame(table(sce[[celltype.variable]]))
  tp.prop <- tp[,2]; names(tp.prop) <- tp[,1]
  P <- tp.prop/sum(tp.prop)
  Z <- do.call(cbind, lapply(unique.types, function(typei){
    rowMeans(assays(sce)[[assay.name]])
  }))
  ZS <- sweep(Z, 2, S, "*")
  ypb <- t(t(P) %*% t(ZS))
  
  # save new data
  # save indices
  save(lindex, file = lindex.fpath)
  # save pseudobulk
  save(ypb, file = ypb.fpath)
  # save tp
  save(P, file = tp.fpath)
  
  #---------------------
  # write workflow table
  #---------------------
  wt <- do.call(rbind, lapply(methodv, function(methodi){
    wti <- data.frame(iterations_index = seq(num.iter))
    wti$method <- methodi
    wti$sample_id <- unlist(lapply(lindex, function(li){
      paste0(li$samples, collapse = ";")}))
    wti$celltype_variable <- celltype.variable
    wti$assay_name <- assay.name
    # manage filepaths
    cnamev <- c("sce_filepath", "bulk_filepath", "list_index_filepath",
                "true_proportions_filepath")
    fpathv <- c(file.path("data", sce.fname),
                file.path("data", ypb.fname),
                file.path("data", lindex.fname),
                file.path("data", tp.fname))
    for(ii in seq(length(cnamev))){
      wti[,cnamev[ii]] <- paste0('"$launchDir/', fpathv[ii], '"')}
    wti
  }))
  
  # save
  write.csv(wt, file = wt.fpath, row.names = F)
}

#'
#'
#'
#'
analyze_subsample_results <- function(){
  
}

#'
#'
#'
#'
#'
subsample_summary <- function(){
  
}

#'
#'
#'
#'
#'
subsample_plots <- function(){
  
}

