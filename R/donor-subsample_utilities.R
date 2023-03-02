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
#' @param groups.per.iteration Number of groups to sample per iteration. Total 
#' groups should exceed this value, otherwise data is treated as a single group.
#' @param method.vector Vector of valid deconvolution methods to test.
#' @param celltype.variable Variable containing cell type labels.
#' @param assay.name Name of assay in sce to use for pseudobulk.
#' @param group.variable Variable containing group labels for donor/sample/batch.
#' @param cell.proportions Cell type proportions used to calculate number of cells
#' of each type to take, per iteration.
#' @param count.minimum Number of cells to take for the least abundant cell type,
#' according to proportions provdided in `cell.proportions` argument.
#' @param scale.factor Vector of scale factors corresponding to alphebetical order of
#' unique type labels in celltype.variable.
#' @param seed.num Random seed for computational reproducibility.
#' @param save.fnstem Character string stem to append to new saved filenames.
#' @param base.path Path of directory to store data for workflow run.
#' @param which.save Which items to save when running this function. Either "sce" for 
#' SingleCellExperiment, "tp" for true proportions, "ypb" for pseudobulk, or "wt"
#' for the `r-nf` compatible workflow table. These objects correspond to filepaths
#' in save.paths (see details).
#' @param save.names List of new file names to save.
#' @returns List of iterations data, saving new files as side-effect.
#' @details Prepares files for a subsampling experiment to evaluate the impact 
#' of a group-wise bias effect on deconvolution outcomes. Either real or simulated
#' data can be provided to the sce argument. This is then randomly subsampled 
#' where remaining arguments specify details for random iterations. 
#' 
#' The fraction of cells per type and group is specified by the fraction.cells 
#' argument. To ensure the same number of cells are sampled across iterations,
#' the final total of cells by type to sample is taken from the group-wise 
#' minima.
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
#' sce[["Sample"]] <- c(rep("sample1", 10), rep("sample2", 2))
#' list.iter <- prepare_subsample_experiment(sce, groups.per.iteration = 1, 
#'                                          scale.factor = c(type1=1,type2=1), 
#'                                          celltype.variable = "celltype", 
#'                                          assay.name = "counts", which.save = c())
#' @export
prepare_subsample_experiment <- function(sce, scale.factor, iterations = 10, 
                                         groups.per.iteration = 3,
                                         method.vector = c("nnls", "music"),
                                         celltype.variable = "k2", 
                                         group.variable = "Sample",
                                         assay.name = "counts_adj",
                                         cell.proportions = c("glial" = 0.3, 
                                                              "neuron" = 0.7),
                                         count.minimum = 200,
                                         seed.num = 0,
                                         save.fnstem = "",
                                         which.save = c("sce", "tp", "ypb", "wt", "li"),
                                         base.path = "data",
                                         save.names = list(sce.name = "sce.rda",
                                                           wt.name = "workflow-table.csv",
                                                           tp.name = "true-proportions.rda",
                                                           ypb.name = "ypb.rda",
                                                           li.name = "lindex.rda"),
                                         verbose = TRUE){
  require(SummarizedExperiment)
  require(SingleCellExperiment)
  set.seed(seed.num)
  
  # get metadata vectors
  cd <- colData(sce)
  groups.vector <- cd[,group.variable]
  celltype.vector <- cd[,celltype.variable]
  # get unique groups
  unique.groups <- unique(groups.vector)
  # get alphabetized celltype labels
  unique.types <- unique(celltype.vector)
  unique.types <- unique.types[order(unique.types)]
  message("Found ",length(unique.types), " cell type labels: ", 
          paste0(unique.types, collapse = ";"))
  # parse scale factors
  S <- scale.factor
  S <- S[names(S) %in% unique.types]
  S <- S[order(match(names(S), names(unique.types)))]
  message("Found ",length(S), " cell types with scale factors")
  
  # prepare iterations parameters
  if(groups.per.iteration > length(unique.groups)){
    groups.per.iteration <- "NULL"}
  # get the exact numbers of cells of each type to sample
  dft <- as.data.frame(table(celltype.vector, groups.vector))
  
  # parse cell amount param
  num.cells.vector <- get_cell_quantities(cell.proportions, count.minimum)
  
  if(verbose){message("Getting random cell indices for iterations...")}
  lindex <- lapply(seq(iterations), function(ii){
    if(!is(groups.per.iteration, "NULL")){
      random.groups <-  sample(unique.groups, groups.per.iteration)
      cdf <- cd[cd[,group.variable] %in% random.groups,]
    } else{
      cdf <- cd
      random.groups = unique.groups
    }
    cell.index.vector <- unlist(lapply(unique.types, function(typei){
      type.filter <- cdf[,celltype.variable]==typei
      cell.id.vector <- rownames(cdf[type.filter,])
      sample(cell.id.vector, num.cells.vector[typei])
    }))
    vindex <- which(colnames(sce) %in% cell.index.vector)
    list(vindex = vindex, samples = random.groups) # return
  })
  
  if(verbose){message("Getting pseudobulk and true proportions...")}
  tp <- as.data.frame(table(sce[[celltype.variable]]))
  tp.prop <- tp[,2]; names(tp.prop) <- tp[,1]
  P <- tp.prop/sum(tp.prop)
  Z <- do.call(cbind, lapply(unique.types, function(typei){
    type.filter <- sce[[celltype.variable]]==typei
    scef <- sce[,type.filter]
    matrix <- as.matrix(assays(scef)[[assay.name]])
    rowMeans(matrix)
  }))
  ZS <- sweep(Z, 2, S, "*")
  ypb <- t(t(P) %*% t(ZS))
  
  if(verbose){message("Parsing save filenames...")}
  if(!save.fnstem==""){
    save.names <- lapply(save.names, function(namei){
      name.str <- unlist(strsplit(namei, "\\."))
      new.name <- paste(name.str[1:length(name.str)-1], save.fnstem, sep = "_")
      new.name <- paste(new.name, name.str[length(name.str)], sep = ".")
      return(new.name)
    })
  }
  
  if(verbose){message("Getting new workflow table...")}
  wti <- data.frame(iterations_index = seq(iterations))
  wti$sample_id <- unlist(lapply(lindex, function(li){
    paste0(li$samples, collapse = ";")}))
  wti$celltype_variable <- celltype.variable
  wti$assay_name <- assay.name
  # manage filepaths
  cnamev <- c("sce_filepath", "bulk_filepath", "list_index_filepath",
              "true_proportions_filepath")
  fpathv <- c(file.path(base.path, save.names[["sce.name"]]), 
              file.path(base.path, save.names[["ypb.name"]]),
              file.path(base.path, save.names[["li.name"]]), 
              file.path(base.path, save.names[["tp.name"]]))
  for(ii in seq(length(cnamev))){
    wti[,cnamev[ii]] <- paste0('"$launchDir/', fpathv[ii], '"')}
  wti$method <- "NA"
  wt <- do.call(rbind, lapply(method.vector, function(methodi){
    wti$method <- methodi; wti
  }))
  
  if(verbose){message("Saving new data...")}
  if(dir.exists(base.path)){
    if("ypb" %in% which.save){
      ypb.fpath <- file.path(base.path, save.names[["ypb.name"]])
      save(ypb, file = ypb.fpath)
    }
    if("tp" %in% which.save){
      tp.fpath <- file.path(base.path, save.names[["tp.name"]])
      save(P, file = tp.fpath)
    }
    if("sce" %in% which.save){
      sce.fpath <- file.path(base.path, save.names[["sce.name"]])
      save(sce, file = sce.fpath)
    }
    if("li" %in% which.save){
      li.fpath <- file.path(base.path, save.names[["li.name"]])
      save(lindex, file = li.fpath)
    }
    if("wt" %in% which.save){
      wt.fpath <- file.path(base.path, save.names[["wt.name"]])
      write.csv(wt, file = wt.fpath, row.names = F)
    }
  } else{
    message("Warning, base.path (",base.path,") not found. Skipping file saves.")
  }
  
  message("Finished run prep. Returning iterations index list.")
  lr <- list(lindex = lindex, wt = wt)
  return(lr)
}

#' get_cell_quantities
#' 
#' Get the number of each cell type to use across subsample iterations
#'
#' @param proportions Vector of target proportions by cell type.
#' @param num.min Known quantity of cells for minimum proportion cell type(s).
#' @returns Vector of cell counts according to specified proportions.
#' @details Gets cell quantities from a provided proportions on minimum count of
#' least-abundant type according to proportions.
#' @examples 
#' get_cell_quantities(c("glial" = 0.3, "neuron" = 0.7), 100)
#' # glial   neuron 
#' # 100.0000 233.3333 
#' 
#' get_cell_quantities(c("glial" = 0.3, "neuron" = 0.7), 200)
#' # glial   neuron 
#' # 200.0000 466.6667 
#' 
#' get_cell_quantities(c("glial" = 0.3, "neuron" = 5), 200)
#' # glial   neuron 
#' # 200.000 3333.333 
#' 
#' get_cell_quantities(c("glial" = 0.3, "neuron" = 5, "other" = 0.1), 200)
#' # glial neuron  other 
#' # 600  10000    200 
#' 
#' @export
get_cell_quantities <- function(proportions = c("glial" = 0.3, "neuron" = 0.7), 
                                num.min = 100){
  if(sum(proportions) > 1){proportions <- proportions/sum(proportions)}
  return(proportions * num.min/min(proportions))
}






#'
#'
#'
#'
analyze_subsample_results <- function(){
  
}

#' subsample_summary
#'
#' Get summary statistics for a results table from a subsampling experiment.
#' 
#' @param results.table
#' @returns dfs, data.frame containing summary statistics.
#' @examples 
#' path <- file.path("inst","examples","results-table.csv")
#' path <- system.file("./inst/examples/results-table.csv", package = 'lute')
#' results.table <- read.csv(path)
#' dfstat <- subsample_summary(results.table)
#' @export
subsample_summary <- function(results.table){
  methodv <- c(unique(dfr$deconvolution_method), "all")
  funv <- c("median", "sd", "length")
  metricv <- c("bias", "rmse.types")
  
  # prepare results data
  dfrs1 <- dfr
  # add type label
  cnv <- colnames(dfrs1)
  unique.types <- unique(unlist(strsplit(dfr$type_labels, ";")))
  unique.types <- unique.types[order(unique.types)]
  for(typei in unique.types){
    cn.filt <- grepl(paste0("type", which(unique.types==typei)), cnv)
    colnames(dfrs1)[cn.filt] <- paste0(colnames(dfrs1)[cn.filt], ".", typei)
  }
  # get all method category
  dfrs2 <- dfrs1; dfrs2$deconvolution_method <- "all"
  dfrs3 <- rbind(dfrs1, dfrs2) # append all category
  methods.vector <- dfrs3$deconvolution_method
  lvar <- list(method = methods.vector)
  unique.methods <- unique(methods.vector)
  unique.methods <- unique.methods[order(unique.methods)]
  # get new colnames for aggregate
  cnv <- colnames(dfrs3)
  grepl.str <- paste0(metricv, collapse = "|")
  cnvf <- cnv[grepl(grepl.str, cnv)]
  
  # get aggregate statistics
  dfs <- do.call(cbind, lapply(funv, function(fi){
    dfai <- aggregate(dfrs3[,cnvf], lvar, FUN = fi); dfai <- dfai[,2:ncol(dfai)]
    fi.str <- fi; if(fi=="length"){fi.str <- "count"}
    colnames(dfai) <- paste0(fi.str, "_", colnames(dfai)); return(dfai)
  }))
  dfs$method <- unique.methods; cnv <- colnames(dfs)
  dfs <- dfs[,c("method", cnv[!cnv=="method"])]
  
  # save
  fname <- "df-sstat-rnf_intra-donor-subsample_ro1-dlpfc.csv"
  write.csv(dfs, file = file.path(save.dpath, fname))
}

#' subsample_plots
#'
#' Plot results of subsample experiments.
#'
#' @param results.table Either a table or a list of such tables, where the list
#' names correspond to experiment labels to use.
#' @returns List of ggplot objects
#' @examples 
#  # load example results table
#' path <- file.path("inst","examples","results-table.csv")
#' path <- system.file("./inst/examples/results-table.csv", package = 'lute')
#' results.table <- read.csv(path)
#'
#' # plot single results table
#' lgg <- subsample_plots(results.table)
#'
#' # plot multiple results tables
#' lgg <- list("expt1" = results.table, "expt2" = results.table)
#' 
#' @export
subsample_plots <- function(results.table){
  lgg <- list() # start return list
  # parse input data
  dfr <- results.table
  is.list <- is(dfr, "list")
  if(is.list){
    dfr <- do.call(rbind, lapply(seq(length(dfr)), function(ii){
      dfi <- results.table[[ii]]; dfi$experiment <- names(results.table)[ii]
      dfi
    }))
  } else{
    dfr$experiment <- "NA"
  }
  unique.methods <- unique(dfr$deconvolution_method)
  unique.types <- unique(unlist(strsplit(dfr$type_labels, ";")))
  
  message("plotting bias")
  metric.plot <- 'bias'
  if(length(unique.methods) > 1){
    message("getting scatter plot of first two methods...")
    # data by method, with unique id
    ldfp <- lapply(unique.methods, function(methodi){
      dfpi <- dfr[dfr$deconvolution_method==methodi,]
      dfpi$iteration.id <- paste0(dfpi$iterations_index, ";", dfpi$experiment)
      dfpi
    })
    id.order.match <- ldfp[[1]]$iteration.id
    ldfp <- lapply(ldfp, function(dfpi){
      dfpi <- dfpi[order(match(dfpi$iteration.id, id.order.match)),]; dfpi
    })
    message("preparing plot data...")
    type.label.vector <- paste0(".type", seq(length(unique.types)))
    names(type.label.vector) <- unique.types
    dfp <- do.call(rbind, lapply(type.label.vector, function(typei){
      var.str <- paste0(metric.plot, typei)
      dfpi <- do.call(cbind, lapply(unique.methods, function(methodi){
        unlist(lapply(ldfp, function(dfpi){
          if(methodi %in% dfpi$deconvolution_method){
            dfpi[dfpi$deconvolution_method==methodi,var.str]
          } else{
            rep("NA", nrow(dfpi))
          }
        }))
      }))
      dfpi <- as.data.frame(dfpi)
      for(c in seq(2)){dfpi[,c] <- as.numeric(dfpi[,c])}
      colnames(dfpi) <- paste0("method", seq(2))
      dfpi$experiment <- ldfp[[1]]$iteration.id
      dfpi$type <- names(typev[typev==typei]); dfpi
    }))
    method1.name <- unique.methods[1]; method2.name <- unique.methods[2]
    dfp$experiment.id <- gsub(".*;", "", dfp$experiment)
    message("getting plot object...")
    ggpt <- ggplot(dfp, aes(x = method1, y = method2)) + 
      geom_abline(slope = 1, intercept = 0, color = "black") +
      geom_hline(yintercept = 0, color = "gray") + 
      geom_vline(xintercept = 0, color = "gray") +
      geom_point(alpha = 0.5) + theme_bw() + geom_smooth() +
      ggtitle("Bias (true - predicted proportions)") +
      xlab(method1.name) + ylab(method2.name)
    
    # get facet
    lgg[["ggpt.facet.bias.methods"]] <- ggpt + facet_wrap(~type+experiment.id)
    
    # get overlays
    ggpt.ol <- ggplot(dfp, aes(x = method1, y = method2, 
                               color = experiment.id, group = experiment.id)) + 
      geom_abline(slope = 1, intercept = 0, color = "black") +
      geom_hline(yintercept = 0, color = "gray") + 
      geom_vline(xintercept = 0, color = "gray") +
      geom_point(alpha = 0.5) + theme_bw() + geom_smooth() +
      ggtitle("Bias (true - predicted proportions)") +
      xlab(method1.name) + ylab(method2.name)
    lgg[["ggpt.overlay.bias.methods"]] <- ggpt.ol
  }
  
  message("getting violin and box plots...")
  bias.label.vector <- paste0("bias.type", seq(length(unique.types)))
  dfp <- do.call(rbind, lapply(unique.methods[1:2], function(methodi){
    dfi <- dfr[dfr$deconvolution_method==methodi,]
    bias.vector <- unlist(lapply(bias.label.vector, function(labeli){dfi[,labeli]}))
    dfpi <- data.frame(bias = bias.vector, type = rep(unique.types, each = nrow(dfi)))
    dfpi$experiment <- rep(dfi$experiment, length(bias.label.vector))
    dfpi$method <- methodi
    dfpi
  }))
  # get facets
  ggvp <- ggplot(dfp, aes(x = experiment, y = bias)) + theme_bw() +
    geom_hline(yintercept = 0, color = "blue") + geom_jitter(alpha = 1e-1) +
    geom_violin(draw_quantiles = 0.5, color = "black", alpha = 0) +
    ylab("Bias (true - predicted proportions)") + xlab("Experiment")
  ggbox <- ggplot(dfp, aes(x = experiment, y = bias)) + theme_bw() +
    geom_hline(yintercept = 0, color = "black") + geom_jitter(alpha = 5e-1) +
    geom_boxplot(color = "blue", alpha = 0) + 
    ylab("Bias (true - predicted proportions)") + xlab("Experiment")
  if(length(unique.methods) > 1){
    ggvp <- ggvp + facet_wrap(~type+method)
    ggbox <- ggbox + facet_wrap(~type+method)
  } else{
    ggvp <- ggvp + facet_wrap(~type)
    ggbox <- ggbox + facet_wrap(~type)
  }
  lgg[["ggviolin.facet.bias.methods"]] <- ggvp
  lgg[["ggboxplot.facet.bias.methods"]] <- ggbox
  message("Finished all plots. Returning...")
  return(lgg)
}
