#!/usr/bin/env R

# Author: Sean Maden
#
# Methods for SummarizedExperimentType objects. 
#
# Sections:
#
#   1. conversions: Between-object class conversions.
#
#   2. methods for adjustments: Adjustment and weighting methods.
#
#   3. set simulations: Methods supporting simultions using 
#       SummarizedExperimentTypes objects.
#

#---------------
# 1. conversions
#---------------
# convert between two classes or two subclasses.

#' set_from_sce
#'
#' Make an object of class `SummarizedExperimentTypes` from an object of class
#' `SingleCellExperiment` or `SummarizedExperiment`.
#'
#' @param sce A `SingleCellExperiment` object.
#' @param method Method to summarize assays on types.
#' @param type.variable Variable containing the type-level labels.
#' @param assayname Name of assays data in assays(sce).
#' @param verbose Whether to show verbose status messages.
#' groups (see `?sce_groupstat` for details).
#' @returns `SummarizedExperimentTypes` object.
#' @examples
#' set <- set_from_sce(random_sce()); show(set)
#' @seealso sce_groupstat
#' @export
set_from_sce <- function(sce, method = "mean", type.variable = "celltype", 
                         assayname = "counts", verbose = FALSE){
  # run checks
  if(!(is(sce, "SingleCellExperiment")|is(sce, "SummarizedExperiment"))){
    stop("sce must be of class SingleCellExperiment or SummarizedExperiment.")}
  
  # get new assays matrix
  typev <- unique(sce[[type.variable]])
  ma <- do.call(cbind, lapply(typev, function(typei){
    if(verbose){message("Summarizing type: ", typei, "...")}
    type.filt <- sce[[type.variable]]==typei
    scef <- sce[,type.filt]
    exprf <- assays(scef)[[assayname]]
    mai <- matrix(rowMeans(exprf), ncol = 1)
    colnames(mai) <- typei
    return(mai)
  }))
  rownames(ma) <- rownames(sce)
  lma <- list(ma); names(lma) <- paste0("summarized_", assayname)
  
  # parse metadata
  # get rowdata
  rd <- sce_groupstat(sce, group.variable = type.variable, 
                      summarytype = "rowData", return.tall = FALSE)
  marker.filt <- grepl(".*;marker$", colnames(rd))
  rownames(rd) <- rd[,which(marker.filt)[1]]
  rd$marker <- rd[,which(marker.filt)[1]]
  rd <- rd[,!marker.filt]
  rd <- rd[order(match(rd$marker, rownames(ma))),]
  cond.rd <- identical(rd$marker, rownames(ma))
  if(!cond.rd){
    message("Warning, couldn't match rowdata markers to ma markers.")
    rd <- matrix(nrow = nrow(ma), ncol = 0)
  }
  # get coldata
  cd <- sce_groupstat(sce, group.variable = type.variable, 
                      summarytype = "colData", return.tall = TRUE)
  rownames(cd) <- cd$group; cd$type <- cd$group
  cd <- cd[,!colnames(cd)=="group"]
  cd <- cd[order(match(cd$type, colnames(ma))),]
  cond.cd <- identical(cd$type, colnames(ma))
  if(!cond.cd){
    message("Warning, couldn't match coldata types to ma types.")
    cd <- matrix(nrow = ncol(ma), ncol = 0)
  }
  
  # metadata
  lmd <- list(assay.info = list(
    stat.method = method, 
    sce.assayname = assayname, 
    type.variable = type.variable
  ))
  
  # make new set object
  set <- SummarizedExperimentTypes(assays = lma, rowData = rd, 
                                   colData = cd, metadata = lmd)
  
  return(set)
}

#' set_from_set
#'
#' Expects a set object with type;group nesting. Returns a set object summarized
#' by each type across any available groups.
#'
#' @param set An object of type SummarizedExperimentTypes or similar.
#' @param group.variable Name of variable containing group level information.
#' @param type.variable Name of variable containing type level information.
#' @param method Statistical method to summarize types on group levels for new
#' assays. Can be either of "mean" or "median".
#' @param assayname Name of assays data in provided set object to summarize.
#' @param make.set.plots Whether to make standard set plots (e.g. heatmap, PCA).
#' @param verbose Whether to show verbose status messages.
#' @returns New set object with type data summarized on provided groups.
#' @examples 
#' sce <- random_sce()
#' sce[["donor"]] <- c(rep("donor1", 2), rep("donor2", 8))
#' sce[["type.variable"]] <- paste0(sce[["celltype"]], ";", sce[["donor"]])
#' set <- set_from_sce(sce, group.variable = "donor", type.variable = "celltype")
#' set2 <- set_from_set(set, type.variable = "celltype", group.variable = "donor")
#' @export
set_from_set <- function(set, group.variable = "donor", type.variable = "celltype", 
                         method = "mean", assayname = "summarized_counts",  
                         make.set.plots = TRUE, verbose = FALSE, ...){
  # check set class
  setclass.cond <- (is(set, "SummarizedExperimentTypes")|
                      is(set, "RangedSummarizedExperimentTypes"))
  if(!setclass.cond){
    stop("set must be of class SummarizedExperimentTypes or similar.")}
  
  # get params
  typev <- unique(sce[[type.variable]])
  expr.sce <- assays(sce)$counts
  if(!is(group.variable, "NULL")){
    ugroupv <- unique(sce[[group.variable]]) # get all possible group lvls
  }
  
  # make new assays data
  cnv <- colnames(set)
  ldat <- lapply(typev, function(typei){
    # filter set data
    type.filt <- grepl(paste0("^", typei, "(;.*|$)"), cnv)
    mei <- assays(set[,type.filt,drop = F])[[assayname]]
    mei <- as.matrix(mei)
    # get data summary variables
    group.var <- group.mean <- NA
    unique.group.lvl <- unique(gsub(paste0("^",typei,";"), "", colnames(mei)))
    num.group <- length(unique.group.lvl)
    unique.group.lvl <- paste0(unique.group.lvl, collapse = ";")
    if(num.group > 1){group.var <- rowVars(mei);group.mean <- rowMeans(mei)}
    # get summarized expr data
    mei <- make_new_assaydata(mei, method = method, 
                              na.rm = T, verbose = verbose)
    # return as list
    list(rd = data.frame(group.var = group.var,
                         group.mean = group.mean),
         cd = data.frame(type = typei, num.group = num.group, 
                         group.lvl = unique.group.lvl),
         mexpr = mei)
  })
  names(ldat) <- typev
  
  # assemble set components
  # assay expr table
  new.assay <- do.call(cbind, lapply(ldat, function(li){li$mexpr}))
  colnames(new.assay) <- names(ldat)
  # coldata
  new.cd <- do.call(rbind, lapply(ldat, function(li){li$cd}))
  rownames(new.cd) <- colnames(new.assay)
  # rowdata
  new.rd <- do.call(cbind, lapply(seq(length(ldat)), function(ii){
    typei <- names(ldat)[ii]; rdi <- ldat[[ii]]$rd
    colnames(rdi) <- paste0(typei, ";", colnames(rdi))
    rdi
  }))
  rownames(new.rd) <- rownames(set)
  # metadata
  lmd <- list(assay.info = list(stat.method = method,
                                   type.variable = type.variable),
                 set.original = set)
  # make new set object
  la <- list(assayname = new.assay)
  names(la) <- paste0("summarized_", gsub(".*_", "", assayname))
  set.new <- SummarizedExperimentTypes(assays = la, rowData = new.rd, 
                                       colData = new.cd, metadata = lmd)
  # parse standard plot options
  if(make.set.plots){
     lp <- get_set_plots(set = set.new, group.variable = group.variable,
                         type.variable = type.variable, verbose = verbose, ...)
     metadata(set.new)[["set_plots"]] <- lp
  }
  return(set.new)
}

#' make_new_assaydata
#'
#' Efficiently makes the new assay data by summarizing sce assays with context
#' sensitivity. Currently supports summary functions using either 
#' DelayedMartrixStats (if a `DelayedArray` object provided) and Matrix (all 
#' other inputs).
#'
#' @param ma Assay data matrix (cols = samples/cells/donors, rows = genes/loci).
#' Accepts tables of type data.frame, matrix, or DelayedArray.
#' @param method Statistical method to perform summary.
#' @param na.rm Option for NA removal.
#' @param verbose Whether to show verbose status messages.
#' @returns Assay data matrix.
#' @export
make_new_assaydata <- function(ma, method = "mean", na.rm = TRUE, verbose = FALSE){
  if(verbose){message("Parsing ma class...")}
  ma.isda <- is(ma, "DelayedArray")
  ma.ismat <- is(ma, "matrix") & !ma.isda
  ma.isdf <- is(ma, "data.frame") & !ma.isda
  if(!(ma.isda|ma.ismat|ma.isdf)){stop("Error: invalid ma class.")}
  if(verbose){message("Getting ma in valid class...")}
  ma.form <- ma
  if(ma.isdf){ma.form <- as.matrix(ma)}
  ma.form.isda <- is(ma.form, "DelayedArray")
  ma.form.ismat <- is(ma.form, "matrix") & !ma.isda
  if(verbose)(message("Setting the R library..."))
  lib.str <- ifelse(ma.form.isda, "DelayedMatrixStats", "Matrix")
  if(verbose)(message("Setting the summary method..."))
  method.str <- "NA"
  if(method == "mean"){
    method.str <- "rowMeans"
  } else if(method == "median"){
    method.str <- "rowMedians"
  } else{
    stop("Error: invalid summary method.")
  }
  if(verbose){message("Getting summary assay data...")}
  call.str <- paste0(lib.str, "::", method.str,"(ma, na.rm = ",na.rm,")")
  ma.new <- eval(parse(text = call.str))
  return(ma.new)
}

#---------------------------
# 3. methods for adjustments
#---------------------------

#' append_gvariable_adj
#'
#' Perform group-wise adjustments on a SummarizedExperimentTypes object. Appends 
#' group adjustments to rowData and group-adjusted signals to assays.
#'
#' @param set A SummarizedExperimentTypes object
#' @param gvariable_pattern String pattern to identify group variables from 
#' rowData, if present.
#' @param assayname Name of original assays data to adjust, as identifiable from
#' names(assays(set)).
#' @param rd.method.str Character string to append for new assays and rowData 
#' matrices.
#' @param type Adjustment type (options: "mgvdenom").
#' @param verbose Whether to show verbose status messages.
#' @param ... Additional arguments passed for specific adjustment type function.
#' @return set object with adjusted marker signal assay data appended.
#' @examples
#' sce <- random_sce()
#' sce[["donor"]] <- c(rep("donor1", 2), rep("donor2", 8))
#' set <- set_from_sce(sce, group.variable = "donor")
#' set <- append_gvariable_adj(set)
#' @export
append_gvariable_adj <- function(set, gvariable_pattern = "donor.*", 
                            assayname = "counts", type = "mgvdenom",
                            verbose = FALSE, ...){
  if(!(is(set, "SummarizedExperimentTypes"))){
    stop("set must be of class SummarizedExperimentTypes.")}
  if(type == "mgvdenom"){
    set <- groupadj_mgvdenom_fromrd(set = set, 
                                    gvariable_pattern = gvariable_pattern,
                                    rd.method.str = type, assayname = assayname,
                                    verbose = verbose, ...)
  } else{stop("Error: invalid adjustment type.")}
  return(set)
}

#' groupadj_mgvdenom_fromrd
#'
#' Get adjusted assays and rowadata, calculated as the counts over the mean 
#' group-wise variances ("mgvdenom"), or: 
#' 
#' $$counts/mgvdenom$$
#' 
#' @param set A SummarizedExperimentTypes
#' @param gvariable_pattern String pattern to identify group-level summaries from
#' rowData column labels with grepl.
#' @param assayname Name of original assays data to adjust, as identifiable from
#' names(assays(set)).
#' @param rd.method.str Character string to append for new assays and rowData 
#' matrices.
#' @param verbose Whether to return verbose status messages.
#' @returns set with updated rowData columns and new assays matrix.
#' @export
groupadj_mgvdenom_fromrd <- function(set, 
                                     gvariable_pattern = "donor.*",
                                     assayname = "counts",
                                     rd.method.str = "mgvdenom", 
                                     verbose = FALSE){
  if(verbose){message("Checking for groups in rowData...")}
  typev <- colnames(set); rd <- rowData(set); rd.cnv <- colnames(rd)
  str.patt <- paste0("^type[0-9];", gvariable_pattern, ";var$")
  rd.cnvf <- rd.cnv[grepl(str.patt, rd.cnv)]
  if(length(rd.cnvf)==0){stop("Error, no groups found in rowData.")}
  if(verbose){message("Group data found. Getting marker adjustments...")}
  rd.new <- do.call(cbind, lapply(typev, function(typei){
    type.str.patt <- paste0("^", typei, ";.*")
    which.rd.filt <- which(grepl(type.str.patt, rd.cnvf))
    rd.cnvf.filt <- rd.cnvf[which.rd.filt]
    rdf <- as.data.frame(rd[,rd.cnvf.filt])
    rowMeans(rdf, na.rm = T)
  }))
  colnames(rd.new) <- paste0(typev, ";", rd.method.str)
  if(verbose){message("Updating rowData...")}
  rowData(set) <- cbind(rd, rd.new) # update rowdata
  if(verbose){message("Getting new adjusted assay data...")}
  new.assay.name <- paste0(assayname,"_adj_",rd.method.str)
  assays(set)[[new.assay.name]] <- assays(set)[[assayname]]/rd.new
  return(set)
}

#-------------------
# 4. set simulations
#-------------------

#-------------
# 5. set plots
#-------------
# Make standard plots for set objects

#' get_set_plots
#'
#' Makes standard plots for SummarizedExperimentTypes and similar objects.
#'
#' @param set A SummarizedExperimentTypes object or similar.
#' @param lplots List of plot types to make. Accepts "pca" and "hm".
#' @param assayname Name of assays matrix in set object to plot.
#' @param type.variable Variable name for type labels in set colData.
#' @param group.variable Variable name for type labels in set colData.
#' @param mtype.variable Variable name for marker labels in rowData.
#' @param randcol.seednum Number for random seed to make random colors.
#' @param hm.topanno Optional object produced using HeatmapAnnotation(). If NULL,
#' makes this annotation using set colData and other specified arguments.
#' @param hm.leftanno Optional object produced using rowAnnotation(). If NULL,
#' makes this annotatiomn using set rowData and other specified arguments.
#' @param scale.hmdata Whether to rescale heatmap data with scale().
#' @param verbose Whether to show verbose status messages.
#' @returns List of plot objects
#' @export
get_set_plots <- function(set, lplots = c("hm"), 
                          assayname = "summarized_counts", 
                          type.variable = NULL, group.variable = NULL,
                          mtype.variable = NULL, randcol.seednum = 0, 
                          hm.topanno = NULL, hm.leftanno = NULL,
                          scale.hmdata = TRUE, verbose, ...){
  lp <- list()
  if("hm" %in% lplots|"heatmap" %in% lplots){
    if(verbose){message("Making new heatmap...")}
    lp[["heatmap"]] <- get_set_heatmap(set = set, 
                                       assayname = assayname,
                                       type.variable = type.variable, 
                                       group.variable = group.variable,
                                       mtype.variable = mtype.variable,
                                       randcol.seednum = randcol.seednum,
                                       hm.topanno = hm.topanno,
                                       hm.leftanno = hm.leftanno,
                                       scale.hmdata = scale.hmdata,
                                       verbose = verbose)
  }
  if("pca" %in% lplots|"principalcomponentanalysis" %in% lplots){
    if(verbose){message("Making new PCA...")}
  }
  return(lp)
}

#' get_set_heatmap
#'
#' Make a standard heatmap of set marker expression with row (marker) and 
#' column (type) annotations.
#' 
#' @param set A SummarizedExperimentTypes object or similar.
#' @param hm.topanno Optional object produced using HeatmapAnnotation(). If NULL,
#' makes this annotation using set colData and other specified arguments.
#' @param hm.leftanno Optional object produced using rowAnnotation(). If NULL,
#' makes this annotation using set rowData and other specified arguments.
#' @param assayname Name of assays matrix in set object to plot.
#' @param type.variable Variable name for type labels in set colData.
#' @param group.variable Variable name for type labels in set colData.
#' @param mtype.variable Variable name for marker labels in rowData.
#' @param randcol.seednum Number for random seed to make random colors.
#' @param scale.hmdata Whether to rescale heatmap data with scale().
#' @param verbose Whether to show verbose status messages.
#' @returns Returns heatmap object returned from ComplexHeatmap::Heatmap()
#' @examples 
#' sce <- random_sce()
#' sce[["donor"]] <- c(rep("donor1", 2), rep("donor2", 8))
#' sce[["typevar"]] <- paste0(sce[["celltype"]], ";", sce[["donor"]])
#' set <- set_from_sce(sce, group.variable = "donor", type.variable = "typevar")
#' metadata(set)[["set_plots"]]$heatmap
#' @export
get_set_heatmap <- function(set, assayname = "logcounts_bytype",
                            type.variable = NULL, group.variable = NULL, 
                            mtype.variable = NULL, randcol.seednum = 0, 
                            scale.hmdata = TRUE, hm.topanno = NULL, 
                            hm.leftanno = NULL, verbose = FALSE){
  suppressPackageStartupMessages(require(ComplexHeatmap))
  if(!is(set,  "SummarizedExperimentTypes")){
    stop(paste0("Error: set must be an object of class ",
                "SummarizedExperimentTypes or similar."))
  }
  # parse heatmap assays data
  if(!assayname %in% names(assays(set))){
    if(verbose){message("Warning: assayname '",assayname,
                        "' not found in set assays.")}
    if(verbose){message("Checking available assays...")}
    assayname <- names(assays(set))[1]
    if(is(assayname, "NULL")){stop("Error: no assay names in set.")} else{
      if(verbose){message("Using assay '",assayname,"'...")}
    }
  }
  hm.data <- assays(set)[[assayname]]
  # parse top annotation options
  if(is(hm.topanno, "NULL")){
    if(verbose){message("Making new top annotation from arguments...")}
    set.seed(randcol.seednum); cd <- colData(set); topanno.str <- NULL
    if(verbose){message("Checking variables...")}
    if(is(type.variable, "NULL")){
      if(verbose){message("Taking assay colnames as type variable.")}
      type.variable <- "type"; set[[type.variable]] <- colnames(set)
    }
    if(is(group.variable, "NULL")){
      if(verbose){message("Proceeding without specifying group variable...")}
      topanno.str <- paste0("HeatmapAnnotation(", 
                            type.variable, " = set[[type.variable]],", 
                            "annotation_name_side = 'left')")
    } else{
      if(group.variable %in% colnames(cd)){
        if(verbose){message("Proceeding with specified group variable...")}
        topanno.str <- paste0("HeatmapAnnotation(", 
                              type.variable, " = set[[type.variable]],",
                              group.variable," = set[[group.variable]], ",
                              "annotation_name_side = 'left')")
      } else{
        if(verbose){message("Warning: group variable '", 
                            group.variable,"' not found in colData.")}
      }
    }
    topanno <- eval(parse(text = paste0(topanno.str))) # parse string as command
  } else{
    if(verbose){message("Using provided top annotation...")}
    topanno <- hm.topanno
  }
  # parse left anno
  if(is(hm.leftanno, "NULL")){
    if(verbose){message("Making left annotation from arguments...")}
    leftanno <- NULL
    if(!is(mtype.variable, "NULL")){
      rd <- rowData(set)
      if(mtype.variable %in% colnames(rd)){
        if(verbose){message("Getting marker type variable from rowData...")}
        leftanno <- rowAnnotation(marker_type = rd[,mtype.variable])
      } else{
        if(verbose){message("Warning: marker type variable not in rowData.")}
      }
    }
  } else{
    if(verbose){message("Using provided left annotation...")}
    leftanno <- hm.leftanno
  }
  # get legend character string
  # parse scale option
  legend.str <- gsub("_.*", "", assayname)
  if(scale.hmdata == TRUE){
    if(verbose){message("Scaling heatmap data...")}
    hm.data <- scale(hm.data); legend.str <- paste0("scaled\n", legend.str)
  }
  # parse assays summary metadata
  if(verbose){message("Formatting legend/heatmap name...")}
  if("assay.info" %in% names(metadata(set))){
    ai <- metadata(set)$assay.info  
    if("stat.method" %in% names(ai)){
      legend.str <- paste0(ai[["stat.method"]], "\n", legend.str)
    }
  }
  # make first char uppercase
  lv <- unlist(strsplit(legend.str, "")); lv[1] <- toupper(lv[1])
  legend.final <- paste0(lv, collapse = "")
  if(verbose){message("Making new heatmap object...")}
  Heatmap(hm.data, name = legend.final, show_column_dend = FALSE, 
          top_annotation = topanno, left_annotation = leftanno)
}

#' get_set_pca
#'
#' @param set A SummarizedExperimentTypes object or similar.
#' @param assayname Name of assays matrix in set object to plot.
#' @param type.variable Variable name for type labels in set colData.
#' @param group.variable Variable name for type labels in set colData.
#' @param verbose Whether to show verbose status messages.
#' @returns Returns PCA results with plots.
#' @examples 
#' sce <- random_sce()
#' sce[["donor"]] <- c(rep("donor1", 2), rep("donor2", 8))
#' sce[["typevar"]] <- paste0(sce[["celltype"]], ";", sce[["donor"]])
#' set <- set_from_sce(sce, group.variable = "donor", type.variable = "typevar")
#' get_set_pca(set, assayname = "summarized_counts", type.variable = "type")
#' @export
get_set_pca <- function(set, assayname = "logcounts_bytype",
                        type.variable = NULL, group.variable = NULL){
  
}
