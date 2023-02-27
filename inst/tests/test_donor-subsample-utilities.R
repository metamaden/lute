require(lute)

sce <- random_sce()
sce[["Sample"]] <- c(rep("sample1", 10), rep("sample2", 2))
list.iter <- prepare_subsample_experiment(sce, groups.per.iteration = 1,
                                          scale.factor = c(type1=1,type2=1),
                                          celltype.variable = "celltype",
                                          assay.name = "counts", which.save = c())

path <- file.path("inst","examples","results-table.csv")
path <- file.path(system.file(package = 'lute'), path)
results.table <- read.csv(path)
dfstat <- subsample_summary(results.table)
lgg <- list("expt1" = results.table, "expt2" = results.table)