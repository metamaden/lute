#!/usr/bin/env R

#
# testing data class conversion functions
#

# functions to get obtain test data:
# lexample <- lute:::.get_decon_example_data_bisque()
# lexample <- lute:::.get_decon_example_data()
# lexample <- lute:::.get_decon_example_data_music2()
# sce <- random_sce()

#-----------------
# eset conversions
#-----------------
# get example eset data
lexample <- lute:::.get_decon_example_data_bisque()
sc.eset <- lexample$sc.eset
names(attributes(sc.eset))

# convert eset to se
assay.name <- "counts"
se.new <- SummarizedExperiment(assays = list(assay.name = exprs(sc.eset)))
colData(se.new) <- DataFrame(as.matrix(pData(sc.eset)))

# convert eset to sce
sce.new <- SingleCellExperiment(assays = list(assay.name = exprs(eset)))
colData(sce.new) <- DataFrame(as.matrix(pData(eset)))
metadata(sce.new) <- metadata(eset)

#----------------
# sce conversions
#----------------
sce <- random_sce()

# convert to eset
assay.name <- "counts"
eset <- ExpressionSet(assayData = assays(sce)[[assay.name]])
pData(eset) <- as.data.frame(colData(sce))
metadata(eset) <- metadata(sce)

# convert to se
se <- SummarizedExperiment(assays = assays(sce))
colData(se) <- colData(sce)
metadata(se) <- metadata(sce)

#---------------
# se conversions
#---------------
# convert to sce
sce <- SingleCellExperiment(assays = assays(se))
colData(sce) <- colData(se)
metadata(sce) <- metadata(se)


# convert to eset
assay.name <- "counts"
eset <- ExpressionSet(assayData = assays(se)[[assay.name]])
pData(eset) <- as.data.frame(colData(se))

