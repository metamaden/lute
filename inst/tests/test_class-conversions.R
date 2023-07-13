#!/usr/bin/env R

#
# testing data class conversion functions
#

# functions to get obtain test data:
# lexample <- lute:::.get_decon_example_data_bisque()
# lexample <- lute:::.get_decon_example_data()
# lexample <- lute:::.get_decon_example_data_music2()

#-----------------
# eset conversions
#-----------------
# get example eset data
lexample <- lute:::.get_decon_example_data_bisque()
sc.eset <- lexample$sc.eset
names(attributes(sc.eset))

# convert eset to se
assay.name <- "counts"
sce.new <- SingleCellExperiment(assays = list(assay.name = exprs(sc.eset)))
colData(sce.new) <- DataFrame(as.matrix(pData(sc.eset)))
metadata(sce.new) <- metadata(sc.eset)

# convert eset to sce

#----------------
# sce conversions
#----------------
# lexample <- lute:::.get_decon_example_data()
lexample <- lute:::.get_decon_example_data_music2()