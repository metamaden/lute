#!/usr/bin/env R

# Author: Sean Maden
#
# Saves the example data used in the vignette "lute_pseudobulk_example.Rmd".
#
#
#

libv <- c("scRNAseq", "lute")
sapply(libv, library, character.only = TRUE)
data <- DarmanisBrainData()

sample.id.variable <- "experiment_sample_name"
old.types <- "cell.type"; new.types <- "k2"

# remove non-k2 types
filter.k2 <- data[[old.types]] %in% 
  c("neurons", "oligodendrocytes", "astrocytes", "OPC", "microglia")
data <- data[,filter.k2]
# define new k2 variable
data[[new.types]] <- ifelse(data[[old.types]]=="neurons", "neuron", "glial")
data[[new.types]] <- factor(data[[new.types]])

# get markers
markers <- lute(sce = data, celltype.variable = "k2", markers.per.type = 100,
                deconvolution.algorithm = NULL)

# subset on markers
data.example <- data[markers$typemarker.results,]
dim(data.example)

# save
save(data.example, file = "./inst/extdata/scRNAseq/darmanis_example.rda")