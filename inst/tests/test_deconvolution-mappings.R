#!/usr/bin/env R

# Test mappings to deconvolution functions in lute.

require(lute)

sce <- random_sce()
typev <- unique(sce[["celltype"]])
Z <- do.call(cbind, lapply(typev, function(typei){
rowMeans(counts(sce[,sce[["celltype"]]==typei]))
}))
colnames(Z) <- c('type1', 'type2')
Y <- matrix(rowMeans(counts(sce)), ncol = 1)

# run nnls
ldecon <- run_deconvolution(method = "nnls", Y = Y, Z = Z)
# inspect results
ldecon

# run music
method <- "music"
arguments <- list("sigma" = NULL)
sigma <- matrix(0, ncol = 1, nrow = nrow(Z))
arguments <- list("sigma = NULL")
command.list <- map_music(arguments)

# run deconrnaseq
method <- "deconrnaseq"
Z <- as.data.frame(Z)
res <- run_deconvolution(Z = Z, Y = Y, method = method)

# run epic
#devtools::install_github("GfellerLab/EPIC")
library(EPIC)

method <- "epic"


res1 <- EPIC(melanoma_data$counts)
res1$cellFractions
res2 <- EPIC(melanoma_data$counts, TRef)
res3 <- EPIC(bulk=melanoma_data$counts, reference=TRef)
res4 <- EPIC(melanoma_data$counts, reference="TRef")
res5 <- EPIC(melanoma_data$counts, mRNA_cell_sub=c(Bcells=1, otherCells=5))

