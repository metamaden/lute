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
Yi <- as.data.frame(cbind(Y, Y))
colnames(Yi) <- paste0("sample", seq(ncol(Yi)))

res <- run_deconvolution(Z = as.data.frame(Z), Y = Yi, method = method)





