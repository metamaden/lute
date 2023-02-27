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
Yi <- cbind(Y, Y)
colnames(Yi) <- paste0("sample", seq(2))
res <- DeconRNASeq(as.data.frame(Yi), as.data.frame(Z), use.scale = FALSE)

source("DeconRNASeq.R")

data(multi_tissue)
datasets <- x.data[,2:11]
signatures <- x.signature.filtered.optimal[,2:6]
proportions <- fraction
DeconRNASeq(datasets, signatures, proportions, 
            checksig=FALSE, known.prop = TRUE, use.scale = TRUE)
