#!/usr/bin/env R

# Author: Sean Maden
#
# Test donor simulations.
#

require(lute)

df <- rand_donor_marker_table(ndonor = 10, gindexv = c(1, 2), 
                              sd.offset.pos = 5,sd.offset.neg = 5)


# adjust for donor bias
# use combat
require(sva)
# make expr matrix
filt.donor <- grepl("donor\\d", colnames(df))
mexpr <- do.call(rbind, lapply(unique(df$marker), function(mi){
  dff <- df[df$marker==mi, ]
  unlist(lapply(unique(dff[dff$marker==mi,]$type), function(ti){
    datv <- dff[dff$type==ti, filt.donor]
    names(datv) <- paste0(colnames(dff[,filt.donor]), ";", ti)
    return(datv)
  }))
}))
rownames(mexpr) <- unique(df$marker)
# make pheno 
cnv <- colnames(mexpr)
pheno <- data.frame(donor = gsub(";.*", "", cnv),
                    type = gsub(".*;", "", cnv))
# get combat vars
mod <- model.matrix(~1, data = pheno)
batch <- pheno$donor
madj <- ComBat(dat = mexpr, batch = batch, mod = mod,
               par.prior = TRUE, prior.plots = FALSE)
