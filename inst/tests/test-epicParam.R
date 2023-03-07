source("./R/lute_generics.R")
source("./R/lute_utilities.R")
source("./R/deconParam-class.R")
source("./R/epicParam-class.R")

# example
lexample <- .get_decon_example_data()
param <- epicParam(s = lexample[["s"]], y = lexample[["y"]], z = lexample[["z"]])

# return only predicted proportions
deconvolution(param)
# type1     type2 
# 0.9819837 0.0180163 

# return full results
param@return.info <- T
names(deconvolution(param))
# [1] "predictions" "result.info" "metadata"

#-----------------
# test s transform
#-----------------
library(lute)

source("./R/lute_utilities.R")

lexample <- .get_decon_example_data()
z <- lexample[["z"]]
y <- lexample[["y"]]
# s <- lexample[["s"]]
s <- c("type1" = 3, "type2" = 1000)

names(s) <- colnames(z)
p <- c(0.3, 0.7)

zs <- sweep(z, 2, s, FUN = "*")
ys <- t(t(p) %*% t(zs))

# get random z.var
z.var <- matrix(rnbinom(nrow(z)*ncol(z), size = 10, mu = 10),
                nrow = nrow(z))
rownames(z.var) <- rownames(z)
colnames(z.var) <- colnames(z)
prop.adj.var <- deconvolution(epicParam(s = s, y = ys, z = zs, z.var = z.var, return.info = T))
prop.adj.novar <- deconvolution(epicParam(s = s, y = ys, z = zs, z.var = NULL, return.info = T))
prop.adj.var$predictions
prop.adj.novar$predictions

prop.unadj <- deconvolution(epicParam(s = NULL, y = ys, z = z, return.info = T))

prop.adj$predictions
prop.unadj$predictions


identical(prop.adj$metadata$z, prop.unadj$metadata$z)

bias.adj <- deconvolution(epicParam(s = s, y = ys, z = zs))[1:2]-p
bias.unadj <- deconvolution(epicParam(s = s, y = ys, z = z))[1:2]-p





