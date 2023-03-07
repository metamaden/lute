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
lexample <- .get_decon_example_data()
z <- lexample[["z"]]
y <- lexample[["y"]]
s <- lexample[["s"]]
names(s) <- colnames(z)
p <- c(0.3, 0.7)

zs <- sweep(z, 2, s, FUN = "*")
ys <- t(t(p) %*% t(zs))

prop.adj <- deconvolution(epicParam(s = s, y = ys, z = zs, return.info = T))
prop.unadj <- deconvolution(epicParam(s = s, y = ys, z = z, return.info = T))

bias.adj <- deconvolution(epicParam(s = s, y = ys, z = zs))[1:2]-p
bias.unadj <- deconvolution(epicParam(s = s, y = ys, z = z))[1:2]-p





