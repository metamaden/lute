library(lute)
source("./R/lute_utilities.R")


#----------------------------------
# make the basic class from scratch
#----------------------------------
lexample <- .get_decon_example_data()

lparam <- new("deconParam", s = lexample[["s"]], 
             y = lexample[["y"]], z = lexample[["z"]])

# instantiate objects
y <- lparam[["y"]]
z <- lparam[["z"]]
s <- lparam[["s"]]
# format objects
y <- as.matrix(y)
z <- as.matrix(z)
s <- as.numeric(s)

y.eset <- ExpressionSet(assayData = y)
z.eset <- ExpressionSet()

BisqueRNA::ReferenceBasedDecomposition(bulk.eset = y, sc.eset = z)