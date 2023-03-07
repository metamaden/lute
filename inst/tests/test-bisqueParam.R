library(lute)
source("./R/lute_utilities.R")


#----------------------------------
# make the basic class from scratch
#----------------------------------
lexample <- .get_decon_example_data()

lexample <- .get_decon_example_data_bisque()

results <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset = y.eset, 
                                                  sc.eset = z.eset)

proportions <- results$bulk.props
proportions.sample1 <- proportions$bulk.props[,1]

#-------------------------------------
# example from BisqueRNA documentation
#-------------------------------------
set.seed(0)
cell.types <- c("Neurons", "Astrocytes", "Oligodendrocytes", "Microglia", "Endothelial Cells")
avg.props <- c(.5, .2, .2, .07, .03)
sim.data <- SimulateData(n.ind=10, n.genes=100, n.cells=500, cell.types=cell.types, avg.props=avg.props)
sc.eset <- sim.data$sc.eset[,sim.data$sc.eset$SubjectName %in% as.character(6:10)]
bulk.eset <- sim.data$bulk.eset
true.props <- sim.data$props
markers <- sim.data$markers
res <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset, sc.eset, markers=NULL, use.overlap=TRUE)
