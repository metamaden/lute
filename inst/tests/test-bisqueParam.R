library(lute)
source("./R/lute_utilities.R")


#----------------------------------
# make the basic class from scratch
#----------------------------------
lexample <- .get_decon_example_data()

lexample <- .get_decon_example_data_bisque()

# get sc eset
sce <- random_sce(num.genes = 10, num.cells = 100, num.types = 2)
df.z.pheno <- data.frame(cellType = sce[["celltype"]], 
                         SubjectName = paste0("sample", seq(ncol(sce))))
rownames(df.z.pheno) <- colnames(sce)
z.eset <- ExpressionSet(assayData = counts(sce),
                        phenoData = AnnotatedDataFrame(df.z.pheno))
rownames(z.eset) <- rownames(y.eset)

# make eset with multiple samples from y
y <- lparam[["y"]]
y <- cbind(y, y, y, y, y, y)
colnames(y) <- c(paste0("sample", seq(2)), paste0("bulk",seq(4)))
df.y.pheno <- data.frame(SubjectName = colnames(y))
rownames(df.y.pheno) <- colnames(y)
y.eset <- ExpressionSet(assayData = y,
                        phenoData = AnnotatedDataFrame(df.y.pheno))

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
