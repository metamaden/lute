# get example bulk data
y <- lute:::.get_decon_example_data()$y

# get example sce
sce <- random_sce()[seq(10),]

# get framework results
experiment.results <- lute(sce = sce, y = y)