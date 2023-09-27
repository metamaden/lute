#!/usr/bin/R

# Author: Sean Maden
#
#

require(lute)


#----------------
# test biasexpt()
#----------------
donordf = ddf # NULL
method = "nbinom"
lambda.pos = 20
lambda.neg = 2
lambda.sdoff.pos = 5
lambda.sdoff.neg = 2
gamma.pos = 20
gamma.neg = 2
P = c(0.25, 0.75)
donor.adj.method = "combat"
plot.biasadj = TRUE
plot.pca = TRUE
cname.donorsummary = "donor.mean"
gindexv = c(1, 2)
ndonor = 10
seed.num = 0
verbose = TRUE


# get ypb
if(!is(donordf, "NULL")){
  ddff <- ddf[!duplicated(ddf$marker),]
  ddff$type <- as.factor(ddff$type)
  gindexv <- unlist(lapply(levels(ddff$type), function(ti){
    rep(ti, length(which(ddff$type==ty)))
  }))
}

if(!nrow(Ypb) == length(gindexv)){
  stop("Error, couldn't make gindexv from Ypb.")}

df <- random_donordf(ndonor = 1, gindexv = gindexv,
                     lambda.sdoff.pos = 0, lambda.sdoff.neg = 0)
ktotal <- length(P)
Z <- matrix(df[,"donor1"], ncol = ktotal)
Ypb <- ypb_fromtypes(Z = Z, P = P)

# params
df = ddf
Ypb = Ypb 
P = P
donor.unadj = NULL
donor.adj.method = "combat"
plot.biasadj = TRUE
verbose = FALSE

# function
lr <- list() # begin return list
if(!check_donordf(df)){
  stop("Data.frame is not a valid simulated donor signals data.frame.")}
if(is(donor.unadj, "NULL")){donor.unadj <- df[,"donor.mean"]}
ktotal <- length(unique(df$type))
Zunadj <- matrix(donor.unadj, ncol = ktotal)

#punadj <- predtype(Z = Zunadj, Y = Ypb, strict_method = "nnls",
#                   proportions = TRUE, verbose = verbose)

Z = Zunadj
Y = Ypb
strict_method = "nnls"
proportions = TRUE
verbose = verbose

if (strict_method == "nnls") {
  p <- nnls::nnls(Z, Y)$x
} else {
  stop("Error, method not supported. Choose one of either: ", 
       paste0(names(supported_strict_methods()), collapse = ","))
}
p <- as.numeric(p)

if (proportions) {
  if (verbose) {
    message("Computing proportions from outputs.")
  }
  p <- p/sum(p)
} else {
  if (verbose) {
    message("Returning unmodified point prediction outputs.")
  }
}


# initial variable defs
prop.typev <- rep("punadj", ktotal)
ppredv <- punadj; ptruev <- P
type.indexv <- seq(ktotal)
lr[["donor.unadj"]] <- donor.unadj
if(!is(donor.adj.method, "NULL")){
  donor.adjv <- donoradj(df = df, donor.unadj = donor.unadj,
                         donor.adj.method = donor.adj.method, ...)
  lr[["donor.adj"]] <- donor.adjv
  Zadj <- matrix(donor.adjv, ncol = ktotal)
  padj <- predtype(Z = Zadj, Y = Ypb, strict_method = "nnls",
                   proportions = TRUE, verbose = verbose)
  # append to variable defs
  ptruev <- c(ptruev, P); ppredv <- c(ppredv, padj)
  prop.typev <- c(prop.typev, rep("padj", ktotal))
  type.indexv <- rep(type.indexv, 2)
  # parse plot arg
  if(plot.biasadj){
    dfp <- data.frame(unadj = lr[["donor.unadj"]], adj = lr[["donor.adj"]],
                      marker = df$marker, type = df$type)
    lr[["ggpt.biasadj"]] <- ggpt_donorbias(dfp, method.str = donor.adj.method)
  }
}
# append results
biasv <- ptruev - ppredv
lr[["dfres"]] <- data.frame(prop.type = prop.typev, prop.pred = ppredv,
                            prop.true = ptruev, bias = biasv,
                            type.index = type.indexv)

#---------------------------
# test run_donor_bias_expt()
#---------------------------

donordf = ddf # NULL
method = "nbinom"
lambda.pos = 20
lambda.neg = 2
lambda.sdoff.pos = 5
lambda.sdoff.neg = 2
gamma.pos = 20
gamma.neg = 2
P = c(0.25, 0.75)
donor.adj.method = "combat"
plot.biasadj = TRUE
plot.pca = TRUE
cname.donorsummary = "donor.mean"
gindexv = c(1, 2)
ndonor = 10
seed.num = 0
verbose = TRUE




set.seed(seed.num); lr <- list()
if(verbose){message("Making pseudobulk sample from types matrix...")}
df <- random_donordf(ndonor = 1, gindexv = gindexv,
                     lambda.sdoff.pos = 0, lambda.sdoff.neg = 0)
ktotal <- length(P)
Z <- matrix(df[,"donor1"], ncol = ktotal)
Ypb <- ypb_fromtypes(Z = Z, P = P)

if(is(donordf, "NULL")){
  if(verbose){message("Getting randomized donor marker data...")}  
  donordf <- random_donordf(ndonor = ndonor, method = method,
                            gindexv = gindexv, 
                            lambda.sdoff.pos = lambda.sdoff.pos,
                            lambda.sdoff.neg = lambda.sdoff.neg,
                            gamma.pos = gamma.pos, gamma.neg = gamma.neg,
                            seed.num = seed.num)
} else{
  if(verbose){message("Checking if provided donordf passes checks...")}
  if(!check_donordf(donordf)){
    stop("Error, provided donordf is invalid. ",
         "Check that it contains all required columns.")}
}

if(verbose){message("Getting type predictions...")}
type.indexv <- seq(ktotal)
donor.unadj <- donordf[,cname.donorsummary] # get donor summary datas
lbias <- biasexpt(df = donordf, Ypb = Ypb, P = P, donor.unadj = donor.unadj,
                  donor.adj.method = donor.adj.method,
                  plot.biasadj = plot.biasadj, verbose = verbose)

# get return object
lr[["dfres"]] <- lbias$dfres # get results df
lr[["donordf"]] <- donordf
lr[["Ypb"]] <- Ypb
lr[["adj.method"]] <- method
# make new plots
if(plot.pca){lr[["lpca"]] <- pcaplots_donor(donordf = donordf)}
if(plot.biasadj){lr[["ggpt.biasadj"]] <- lbias$ggpt.biasadj}