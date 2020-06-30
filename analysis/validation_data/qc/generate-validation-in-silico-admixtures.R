usePackage <- function(p) 
{
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, repos = "http://cran.us.r-project.org", dep = TRUE)
  require(p, character.only = TRUE)
}

usePackage("pacman")

home_dir <- "../../"
cur_dir <- getwd()
cache_dir <- paste0(home_dir, "/admixtures/new-admixtures/")

source("../../define-biological-models.R")
source("../../utils.R")

setwd(cache_dir)
source("utils.R")
source("load-data.R")
setwd(cur_dir)

## Spike-in dilution experiment. 

suppressPackageStartupMessages(p_load(synapser))
suppressPackageStartupMessages(p_load(plyr))
suppressPackageStartupMessages(p_load(dplyr))
suppressPackageStartupMessages(p_load(reshape2))
## install.packages("dirichlet", repos="http://R-Forge.R-project.org")
suppressPackageStartupMessages(p_load(dirichlet))
suppressPackageStartupMessages(p_load(hitandrun))

suppressPackageStartupMessages(p_load("foreach"))
suppressPackageStartupMessages(p_load("parallel"))

num.cores <- detectCores()
if(!is.na(num.cores) && (num.cores > 1)) {
  suppressPackageStartupMessages(p_load("doMC"))
  cat(paste("Registering ", num.cores-1, " cores.\n", sep=""))
  registerDoMC(cores=(num.cores-1))
}

synLogin()

## Load the TPM validation data
##synId <- "syn21574299"
##obj <- synGet(synId, downloadFile = TRUE)
##cpm.expr <- read.table(obj$path, sep = "\t", header = TRUE)

## this is symbol_tpm.csv
synId <- "syn21576632"
obj <- synGet(synId, downloadFile = TRUE)
cpm.expr <- read.table(obj$path, sep = ",", header = TRUE)

samples <- colnames(cpm.expr)
## Exclude the admixtures samples
flag <- grepl(samples, pattern="RM")
samples <- samples[!flag]
flag <- grepl(samples, pattern="BM")
samples <- samples[!flag]
samples <- samples[!(samples == "Gene")]
names(samples) <- samples

## Map populations to samples
populations <- samples
populations <- gsub(populations, pattern="_1", replacement="")
populations <- gsub(populations, pattern="_2", replacement="")
pop.df <- data.frame(population = unname(populations), sample = names(populations))

## Translate the population names to those we use in the challenge
fine.grained.map <- list(
    "Breast" = "Breast",
    "CRC" = "CRC",
    "Dendritic_cells" = "myeloid.dendritic.cells",
    "Endothelial_cells" = "endothelial.cells",
    "Fibroblasts" = "fibroblasts",
    "Macrophages" = "macrophages",
    "Memory_CD4_T_cells" = "memory.CD4.T.cells",
    "Memory_CD8_T_cells" = "memory.CD8.T.cells",
    "Monocytes" = "monocytes",
    "NK_cells" = "NK.cells",
    "Naive_B_cells" = "naive.B.cells",
    "Naive_CD4_T_cells" = "naive.CD4.T.cells",
    "Naive_CD8_T_cells" = "naive.CD8.T.cells",
    "Neutrophils" = "neutrophils",
    "Tregs" = "regulatory.T.cells")

coarse.grained.map <- list(
    "Breast" = "Breast",
    "CRC" = "CRC",
    "Dendritic_cells" = "monocytic.lineage",
    "Endothelial_cells" = "endothelial.cells",
    "Fibroblasts" = "fibroblasts",
    "Macrophages" = "monocytic.lineage",
    "Memory_CD4_T_cells" = "CD4.T.cells",
    "Memory_CD8_T_cells" = "CD8.T.cells",
    "Monocytes" = "monocytic.lineage",
    "NK_cells" = "NK.cells",
    "Naive_B_cells" = "B.cells",
    "Naive_CD4_T_cells" = "CD4.T.cells",
    "Naive_CD8_T_cells" = "CD8.T.cells",
    "Neutrophils" = "neutrophils",
    "Tregs" = "CD4.T.cells")
    

fine.grained.map.df <- data.frame(population = names(fine.grained.map), challenge.population = unname(unlist(fine.grained.map)))
fine.grained.pop.df <- merge(pop.df, fine.grained.map.df)

coarse.grained.map.df <- data.frame(population = names(coarse.grained.map), challenge.population = unname(unlist(coarse.grained.map)))
coarse.grained.pop.df <- merge(pop.df, coarse.grained.map.df)

## Separate the samples into two batches
fine.grained.pop.df <- fine.grained.pop.df[order(fine.grained.pop.df$sample),]
fine.grained.pop1.df <- ddply(fine.grained.pop.df, .variables = c("population"), .fun = function(df) df[1,,drop=F])
fine.grained.pop2.df <- ddply(fine.grained.pop.df, .variables = c("population"), .fun = function(df) df[min(2,nrow(df)),,drop=F])

coarse.grained.pop.df <- coarse.grained.pop.df[order(coarse.grained.pop.df$sample),]
coarse.grained.pop1.df <- ddply(coarse.grained.pop.df, .variables = c("population"), .fun = function(df) df[1,,drop=F])
coarse.grained.pop2.df <- ddply(coarse.grained.pop.df, .variables = c("population"), .fun = function(df) df[min(2,nrow(df)),,drop=F])

immune.fine.grained.populations <- unique(as.character(fine.grained.pop.df$challenge.population))
immune.fine.grained.populations <- immune.fine.grained.populations[immune.fine.grained.populations != "CRC"]
immune.fine.grained.populations <- immune.fine.grained.populations[immune.fine.grained.populations != "Breast"]
num.immune.fine.grained.populations <- length(immune.fine.grained.populations)

immune.coarse.grained.populations <- unique(as.character(coarse.grained.pop.df$challenge.population))
immune.coarse.grained.populations <- immune.coarse.grained.populations[immune.coarse.grained.populations != "CRC"]
immune.coarse.grained.populations <- immune.coarse.grained.populations[immune.coarse.grained.populations != "Breast"]
num.immune.coarse.grained.populations <- length(immune.coarse.grained.populations)

create.rand.admixture <- function(populations, n.admixs) {
    num.pops <- length(populations)
    df <- generate.random.uniform.admixtures(populations, n = n.admixs, min.prop = 0, max.prop = 1, step = 0.001)
    admixs <- t(df)
}

generate.random.admixtures <- function(populations, population.maps, n.admixs) {

    rand.admixtures <-
        ldply(population.maps,
              .fun = function(pop.df) {
                  ret <- create.rand.admixture(populations, n.admixs)
                  long <- melt(as.matrix(ret))
                  colnames(long) <- c("admixture", "challenge.population", "prop")
                  long$challenge.population <- as.character(long$challenge.population)
                  col <- "sample"
                  long <- merge(long, pop.df[, c("challenge.population", col)], all.x = TRUE)
                  long <-
                      ddply(long, .variables = c("challenge.population", "admixture"),
                            .fun = function(df) {
                                if(nrow(df) == 1) { return(df) }
                                ## Distribute the proportion across the samples mapped to this
                                ## challenge population
                                ## (e.g., in the coarse-grained challenge memory CD4, naive CD4,
                                ## and Tregs map to CD4)
                                props <- as.numeric(rdirichlet(n = 1, alpha = rep(1/nrow(df), nrow(df))))
                                df$prop <- df$prop * props
                                df
                            })
                  long <- long[, c("admixture", col, "prop")]
                  long[,col] <- as.character(long[,col])
                  long[,"prop"] <- as.numeric(long[,"prop"])
                  long
              })
    
    rand.admixtures$admixture <- paste0(rand.admixtures$.id, rand.admixtures$admixture)
    rand.admixtures <- rand.admixtures[, !(colnames(rand.admixtures) == ".id")]

    ## Confirm all admixtures sum to one
    sums <-
        ddply(rand.admixtures, .variables = c("admixture"),
              .fun = function(df) sum(df$prop))
    eps <- 10^-5
    flag <- abs(sums$V1 - 1) < eps 
    if(!(all(flag))) {
        print(head(sums[!flag,]))
        stop("Some admixtures don't sum to one\n")
    }

    rand.admixtures
}

generate.biological.admixtures_ <- function(model, populations, n = 15) {
    model <- model[, populations]
    ndim <- ncol(model)
    n.tot <- 1000
    ##N <- 1000*ndim^3
    ## hitandrun documentation says that, since har uses MCMC which is time-correlated,
    ## need to run for O(n^3) steps.
    prefactor <- 100
    if(prefactor < ndim) { stop("Expected prefator >= ndim to ensure O(n^3) steps\n") }
    N <- prefactor*(ndim)^3
    N.thin <- N
    n.rand <- 5
    N <- 10*ndim^3
    N.thin <- 10*ndim^3
    ## despite what documentation says, output has length N / N.thin (not N * N.thin)
    N <- 1000*(ndim)^3
    N.thin <- (ndim)^3

    eq.constr <- list(
        constr = matrix(rep(1,ndim),nrow=1,ncol=(ndim)),
        dir = '=',
        rhs = 1)

    ineq.constr <- list(
        constr = rbind(-diag(ndim),diag(ndim)),
        dir = rep('<=',(ndim)*2),
        rhs = c(-as.numeric(model["min",]),as.numeric(model["max",])))
    
    basis <- solution.basis(eq.constr)
    transform <- createTransform(basis)
    constr <- transformConstraints(transform,ineq.constr)
    
    x0 <- createSeedPoint(constr,homogeneous=TRUE,randomize=TRUE)
    df <- har(x0,constr,N,thin=N.thin,transform=transform, homogeneous=TRUE)$samples
    eps <- 10^-4
    flag <- abs( rowSums(df) - 1 ) < eps
    df <- df[flag, ,drop=F]
    df <- df[!duplicated(df),,drop=F]
    colnames(df) <- colnames(model)
    admixtures <- df

    t(find.extremal.points(t(admixtures), n))
}

generate.biological.admixtures <- function(model, populations, population.maps, n = 15) {

    admixtures <-
        ldply(population.maps,
              .fun = function(pop.df) {
                  ret <- generate.biological.admixtures_(model, populations, n = n)
                  long <- melt(as.matrix(ret))
                  colnames(long) <- c("admixture", "challenge.population", "prop")
                  long$challenge.population <- as.character(long$challenge.population)
                  pop.df$challenge.population <- as.character(pop.df$challenge.population)
                  col <- "sample"
                  long <- merge(long, pop.df[, c("challenge.population", col)], all.x = TRUE)
                  long <-
                      ddply(long, .variables = c("challenge.population", "admixture"),
                            .fun = function(df) {
                                if(nrow(df) == 1) { return(df) }
                                ## Distribute the proportion across the samples mapped to this
                                ## challenge population
                                ## (e.g., in the coarse-grained challenge memory CD4, naive CD4,
                                ## and Tregs map to CD4)
                                props <- as.numeric(rdirichlet(n = 1, alpha = rep(1/nrow(df), nrow(df))))
                                old <- df$prop[1]
                                if(!all(df$prop == df$prop[1])) {
                                    stop(paste0("Not all props are the same for pop ", df$challenge.population[1], "\n"))
                                }
                                eps <- 10^-4
                                df$prop <- df$prop * props
                                if(abs(sum(df$prop) - old) > eps) {
                                    print(props)
                                    stop(paste0(sum(df$prop), " != ", old, "\n"))
                                }
                                df
                            })
                  long <- long[, c("admixture", col, "prop")]
                  long[,col] <- as.character(long[,col])
                  long[,"prop"] <- as.numeric(long[,"prop"])
                  long
              })
    colnames(admixtures)[1] <- "pop.id"
    
    admixtures$admixture <- paste0(admixtures$pop.id, admixtures$admixture)
    admixtures <- admixtures[, !(colnames(admixtures) %in% c("pop.id"))]

    ## Confirm all admixtures sum to one
    sums <-
        ddply(admixtures, .variables = c("admixture"),
              .fun = function(df) sum(df$prop))
    eps <- 10^-5
    flag <- abs(sums$V1 - 1) < eps 
    if(!(all(flag))) {
        print(head(sums[!flag,]))
        stop("Some admixtures don't sum to one\n")
    }

    admixtures
}

flat.model <- define.biological.model(unname(unlist(populations)), cache.dir = cache_dir)
tmp <- subset(fine.grained.map.df, population %in% colnames(flat.model))
old.cols <- c(as.character(tmp$population), "tumor.fraction")
new.cols <- c(as.character(tmp$challenge.population), "tumor.fraction")
flat.model <- flat.model[, old.cols]
colnames(flat.model) <- new.cols

flat.model.coarse <- flat.model
flat.model.coarse <- cbind(flat.model.coarse, "monocytic.lineage" =
    flat.model.coarse[, "myeloid.dendritic.cells"] +
    flat.model.coarse[, "macrophages"] +
    flat.model.coarse[, "monocytes"])
flat.model.coarse <- cbind(flat.model.coarse, "CD4.T.cells" =
    flat.model.coarse[, "memory.CD4.T.cells"] +
    flat.model.coarse[, "naive.CD4.T.cells"] +
    flat.model.coarse[, "regulatory.T.cells"])
flat.model.coarse <- cbind(flat.model.coarse, "CD8.T.cells" =
    flat.model.coarse[, "memory.CD8.T.cells"] +
    flat.model.coarse[, "naive.CD8.T.cells"])
flat.model.coarse <- cbind(flat.model.coarse,
                           "B.cells" = flat.model.coarse[, "naive.B.cells"])
## all of the population minimums were set to 0.01, so when we add those we get something higher.
## reset to 0.01
flat.model.coarse["min", "monocytic.lineage"] <- 0.01
flat.model.coarse["min", "CD4.T.cells"] <- 0.01
flat.model.coarse["min", "CD8.T.cells"] <- 0.01

tmp <- subset(coarse.grained.map.df, challenge.population %in% colnames(flat.model.coarse))
flat.model.coarse <- flat.model.coarse[, unique(as.character(tmp$challenge.population), "tumor.fraction")]

fine.grained.pop.maps <- list("A" = fine.grained.pop1.df, "B" = fine.grained.pop2.df)
coarse.grained.pop.maps <- list("A" = coarse.grained.pop1.df, "B" = coarse.grained.pop2.df)

set.seed(100)

n.admixs <- 15
cancers <- c("Breast", "CRC")
names(cancers) <- cancers

coarse.grained.rand.admixtures <- 
    llply(cancers,
          .fun = function(cancer) {
              ret <- generate.random.admixtures(c(immune.coarse.grained.populations, cancer), coarse.grained.pop.maps, n.admixs = n.admixs)
              ret$admixture <- paste0("ISCR", ret$admixture)
              ret
          })

fine.grained.rand.admixtures <- 
    llply(cancers,
          .fun = function(cancer) {
              ret <- generate.random.admixtures(c(immune.fine.grained.populations, cancer), fine.grained.pop.maps, n.admixs = n.admixs)
              ret$admixture <- paste0("ISFR", ret$admixture)
              ret
              
          })

## Generate fine-grained biological admixtures w/ cancer
## For each cancer, generate 10 replicates for each of two sets of samples

num.w.cancer.datasets <- 5
nums <- 1:num.w.cancer.datasets
names(nums) <- nums

cat("Generating fine-grained biological admixtures\n")
fine.grained.bio.admixtures <-
    llply(cancers,
          function(cancer) {
              tmp <-
                  ldply(nums,
                        .parallel = TRUE,
                        .fun = function(i) {
                            mdl <- flat.model
                            colnames(mdl)[colnames(mdl) == "tumor.fraction"] <- cancer
                            populations <- colnames(mdl)
                            ret <- generate.biological.admixtures(mdl, populations,
                                                                  fine.grained.pop.maps, n = 4)
                        })
              
              if(!(".id" %in% colnames(tmp))) {
                  stop(paste0(".id not in tmp\n"))
              }
              tmp$admixture <- paste0("ISFB", tmp$.id, tmp$admixture)
              tmp <- tmp[, !(colnames(tmp) %in% ".id")]
              tmp
          })

cat("Generating coarse-grained biological admixtures\n")
coarse.grained.bio.admixtures <-
    llply(cancers,
          function(cancer) {
              tmp <-
                  ldply(nums,
                        .parallel = TRUE,
                        .fun = function(i) {
                            mdl <- flat.model.coarse
                            colnames(mdl)[colnames(mdl) == "tumor.fraction"] <- cancer
                            populations <- colnames(mdl)
                            ret <- generate.biological.admixtures(mdl, populations,
                                                                  coarse.grained.pop.maps, n = 4)
                        })
              
              if(!(".id" %in% colnames(tmp))) {
                  stop(paste0(".id not in tmp\n"))
              }
              tmp$admixture <- paste0("ISCB", tmp$.id, tmp$admixture)
              tmp <- tmp[, !(colnames(tmp) %in% ".id")]
              tmp
          })

nxt.letter <- 1
datasets <- list()

my.letters <- c(paste0("A", LETTERS))

l <- fine.grained.rand.admixtures
for(nm in names(l)) {
    ds <- l[[nm]]
    tt <- nm
    if(tt == "Breast") { tt <- "BRCA" }
    datasets[[my.letters[nxt.letter]]] <- list("data" = ds, "mixture.type" = "Random", "subchallenge" = "fine", "tumor.type" = tt)
    nxt.letter <- nxt.letter + 1
}

l <- coarse.grained.rand.admixtures
for(nm in names(l)) {
    ds <- l[[nm]]
    tt <- nm
    if(tt == "Breast") { tt <- "BRCA" }
    datasets[[my.letters[nxt.letter]]] <- list("data" = ds, "mixture.type" = "Random", "subchallenge" = "coarse", "tumor.type" = tt)
    nxt.letter <- nxt.letter + 1
}

l <- fine.grained.bio.admixtures
for(nm in names(l)) {
    ds <- l[[nm]]
    tt <- nm
    if(tt == "Breast") { tt <- "BRCA" }
    datasets[[my.letters[nxt.letter]]] <- list("data" = ds, "mixture.type" = "Biological", "subchallenge" = "fine", "tumor.type" = tt)
    nxt.letter <- nxt.letter + 1
}

l <- coarse.grained.bio.admixtures
for(nm in names(l)) {
    ds <- l[[nm]]
    tt <- nm
    if(tt == "Breast") { tt <- "BRCA" }
    datasets[[my.letters[nxt.letter]]] <- list("data" = ds, "mixture.type" = "Biological", "subchallenge" = "coarse", "tumor.type" = tt)
    nxt.letter <- nxt.letter + 1
}

metadata <- ldply(datasets, .fun = function(entry) as.data.frame(entry[names(entry) != "data"], stringsAsFactors = FALSE))
colnames(metadata)[1] <- "dataset.name"
coarse.datasets <- subset(metadata, subchallenge == "coarse")$dataset.name
fine.datasets <- subset(metadata, subchallenge == "fine")$dataset.name

metadata.file <- "in-silico-validation-metadata.csv"
write.table(file = metadata.file, metadata, row.names = FALSE, col.names = TRUE, sep = ",", quote = FALSE)

## Create the gold standards
## gold standard files should be csvs with columns: dataset.name,sample.id,cell.type,measured

## NB: these need to be sum'ed across cell populations that map to the same challenge subtype
## NB: ideally, include both in the ground truth

coarse.gs <- ldply(datasets, .fun = function(df) df$data)
colnames(coarse.gs) <- c("dataset.name", "sample.id", "sample", "measured")
coarse.gs <- subset(coarse.gs, dataset.name %in% coarse.datasets)
old.sz <- nrow(coarse.gs)
coarse.gs <- merge(coarse.gs, coarse.grained.pop.df[, c("sample", "challenge.population")], by = "sample")
new.sz <- nrow(coarse.gs)
if(old.sz != new.sz) { stop(paste0("Size changed from ", old.sz, " to ", new.sz, "\n")) }
coarse.gs <- coarse.gs[, c("dataset.name", "sample.id", "challenge.population", "sample", "measured")]
colnames(coarse.gs) <- c("dataset.name", "sample.id", "sample", "vendor.sample", "measured")

fine.gs <- ldply(datasets, .fun = function(df) df$data)
colnames(fine.gs) <- c("dataset.name", "sample.id", "sample", "measured")
fine.gs <- subset(fine.gs, dataset.name %in% fine.datasets)
old.sz <- nrow(fine.gs)
fine.gs <- merge(fine.gs, fine.grained.pop.df[, c("sample", "challenge.population")], by = "sample")
new.sz <- nrow(fine.gs)
if(old.sz != new.sz) { stop(paste0("Size changed from ", old.sz, " to ", new.sz, "\n")) }
fine.gs <- fine.gs[, c("dataset.name", "sample.id", "challenge.population", "sample", "measured")]
colnames(fine.gs) <- c("dataset.name", "sample.id", "sample", "vendor.sample", "measured")


## this is symbol_tpm.csv
synId <- "syn21576632"
obj <- synGet(synId, downloadFile = TRUE)
cpm.expr <- read.table(obj$path, sep = ",", header = TRUE)

rownames(cpm.expr) <- as.character(cpm.expr$Gene)
cpm.expr <- cpm.expr[, !(colnames(cpm.expr) %in% c("Gene"))]

## this is ensg_tpm.csv
synId <- "syn21576631"
obj <- synGet(synId, downloadFile = TRUE)
ensg.cpm.expr <- read.table(obj$path, sep = ",", header = TRUE)

rownames(ensg.cpm.expr) <- as.character(ensg.cpm.expr$Gene)
ensg.cpm.expr <- ensg.cpm.expr[, !(colnames(ensg.cpm.expr) %in% c("Gene"))]

## Load the counts
## this is symbol_counts.csv
synId <- "syn21576630"
obj <- synGet(synId, downloadFile = TRUE)
cnts.expr <- read.table(obj$path, sep = ",", header = TRUE)

rownames(cnts.expr) <- as.character(cnts.expr$Gene)
cnts.expr <- cnts.expr[, !(colnames(cnts.expr) %in% c("Gene"))]

## this is ensg_counts.csv
synId <- "syn21576629"
obj <- synGet(synId, downloadFile = TRUE)
ensg.cnts.expr <- read.table(obj$path, sep = ",", header = TRUE)

rownames(ensg.cnts.expr) <- as.character(ensg.cnts.expr$Gene)
ensg.cnts.expr <- ensg.cnts.expr[, !(colnames(ensg.cnts.expr) %in% c("Gene"))]

coarse.gs <- ldply(datasets, .fun = function(df) df$data)
colnames(coarse.gs) <- c("dataset.name", "sample.id", "sample", "measured")
coarse.gs <- subset(coarse.gs, dataset.name %in% coarse.datasets)

fine.gs <- ldply(datasets, .fun = function(df) df$data)
colnames(fine.gs) <- c("dataset.name", "sample.id", "sample", "measured")
fine.gs <- subset(fine.gs, dataset.name %in% fine.datasets)

all.gs <- rbind(coarse.gs, fine.gs)

## Create the expression files
## First column is Gene

if(FALSE) {
admixtures <-
    dlply(all.gs,
          .parallel = TRUE,
          .variables = c("dataset.name"),
          .fun = function(df) {
              mxs <-
                  dlply(df,
                        .variables = c("sample.id"),
                        .fun = function(df) {
                            cols <- as.character(df$sample)
                            fracs <- as.numeric(df$measured)
                            if(!all(cols %in% colnames(cpm.expr))) {
                                print(cols)
                                stop("Missing some columns\n")
                            }
                            mat <- cpm.expr[, cols]
                            ret <- as.matrix(mat) %*% fracs
                            colnames(ret) <- df$sample.id[1]
                            ret
                        })
              mxs <- do.call(cbind, mxs)
              mxs <- cbind(Gene = rownames(cpm.expr), mxs)
              mxs
          })
} ## FALSE

admixtures <-
    dlply(all.gs,
          .parallel = TRUE,
          .variables = c("dataset.name"),
          .fun = function(df) {
              mat <- as.matrix(cpm.expr)
              mxs <- create.in.silico.admixtures(mat, df,
                                                 sample.col = "sample.id", cell.type.col = "sample",
                                                 fraction.col = "measured")
              mxs <- cbind(Gene = rownames(mat), mxs)
              mxs
          })

ensg.admixtures <-
    dlply(all.gs,
          .parallel = TRUE,
          .variables = c("dataset.name"),
          .fun = function(df) {
              mat <- as.matrix(ensg.cpm.expr)
              mxs <- create.in.silico.admixtures(mat, df,
                                                 sample.col = "sample.id", cell.type.col = "sample",
                                                 fraction.col = "measured")
              mxs <- cbind(Gene = rownames(mat), mxs)
              mxs
          })

## Create the same admixtures, but for count matrices
purified.samples <- sort(as.character(unique(all.gs$sample)))

purified.cnts.expr <- cnts.expr[, purified.samples]
rownames(purified.cnts.expr) <- cnts.expr$Gene
tot.purified.cnts <- colSums(purified.cnts.expr)
median.tot.purified.cnts <- median(tot.purified.cnts)

norm.purified.cnts.expr <- sweep(purified.cnts.expr, 2, colSums(purified.cnts.expr),`/`) * median.tot.purified.cnts
rownames(norm.purified.cnts.expr) <- rownames(cnts.expr)

ensg.purified.cnts.expr <- ensg.cnts.expr[, purified.samples]
rownames(ensg.purified.cnts.expr) <- ensg.cnts.expr$Gene
ensg.tot.purified.cnts <- colSums(ensg.purified.cnts.expr)
ensg.median.tot.purified.cnts <- median(ensg.tot.purified.cnts)

ensg.norm.purified.cnts.expr <- sweep(ensg.purified.cnts.expr, 2, colSums(ensg.purified.cnts.expr),`/`) * ensg.median.tot.purified.cnts
rownames(ensg.norm.purified.cnts.expr) <- rownames(ensg.cnts.expr)

cnts.admixtures <-
    dlply(all.gs,
          .parallel = TRUE,
          .variables = c("dataset.name"),
          .fun = function(df) {
              mat <- as.matrix(norm.purified.cnts.expr)
              mxs <- create.in.silico.admixtures(mat, df,
                                                 sample.col = "sample.id", cell.type.col = "sample",
                                                 fraction.col = "measured")
              mxs <- cbind(Gene = rownames(mat), mxs)
              mxs
          })

ensg.cnts.admixtures <-
    dlply(all.gs,
          .parallel = TRUE,
          .variables = c("dataset.name"),
          .fun = function(df) {
              mat <- as.matrix(ensg.norm.purified.cnts.expr)
              mxs <- create.in.silico.admixtures(mat, df,
                                                 sample.col = "sample.id", cell.type.col = "sample",
                                                 fraction.col = "measured")
              mxs <- cbind(Gene = rownames(mat), mxs)
              mxs
          })

## Convert cnts to cpms
suppressPackageStartupMessages(p_load(edgeR))
cnts.to.cpm.admixtures <-
    llply(cnts.admixtures,
          .fun = function(df) {
              rownames(df) <- df[, "Gene"]
              df <- df[, !(colnames(df) == "Gene")]
              df <- as.matrix(df)
              class(df) <- "numeric"
              df <- round(df)
              cpm.df <- cpm(df, log = FALSE)
              cbind(Gene = rownames(df), cpm.df)
          })

ensg.cnts.to.cpm.admixtures <-
    llply(ensg.cnts.admixtures,
          .fun = function(df) {
              rownames(df) <- df[, "Gene"]
              df <- df[, !(colnames(df) == "Gene")]
              df <- as.matrix(df)
              class(df) <- "numeric"
              df <- round(df)
              cpm.df <- cpm(df, log = FALSE)
              cbind(Gene = rownames(df), cpm.df)
          })


print(warnings())

## Create the input.csv file

cat("Writing input.csv")

input.tbl <-
    ldply(1:nrow(metadata),
          .fun = function(i) {
              ds <- as.character(metadata[i, "dataset.name"])
              ct <- as.character(metadata[i, "tumor.type"])
              params <- list(
                  "dataset.name" = ds,
                  "cancer.type" = ct,
                  "platform" = "Illumina",
                  "scale" = "Linear",
                  "normalization" = "TPM",
                  "native.probe.type" = "ENSG",
                  "symbol.compression.function" = "colMeans",
                  "ensg.compression.function" = "colMeans",
                  "symbol.to.native.mapping.file" = "native_to_hugo.tsv",
                  "ensg.to.native.mapping.file" = "native_to_ensg.tsv",
                  "hugo.expr.file" = paste0(ds, "_symbol_val_tpm.csv"),
                  "hugo.expr.est.counts.file" = paste0(ds, "_symbol_val_cnts.csv"),
                  "ensg.expr.file" = paste0(ds, "_ensg_val_tpm.csv"),
                  "ensg.expr.est.counts.file" = paste0(ds, "_ensg_val_cnts.csv"),
                  "fastq.samples" = NA,                  
                  "fastq1.files" = NA,
                  "fastq2.files" = NA,
                  "native.expr.file" = paste0(ds, "_ensg_val_tpm.csv"),
                  "native.expr.est.counts.file" = paste0(ds, "_ensg_val_cnts.csv"),
                  "symbol.compression.est.counts.function" = "colMeans",
                  "ensg.compression.est.counts.function" = "colMeans"
              )
              as.data.frame(params)
          })

## Create a dummy table where the symbol tpm files are generated by log'ing the counts.
## Since all of the methods use cpm/tpm -- this will allow me to confirm that the counts
## data are OK.
input.cpm.cnts.tbl <-
    ldply(1:nrow(metadata),
          .fun = function(i) {
              ds <- as.character(metadata[i, "dataset.name"])
              ct <- as.character(metadata[i, "tumor.type"])
              params <- list(
                  "dataset.name" = ds,
                  "cancer.type" = ct,
                  "platform" = "Illumina",
                  "scale" = "Linear",
                  "normalization" = "TPM",
                  "native.probe.type" = "ENSG",
                  "symbol.compression.function" = "colMeans",
                  "ensg.compression.function" = "colMeans",
                  "symbol.to.native.mapping.file" = "native_to_hugo.tsv",
                  "ensg.to.native.mapping.file" = "native_to_ensg.tsv",
                  "hugo.expr.file" = paste0(ds, "_symbol_val_cnts_to_cpm.csv"),
                  "hugo.expr.est.counts.file" = paste0(ds, "_symbol_val_cnts.csv"),
                  "ensg.expr.file" = paste0(ds, "_ensg_val_cnts_to_cpm.csv"),
                  "ensg.expr.est.counts.file" = paste0(ds, "_ensg_val_cnts.csv"),
                  "fastq.samples" = NA,                  
                  "fastq1.files" = NA,
                  "fastq2.files" = NA,
                  "native.expr.file" = paste0(ds, "_ensg_val_cnts_to_cpm.csv"),
                  "native.expr.est.counts.file" = paste0(ds, "_ensg_val_cnts.csv"),
                  "symbol.compression.est.counts.function" = "colMeans",
                  "ensg.compression.est.counts.function" = "colMeans"
              )
              as.data.frame(params)
          })

## "in silico validation phase"
parent.id <- "syn22013292"

file <- "insilico-validation-input.csv"
write.table(file = file, input.tbl, col.names = TRUE, row.names = FALSE, sep = ",", quote = FALSE)
cat(paste0("Storing ", file, " to synapse\n"))
f <- File(file, parentId = parent.id, synapseStore = TRUE)
synStore(f)

file <- "insilico-validation-input-cpm-cnts.csv"
write.table(file = file, input.cpm.cnts.tbl, col.names = TRUE, row.names = FALSE, sep = ",", quote = FALSE)
cat(paste0("Storing ", file, " to synapse\n"))
f <- File(file, parentId = parent.id, synapseStore = TRUE)
synStore(f)

nms <- names(admixtures)
names(nms) <- nms
for(nm in nms) {
    cat(paste0("Writing ", nm, " : ", nrow(admixtures[[nm]]), " x ", ncol(admixtures[[nm]]), "\n"))
    file <- paste0(nm, "_symbol_val_tpm.csv")
    write.table(file = file, admixtures[[nm]], sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)
    cat(paste0("Storing ", nm, " to synapse\n"))
    f <- File(file, parentId = parent.id, synapseStore = TRUE)
    synStore(f)
}

nms <- names(ensg.admixtures)
names(nms) <- nms
for(nm in nms) {
    cat(paste0("Writing ", nm, " : ", nrow(ensg.admixtures[[nm]]), " x ", ncol(ensg.admixtures[[nm]]), "\n"))
    file <- paste0(nm, "_ensg_val_tpm.csv")
    write.table(file = file, ensg.admixtures[[nm]], sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)
    cat(paste0("Storing ", nm, " to synapse\n"))
    f <- File(file, parentId = parent.id, synapseStore = TRUE)
    synStore(f)
}

nms <- names(cnts.admixtures)
names(nms) <- nms
for(nm in nms) {
    cat(paste0("Writing ", nm, " : ", nrow(cnts.admixtures[[nm]]), " x ", ncol(cnts.admixtures[[nm]]), "\n"))
    file <- paste0(nm, "_symbol_val_cnts.csv")
    write.table(file = file, cnts.admixtures[[nm]], sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)
    cat(paste0("Storing ", nm, " to synapse\n"))
    f <- File(file, parentId = parent.id, synapseStore = TRUE)
    synStore(f)
}

nms <- names(ensg.cnts.admixtures)
names(nms) <- nms
for(nm in nms) {
    cat(paste0("Writing ", nm, " : ", nrow(ensg.cnts.admixtures[[nm]]), " x ", ncol(ensg.cnts.admixtures[[nm]]), "\n"))
    file <- paste0(nm, "_ensg_val_cnts.csv")
    write.table(file = file, ensg.cnts.admixtures[[nm]], sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)
    cat(paste0("Storing ", nm, " to synapse\n"))
    f <- File(file, parentId = parent.id, synapseStore = TRUE)
    synStore(f)
}

nms <- names(ensg.cnts.to.cpm.admixtures)
names(nms) <- nms
for(nm in nms) {
    cat(paste0("Writing ", nm, " : ", nrow(ensg.cnts.to.cpm.admixtures[[nm]]), " x ", ncol(ensg.cnts.to.cpm.admixtures[[nm]]), "\n"))
    file <- paste0(nm, "_ensg_val_cnts_to_cpm.csv")
    write.table(file = file, ensg.cnts.to.cpm.admixtures[[nm]], sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)
    cat(paste0("Storing ", nm, " to synapse\n"))
    f <- File(file, parentId = parent.id, synapseStore = TRUE)
    synStore(f)
}

nms <- names(cnts.to.cpm.admixtures)
names(nms) <- nms
for(nm in nms) {
    cat(paste0("Writing ", nm, " : ", nrow(cnts.to.cpm.admixtures[[nm]]), " x ", ncol(cnts.to.cpm.admixtures[[nm]]), "\n"))
    file <- paste0(nm, "_symbol_val_cnts_to_cpm.csv")
    write.table(file = file, cnts.to.cpm.admixtures[[nm]], sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)
    cat(paste0("Storing ", nm, " to synapse\n"))
    f <- File(file, parentId = parent.id, synapseStore = TRUE)
    synStore(f)
}



cat("Done\n")
save.image(".Rdata")
cat("Done saving\n")

file <- metadata.file
f <- File(file, parentId = parent.id, synpaseStore = TRUE)
synStore(f)

## Write out the admixtures and store to synapse
l <- list("in-silico-validation-fine.csv" = fine.gs, "in-silico-validation-coarse.csv" = coarse.gs)
for(nm in names(l)) {
    file <- nm
    write.table(file = file, l[[nm]], sep=",", col.names = TRUE, row.names = FALSE, quote = FALSE)
    f <- File(file, parentId = parent.id, synapseStore = TRUE)
    synStore(f)
}

## Write out the ground truth
tbl <- get.purified.sample.translation.table()
fine.gt <- merge(fine.gs, tbl[, c("sample", "fine.grained.cell.type")], by = "sample")
coarse.gt <- merge(coarse.gs, tbl[, c("sample", "coarse.grained.cell.type")], by = "sample")

coarse.gt <-
    ddply(coarse.gt,
          .variables = c("dataset.name", "sample.id", "coarse.grained.cell.type"),
          .fun = function(df) { data.frame(measured = sum(df$measured)) })

coarse.gt <- coarse.gt[, c("dataset.name", "sample.id", "coarse.grained.cell.type", "measured")]
colnames(coarse.gt) <- c("dataset.name", "sample.id", "cell.type", "measured")
fine.gt <- fine.gt[, c("dataset.name", "sample.id", "fine.grained.cell.type", "measured")]
colnames(fine.gt) <- c("dataset.name", "sample.id", "cell.type", "measured")

l <- list("in-silico-gt-fine.csv" = fine.gt, "in-silico-gt-coarse.csv" = coarse.gt)
for(nm in names(l)) {
    file <- nm
    write.table(file = file, l[[nm]], sep=",", col.names = TRUE, row.names = FALSE, quote = FALSE)    
    f <- File(file, parentId = parent.id, synapseStore = TRUE)
    synStore(f)
}


cat("Exiting\n")
q(status = 0)

## Thurs/Fri:
## 6. run baseline methods against tpm
## 7. run csx methods against tpm
## 8. look at results from baseline(tpm); separated by dataset
## 9. make dummy input file switching to normalized(cnts)
## 10. run baseline methods against normalized(cnts)
## 11. run csx methods against normalized(cnts)
## 12. look at results from baseline(normalized(cnts)); separated by dataset


