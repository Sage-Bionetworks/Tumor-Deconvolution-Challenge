usePackage <- function(p) 
{
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, repos = "http://cran.us.r-project.org", dep = TRUE)
  require(p, character.only = TRUE)
}

usePackage("pacman")

home_dir <- "../"
cur_dir <- getwd()

setwd(home_dir)
source("admixtures/utils.R")
source("admixtures/load-data.R")
source("admixtures/define-biological-models.R")
setwd(cur_dir)

## Spike-in dilution experiment. 

suppressPackageStartupMessages(p_load(synapser))
suppressPackageStartupMessages(p_load(plyr))
suppressPackageStartupMessages(p_load(dplyr))
suppressPackageStartupMessages(p_load(reshape2))
## install.packages("dirichlet", repos="http://R-Forge.R-project.org")
suppressPackageStartupMessages(p_load(dirichlet))


synLogin()

## Load the TPM validation data
##synId <- "syn21574299"
##obj <- synGet(synId, downloadFile = TRUE)
##cpm.expr <- read.table(obj$path, sep = "\t", header = TRUE)

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

set.seed(100)

immune.fine.grained.populations <- unique(as.character(fine.grained.pop.df$challenge.population))
immune.fine.grained.populations <- immune.fine.grained.populations[immune.fine.grained.populations != "CRC"]
immune.fine.grained.populations <- immune.fine.grained.populations[immune.fine.grained.populations != "Breast"]
num.immune.fine.grained.populations <- length(immune.fine.grained.populations)

immune.coarse.grained.populations <- unique(as.character(coarse.grained.pop.df$challenge.population))
immune.coarse.grained.populations <- immune.coarse.grained.populations[immune.coarse.grained.populations != "CRC"]
immune.coarse.grained.populations <- immune.coarse.grained.populations[immune.coarse.grained.populations != "Breast"]
num.immune.coarse.grained.populations <- length(immune.coarse.grained.populations)

create.rand.admixture.wo.cancer <- function(populations, spike.ins) {

    num.pops <- length(populations)
    admixs <-
        ldply(spike.ins,
              .fun = function(spike.in) {
                  dummy <- 1:(num.pops-1)
                  df <- generate.random.uniform.admixtures(dummy, n = 5, min.prop = 0, max.prop = 1 - spike.in, step = 0.001)
                  t(df)
              })
    admixs[,1] <- as.numeric(admixs[,1])

    ret <- ldply(1:num.pops,
                 .fun = function(i) {
                     tmp <- admixs
                     colnames(tmp) <- populations
                     tmp.name <- colnames(tmp)[1]
                     colnames(tmp)[1] <- colnames(tmp)[i]
                     colnames(tmp)[i] <- tmp.name
                     tmp[, populations]
                 })
    return(ret)
}

create.rand.admixture.w.cancer <- function(non.cancer.populations, cancer, spike.ins) {

    num.pops <- length(non.cancer.populations)
    admixs <-
        ldply(spike.ins,
              .fun = function(spike.in) {
                  dummy <- 1:num.pops
                  df <- generate.random.uniform.admixtures(dummy, n = 15, min.prop = 0, max.prop = 1 - spike.in, step = 0.001)
                  t(df)
              })
    admixs[,1] <- as.numeric(admixs[,1])

    colnames(admixs) <- c(cancer, non.cancer.populations)

    admixs
}

generate.random.admixtures.wo.cancer <- function(immune.populations, population.maps, spike.ins) {

    rand.admixtures <-
        ldply(population.maps,
              .fun = function(pop.df) {
                  ret <- create.rand.admixture.wo.cancer(immune.populations, spike.ins)
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
    if(!(all(abs(sums$V1 - 1) < eps))) { stop("Some admixtures don't sum to one\n") }

    rand.admixtures
}

generate.random.admixtures.w.cancer <- function(immune.populations, population.maps, spike.ins, cancers = c("Breast", "CRC")) {

    names(cancers) <- c("C", "D")
    rand.admixtures <-
        ldply(cancers,
              .fun = function(cancer) {
                  tmp <- ldply(population.maps,
                        .fun = function(pop.df) {
                            ret <- create.rand.admixture.w.cancer(immune.populations, cancer, spike.ins)
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
                                          df$prop <- df$prop * props
                                          df
                                      })
                            long <- long[, c("admixture", col, "prop")]
                            long[,col] <- as.character(long[,col])
                            long[,"prop"] <- as.numeric(long[,"prop"])
                            long
                        })
                  colnames(tmp)[1] <- "pop.id"
                  tmp
              })
    colnames(rand.admixtures)[1] <- "cancer.id"
    
    rand.admixtures$admixture <- paste0(rand.admixtures$pop.id, rand.admixtures$cancer.id, rand.admixtures$admixture)
    rand.admixtures <- rand.admixtures[, !(colnames(rand.admixtures) %in% c("cancer.id", "pop.id"))]

    ## Confirm all admixtures sum to one
    sums <-
        ddply(rand.admixtures, .variables = c("admixture"),
              .fun = function(df) sum(df$prop))
    eps <- 10^-5
    if(!(all(abs(sums$V1 - 1) < eps))) { stop("Some admixtures don't sum to one\n") }

    rand.admixtures
}

## The following uses variables defined in define-biological-models.R
## Re-define azizi.hierarchical.model and cibersort.hierarchical.model
## to exclude populations we don't have purified expression profiles for
## (i.e., memory_B_cells and eosinophils)

## "Azizi" hierarchical model
## - iAtlas -> purity, leuko, non-leuko stromal
## - non-leuko stromal -> Tirosh -> endo, fibro

## - leuko -> Azizi -> T, B_naive, neutro, DC, mono, macro, NK
## - T -> 10k -> CD4 memory / naive, CD8 memory / naive, Treg

azizi.hierarchical.model <-
    list(
        list(parent = "root",
             children = c("tumor.fraction", "leukocyte.fraction", "non.leukocyte.stromal.fraction"),
             dist = "unif", params = iatlas.leuko.stromal.unif.pan.cancer.fits),
        list(parent = "non.leukocyte.stromal.fraction",
             children = c("Fibroblasts", "Endothelial_cells"),
             dist = "unif", params = tirosh.non.immune.unif.fits),
        list(parent = "leukocyte.fraction",
             children = c("T_cells", "B_cells", "Neutrophils", "Dendritic_cells",
                          "Monocytes", "Macrophages", "NK_cells", "Eosinophils"),
             dist = "unif", params = azizi.leuko.unif.fits),
        list(parent = "T_cells",
             children = c("Memory_CD4_T_cells", "Naive_CD4_T_cells",
                          "Memory_CD8_T_cells", "Naive_CD8_T_cells", "Tregs"),
             dist = "unif", params = cytof.t.cell.unif.fits),
        list(parent = "B_cells",
             children = c("Memory_B_cells", "Naive_B_cells"),
             dist = "unif", params = cytof.b.cell.unif.fits)
    )

## "Cibersort" hierarchical model
## - iAtlas -> purity, leuko, non-leuko stromal
## - non-leuko stromal -> Tirosh -> endo, fibro

## - leuko -> cibersort -> CD4 memory / naive, CD8, Treg, B memory / naive, neutro, DC, mono, macro, NK, eos
## - CD8 -> 10k -> CD8 memory / naive

cibersort.hierarchical.model <-
    list(
        list(parent = "root",
             children = c("tumor.fraction", "leukocyte.fraction", "non.leukocyte.stromal.fraction"),
             dist = "unif", params = iatlas.leuko.stromal.unif.pan.cancer.fits),
        list(parent = "non.leukocyte.stromal.fraction",
             children = c("Fibroblasts", "Endothelial_cells"),
             dist = "unif", params = tirosh.non.immune.unif.fits),
        list(parent = "leukocyte.fraction",
             children = c("Memory_B_cells", "Naive_B_cells", "Memory_CD4_T_cells", "Naive_CD4_T_cells",
                          "NK_cells", "Tregs", "Dendritic_cells", "Monocytes", "Macrophages",
                          "Neutrophils", "Eosinophils", "CD8_T_cells"),
             dist = "unif", params = cibersort.leuko.unif.pan.cancer.fits),
        list(parent = "CD8_T_cells",
             children = c("Memory_CD8_T_cells", "Naive_CD8_T_cells"),
             dist = "unif", params = cytof.cd8.t.cell.unif.fits)
    )


final.pops <- c(populations, "tumor.fraction", "COAD", "READ", "BRCA")

azizi.flat.model <- flatten.hierarchical.unif.model(azizi.hierarchical.model, pops = final.pops)
cibersort.flat.model <- flatten.hierarchical.unif.model(cibersort.hierarchical.model, pops = final.pops)

tmp <- melt(as.matrix(azizi.flat.model))
tmp <- rbind(tmp, melt(as.matrix(cibersort.flat.model)))
colnames(tmp) <- c("row", "col", "value")

tmp <- ddply(tmp, .variables = c("row", "col"),
             .fun = function(df) {
                 if(df$row[1] == "min") {
                     return(df[which.min(df$value),,drop=F])
                 }
                 df[which.max(df$value),,drop=F]                 
             })

flat.model <- acast(tmp, row ~ col)
flat.model["min",] <- unlist(lapply(flat.model["min",], function(x) max(x, 0.01)))
flat.model["max",] <- unlist(lapply(flat.model["max",], function(x) max(x, 0.01)))
non.tumor.cols <- colnames(flat.model)[!(colnames(flat.model) == "tumor.fraction")]
flat.model["max",non.tumor.cols] <- unlist(lapply(flat.model["max",non.tumor.cols], function(x) min(x, 0.5)))

set.seed(1234)
tumor.mins <- seq(from = 0.4, to = 0.6, by = 0.1)
tumor.mins <- seq(from = 0.2, to = 0.6, by = 0.2)
tumor.maxs <- 0.1 + tumor.mins
## tumor.mins <- seq(from = 0.2, to = 0.4, by = 0.1)
admixtures <- llply(1:length(tumor.mins), .parallel = FALSE,
                    .fun = function(i) {
                        tumor.flat.model <- round(flat.model, digits=2)
                        tumor.flat.model["min", "tumor.fraction"] <- tumor.mins[i]
                        tumor.flat.model["max", "tumor.fraction"] <- tumor.maxs[i]
                        sample.flat.unif.model(10, tumor.flat.model, min.prop = 0.01,
                                               constraint.col = "tumor.fraction")
                    })
admixture.tbl <- ldply(admixtures)

fine.grained.pop.maps <- list("A" = fine.grained.pop1.df, "B" = fine.grained.pop2.df)
coarse.grained.pop.maps <- list("A" = coarse.grained.pop1.df, "B" = coarse.grained.pop2.df)

spike.ins <- c(1/(2^seq(from=2,to=12,by=1)), 0)
names(spike.ins) <- spike.ins

fine.grained.rand.admixtures.wo.cancer <-
    generate.random.admixtures.wo.cancer(immune.fine.grained.populations, fine.grained.pop.maps, spike.ins)

coarse.grained.rand.admixtures.wo.cancer <-
    generate.random.admixtures.wo.cancer(immune.coarse.grained.populations, coarse.grained.pop.maps, spike.ins)

cancer.spike.ins <- c(3/4, 1/2, 1/4, 1/8, 1/16, 1/32, 0)
names(cancer.spike.ins) <- cancer.spike.ins

fine.grained.rand.admixtures.w.cancer <-
    generate.random.admixtures.w.cancer(immune.fine.grained.populations, fine.grained.pop.maps, cancer.spike.ins,
                                        cancers = c("CRC", "Breast"))

coarse.grained.rand.admixtures.w.cancer <-
    generate.random.admixtures.w.cancer(immune.coarse.grained.populations, coarse.grained.pop.maps, cancer.spike.ins,
                                        cancer = c("CRC", "Breast"))
