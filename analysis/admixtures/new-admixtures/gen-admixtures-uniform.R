library(ggplot2)
library(dirichlet)
library(plyr)
library(corrplot)
library(reshape2)
library(ggbeeswarm)

## library(devtools)
## install_github("kylebittinger/polyafit")
library(polyafit)
suppressPackageStartupMessages(library("foreach"))
suppressPackageStartupMessages(library("parallel"))

num.cores <- detectCores()
if(!is.na(num.cores) && (num.cores > 1)) {
  suppressPackageStartupMessages(library("doMC"))
  cat(paste("Registering ", num.cores-1, " cores.\n", sep=""))
  registerDoMC(cores=(num.cores-1))
}

## install.packages("dirichlet", repos=c("http://R-Forge.R-project.org"))
library(dirichlet)

library(reshape2)

library(synapser)
synLogin()


## Define the populations of interest
source("populations.R")

source("load-data.R")
source("utils.R")

## Generate n.admixtures admixtures _per_ tumor type (one of which has only that tumor)
n.admixtures <- 10

## Generate admixtures that include one of the two cancer types,
## all at equal proportions
n.populations <- length(populations)
n.populations.excluding.one.tumor.type <- n.populations - 1

equi.proportion <- 1 / n.populations.excluding.one.tumor.type
proportions <- rep(equi.proportion, n.populations)
names(proportions) <- populations

crc.proportions <- proportions
crc.proportions <- proportions[!(names(crc.proportions) == "BRCA")]

brca.proportions <- proportions
brca.proportions <- proportions[!(names(brca.proportions) == "CRC")]

brca.admixtures <- generate.random.uniform.admixtures(names(brca.proportions), n.admixtures, tumor.type = "BRCA", min.prop = 0.01)
crc.admixtures <- generate.random.uniform.admixtures(names(crc.proportions), n.admixtures, tumor.type = "CRC", min.prop = 0.01)

rand.admixtures <- t(rbind.fill(as.data.frame(t(brca.admixtures)), as.data.frame(t(crc.admixtures))))
rand.admixtures[is.na(rand.admixtures)] <- 0

prefix <- "rand-admixture"

file <- paste0(prefix, "-cell-type-correlations.png")
png(file)
plot.admixture.correlations(t(rand.admixtures), main="Random Admixtures")
d <- dev.off()

g.rand.freq <- plot.frequencies(rand.admixtures)
g.rand.freq <- g.rand.freq + geom_hline(yintercept = 0.01, linetype = "dashed")
g.rand.freq <- g.rand.freq + ggtitle("Random Admixtures")

file <- paste0(prefix, "-cell-type-proportions.png")
png(file)
print(g.rand.freq)
d <- dev.off()

li.crc.cnt.df <- load.li.cell.type.counts()

## Can not compare immune (CD45+) and non-immune (CD45-), as these
## were isolated separately.
## In fact, there was further isolation into CD45-/CD90-/CD31- (deplete fibroblasts and endothelial and enrich for malignant)
## and CD45+/CD3+ (to enrich for T cells).
## Hence, really not use to quantitate relative abundance.
puram.hn.cnt.df <- load.puram.cell.type.counts()
## Subset to Primary samples
flag <- grepl(rownames(puram.hn.cnt.df), pattern="-P")
puram.hn.cnt.df <- puram.hn.cnt.df[flag, ]
puram.hn.cnt.df[is.na(puram.hn.cnt.df)] <- 0
rowSums(puram.hn.cnt.df)

iatlas.tbl <- load.iAtlas.supp()
## The last few columns "Lymphocytes", "Neutrophils", "Eosinophils", "Mast.Cells", "Dendritic.Cells", and "Macrophages"
## are an aggregation of earlier columns--note that "Neutrophils" and "Eosinophils" occur multiple times.
cols.to.keep <- 1:(which(colnames(iatlas.tbl) == "Lymphocytes")-1)

iatlas.tbl <- iatlas.tbl[, cols.to.keep]
iatlas.tbl$tumor.fraction = 1 - iatlas.tbl$Stromal.Fraction
iatlas.tbl$leukocyte.fraction = iatlas.tbl$Leukocyte.Fraction
iatlas.tbl$non.leukocyte.stromal.fraction <- 1 - (iatlas.tbl$tumor.fraction + iatlas.tbl$leukocyte.fraction)

## cibersort columns
cs.cols <- c(
    "B.Cells.Memory",                         
    "B.Cells.Naive",                          
    "Dendritic.Cells.Activated",              
    "Dendritic.Cells.Resting",                
    "Eosinophils",                            
    "Macrophages.M0",                         
    "Macrophages.M1",                         
    "Macrophages.M2",                         
    "Mast.Cells.Activated",                   
    "Mast.Cells.Resting",                     
    "Monocytes",                              
    "Neutrophils",                            
    "NK.Cells.Activated",                     
    "NK.Cells.Resting",                       
    "Plasma.Cells",                           
    "T.Cells.CD4.Memory.Activated",           
    "T.Cells.CD4.Memory.Resting",             
    "T.Cells.CD4.Naive",                      
    "T.Cells.CD8",                            
    "T.Cells.Follicular.Helper",              
    "T.Cells.gamma.delta",                    
    "T.Cells.Regulatory.Tregs"
)

## Define new aggregated cibersort tables
iatlas.tbl$Memory_CD4_T_cells <- as.numeric(iatlas.tbl$T.Cells.CD4.Memory.Activated) +
    as.numeric(iatlas.tbl$T.Cells.CD4.Memory.Resting)
iatlas.tbl$NK_cells <- as.numeric(iatlas.tbl$NK.Cells.Activated) +
    as.numeric(iatlas.tbl$NK.Cells.Resting)
iatlas.tbl$Dendritic_cells <- as.numeric(iatlas.tbl$Dendritic.Cells.Activated) +
    as.numeric(iatlas.tbl$Dendritic.Cells.Resting)
iatlas.tbl$Macrophages <- as.numeric(iatlas.tbl$Macrophages.M0) +
    as.numeric(iatlas.tbl$Macrophages.M1) +
    as.numeric(iatlas.tbl$Macrophages.M2)

old.cs.leuko.cols <- c("B.Cells.Memory", "B.Cells.Naive", "Memory_CD4_T_cells", "T.Cells.CD4.Naive", "NK_cells",
                       "T.Cells.Regulatory.Tregs", "Dendritic_cells", "Monocytes", "Macrophages",
                       "Neutrophils", "Eosinophils", "T.Cells.CD8")
new.cs.leuko.cols <- c("Memory_B_cells", "Naive_B_cells", "Memory_CD4_T_cells", "Naive_CD4_T_cells", "NK_cells",
                       "Tregs", "Dendritic_cells", "Monocytes", "Macrophages",
                       "Neutrophils", "Eosinophils", "CD8_T_cells")


iatlas.tbls <- dlply(iatlas.tbl, .variables = "TCGA.Study", .fun = function(df) df)


iatlas.leuko.stromal.purity.tbls <-
    llply(iatlas.tbls,
          .fun = function(df) {
              ret <- data.frame(tumor.fraction = 1 - df$`Stromal.Fraction`, leukocyte.fraction = df$Leukocyte.Fraction)
              ret$non.leukocyte.stromal.fraction <- 1 - (ret$tumor.fraction + ret$leukocyte.fraction)
              ret
          })

## Can not compare immune (CD45+) and non-immune (CD45-), as these
## were isolated separately.
## Hence, can not compare relative abundance of immune vs non-immune.
## In four tumors (Mel58, 67, 72 and 74), we sequenced primarily the immune infiltrates (CD45+ cells) and there were only zero or one malignant cells by this definition.
## This seems also to be true for CY75, which doesn't have any tumor cells.
tirosh.melanoma.cnt.df <- load.tirosh.cell.type.counts()

tirosh.non.immune.cnt.orig.df <- tirosh.melanoma.cnt.df[, c("CAF", "Endo.", "melanoma")]
exclude <- c("CY58", "CY67", "CY72", "CY74", "CY75")
tirosh.non.immune.cnt.orig.df <- tirosh.non.immune.cnt.orig.df[!(rownames(tirosh.non.immune.cnt.orig.df) %in% exclude),]
tirosh.non.immune.cnt.orig.df[is.na(tirosh.non.immune.cnt.orig.df)] <- 0

old.cols <- c("CAF", "Endo.", "melanoma")
new.cols <- c("Fibroblasts", "Endothelial_cells", "Cancer")
##old.cols <- c("CAF", "Endo.")
##new.cols <- c("Fibroblasts", "Endothelial_cells")
tirosh.non.immune.cnt.df <- tirosh.non.immune.cnt.orig.df[, old.cols]
colnames(tirosh.non.immune.cnt.df) <- new.cols

azizi.cnt.df <- load.azizi.cell.type.counts()
old.cols <- c("T.cell", "B.cell", "Neutrophil", "DC", "Monocyte", "Macrophage", "NK.cell", "Mast.cell")
## NB: We're going to treat eosinophils as noise.  Though Azizi doesn't assay them,
## let's call Mast cells Eosinophils.
new.cols <- c("T_cells", "B_cells", "Neutrophils", "Dendritic_cells",
              "Monocytes", "Macrophages", "NK_cells", "Eosinophils")
azizi.cnt.df <- azizi.cnt.df[, old.cols]
colnames(azizi.cnt.df) <- new.cols

cytof.frac.df <- load.cytof.10k.cell.type.fractions()
cols <- colnames(cytof.frac.df)
tmp <- cytof.frac.df[, cols]
tmp <- tmp / rowSums(tmp)
tmp <- na.omit(tmp)
cytof.frac.df <- tmp


## Get min/max for tumor, leukocyte, non-leukocyte stromal populations
## from iAtlas data

collapse.tumor.min.max <- function(tbls) {

    tbl <- ldply(tbls, .fun = function(df) df)
    colnames(tbl) <- c("tumor", "population", "min", "max")
    
    pan.cancer.tbl <-
        ddply(tbl,
              .variables = c("population"),
              .fun = function(df) {
                  data.frame(min = min(df$min), max = max(df$max))
              })
    pan.cancer.tbl
}

calc.col.min.max <- function(df, min.prob = 0, max.prob = 1) {
    cols <- colnames(df)
    names(cols) <- cols
    ldply(cols,
          .fun = function(col) {
              data.frame(min = as.numeric(quantile(na.omit(df[, col]), probs = c(min.prob))),
                         max = as.numeric(quantile(na.omit(df[, col]), probs = c(max.prob))))
          })
}

clean.up.min.max.params <- function(df) {
    flag <- df$min < 0
    df$min[flag] <- 0
    flag <- df$max > 1
    df$max[flag] <- 1
    df
}

invert.min.max.params <- function(df) {
    pop.col <- c("population", ".id")
    pop.col <- pop.col[pop.col %in% colnames(df)][1]
    pops <- df[, pop.col]
    df <- df[, !(colnames(df) %in% pop.col)]
    df <- t(df)
    colnames(df) <- pops
    df
}
    
min.prob <- pnorm(-2)
max.prob <- pnorm(2)

iatlas.leuko.stromal.unif.fits <-
    llply(iatlas.leuko.stromal.purity.tbls,
          .fun = function(df) {
              cols <- colnames(df)
              names(cols) <- cols
              df <- na.omit(df)
              if(nrow(df) == 0) { return(NULL) }
              ret <- ldply(cols,
                           .fun = function(col) {
                               data.frame(min = as.numeric(quantile(df[, col], probs=c(min.prob))),
                                          max = as.numeric(quantile(df[, col], probs=c(max.prob))))
                           })
              colnames(ret) <- c("population", "min", "max")
              ret
          })

##iatlas.leuko.stromal.unif.pan.cancer.orig.fits <-
##    collapse.tumor.min.max(iatlas.leuko.stromal.unif.fits)

tmp <- ldply(iatlas.leuko.stromal.purity.tbls)
tmp <- tmp[, !(colnames(tmp) == ".id")]
iatlas.leuko.stromal.unif.pan.cancer.orig.fits <- calc.col.min.max(tmp, min.prob = min.prob, max.prob = max.prob)


## Get min/max for leukocyte populations from CIBERSORT
cibersort.leuko.unif.fits <-
    llply(iatlas.tbls,
          .fun = function(df) {
              df <- df[, old.cs.leuko.cols]
              colnames(df) <- new.cs.leuko.cols
              cols <- colnames(df)
              names(cols) <- cols
              df <- na.omit(df)
              if(nrow(df) == 0) { return(NULL) }
              ret <- ldply(cols,
                           .fun = function(col) {
                               data.frame(min = as.numeric(quantile(na.omit(df[, col]), probs = c(min.prob))),
                                          max = as.numeric(quantile(na.omit(df[, col]), probs = c(max.prob))))
                           })
              colnames(ret) <- c("population", "min", "max")
              ret
          })

##cibersort.leuko.unif.pan.cancer.orig.fits <-
##    collapse.tumor.min.max(cibersort.leuko.unif.orig.fits)

tmp <- ldply(iatlas.tbls)
tmp <- tmp[, !(colnames(tmp) == ".id")]
tmp <- tmp[, old.cs.leuko.cols]
colnames(tmp) <- new.cs.leuko.cols
cibersort.leuko.unif.pan.cancer.orig.fits <- calc.col.min.max(tmp, min.prob = min.prob, max.prob = max.prob)

tirosh.non.immune.frac.df <- tirosh.non.immune.cnt.df /
    rowSums(tirosh.non.immune.cnt.df)
rowSums(tirosh.non.immune.frac.df)

tirosh.non.immune.unif.orig.fits <- calc.col.min.max(tirosh.non.immune.frac.df, min.prob = min.prob,
                                                     max.prob = max.prob)
colnames(tirosh.non.immune.unif.orig.fits) <- c("population", "min", "max")

azizi.leuko.frac.df <-
    azizi.cnt.df / rowSums(azizi.cnt.df)

azizi.leuko.unif.orig.fits <- calc.col.min.max(azizi.leuko.frac.df, min.prob = min.prob, max.prob = max.prob)
colnames(azizi.leuko.unif.orig.fits) <- c("population", "min", "max")

cytof.leuko.unif.orig.fits <- calc.col.min.max(cytof.frac.df, min.prob = min.prob, max.prob = max.prob)
colnames(cytof.leuko.unif.orig.fits) <- c("population", "min", "max")

cols <- c("Tregs", "Naive_CD4_T_cells", "Memory_CD4_T_cells",
          "Naive_CD8_T_cells", "Memory_CD8_T_cells")
tmp <- cytof.frac.df[, cols]
cytof.t.cell.unif.orig.fits <- calc.col.min.max(na.omit(tmp / rowSums(tmp)), min.prob = min.prob,
                                                max.prob = max.prob)
colnames(cytof.t.cell.unif.orig.fits) <- c("population", "min", "max")

cols <- c("Naive_CD8_T_cells", "Memory_CD8_T_cells")
tmp <- cytof.frac.df[, cols]
cytof.cd8.t.cell.unif.orig.fits <- calc.col.min.max(na.omit(tmp / rowSums(tmp)), min.prob = min.prob,
                                                    max.prob = max.prob)
colnames(cytof.cd8.t.cell.unif.orig.fits) <- c("population", "min", "max")

cols <- c("Naive_B_cells", "Memory_B_cells")
tmp <- cytof.frac.df[, cols]
cytof.b.cell.unif.orig.fits <- calc.col.min.max(na.omit(tmp / rowSums(tmp)), min.prob = min.prob,
                                           max.prob = max.prob)
colnames(cytof.b.cell.unif.orig.fits) <- c("population", "min", "max")

azizi.leuko.unif.fits <- clean.up.min.max.params(azizi.leuko.unif.orig.fits)
cytof.t.cell.unif.fits <- clean.up.min.max.params(cytof.t.cell.unif.orig.fits)
cytof.cd8.t.cell.unif.fits <- clean.up.min.max.params(cytof.cd8.t.cell.unif.orig.fits)
cytof.b.cell.unif.fits <- clean.up.min.max.params(cytof.b.cell.unif.orig.fits)
cytof.leuko.unif.fits <- clean.up.min.max.params(cytof.leuko.unif.orig.fits)
tirosh.non.immune.unif.fits <- clean.up.min.max.params(tirosh.non.immune.unif.orig.fits)
cibersort.leuko.unif.pan.cancer.fits  <- clean.up.min.max.params(cibersort.leuko.unif.pan.cancer.orig.fits )
iatlas.leuko.stromal.unif.pan.cancer.fits <- clean.up.min.max.params(iatlas.leuko.stromal.unif.pan.cancer.orig.fits)

azizi.leuko.unif.fits <- invert.min.max.params(azizi.leuko.unif.fits)
cytof.t.cell.unif.fits <- invert.min.max.params(cytof.t.cell.unif.fits)
cytof.cd8.t.cell.unif.fits <- invert.min.max.params(cytof.cd8.t.cell.unif.fits)
cytof.b.cell.unif.fits <- invert.min.max.params(cytof.b.cell.unif.fits)
cytof.leuko.unif.fits <- invert.min.max.params(cytof.leuko.unif.fits)
tirosh.non.immune.unif.fits <- invert.min.max.params(tirosh.non.immune.unif.fits)
cibersort.leuko.unif.pan.cancer.fits  <- invert.min.max.params(cibersort.leuko.unif.pan.cancer.fits )
iatlas.leuko.stromal.unif.pan.cancer.fits <- invert.min.max.params(iatlas.leuko.stromal.unif.pan.cancer.fits)


## "Azizi" hierarchical model
## - iAtlas -> purity, leuko, non-leuko stromal
## - non-leuko stromal -> Tirosh -> endo, fibro

## - leuko -> Azizi -> T, B, neutro, DC, mono, macro, NK, mast
## - T -> 10k -> CD4 memory / naive, CD8 memory / naive, Treg
## - B -> 10k -> memory / naive
## - mast -> eosinophils

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

set.seed(1234)
n.admixtures <- 10

foo <- sample.hierarchical.model(n.admixtures, azizi.hierarchical.model, pops = final.pops, min.prop = 0.01)

bar <- sample.hierarchical.model(n.admixtures, cibersort.hierarchical.model, pops = final.pops, min.prop = 0.01)

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

q(status = 0)

set.seed(1234)

tumor.mins <- seq(from = 0.1, to = 0.65, by = 0.05)
admixtures <- llply(tumor.mins,
                    .fun = function(tumor.min) {
                        set.seed(1234)
                        tumor.flat.model <- round(flat.model, digits=2)
                        tumor.flat.model["min", "tumor.fraction"] <- tumor.min
                        tumor.flat.model["max", "tumor.fraction"] <- tumor.min + 0.05
                        print(tumor.min)
                        sample.flat.unif.model(5, tumor.flat.model, min.prop = 0.01,
                                               constraint.col = "tumor.fraction")
                    })

save.image(".Rdata")

broken.stick.admixtures <- generate.random.uniform.admixtures(colnames(flat.model), 50,
                                                              tumor.type = "tumor.fraction",
                                                              min.prop = 0.01)
broken.stick.admixtures <- t(broken.stick.admixtures)
colnames(broken.stick.admixtures) <- colnames(flat.model)

admixture.tbl <- ldply(admixtures)
wb <- createWorkbook("Admixtures")
addWorksheet(wb, "model")
writeData(wb, sheet = 1, flat.model, colNames = TRUE, rowNames = TRUE)
addWorksheet(wb, "uniform-admixtures")
writeData(wb, sheet = 2, admixture.tbl, colNames = TRUE, rowNames = FALSE)
addWorksheet(wb, "broken-stick-admixtures")
writeData(wb, sheet = 3, broken.stick.admixtures, colNames = TRUE, rowNames = FALSE)
saveWorkbook(wb, "031419-admixtures.xlsx", overwrite = TRUE)

prefix <- "uniform"
file <- paste0(prefix, "-cell-type-correlations.png")
png(file)
plot.admixture.correlations(admixture.tbl, main="Uniform cell type correlations")
d <- dev.off()

prefix <- "uniform"
file <- paste0(prefix, "-sample-correlations.png")
png(file)
plot.admixture.correlations(t(admixture.tbl), main="Uniform admixture correlations")
d <- dev.off()

prefix <- "broken-stick"
file <- paste0(prefix, "-cell-type-correlations.png")
png(file)
plot.admixture.correlations(broken.stick.admixtures, main="Broken stick cell type correlations")
d <- dev.off()

prefix <- "broken-stick"
file <- paste0(prefix, "-sample-correlations.png")
png(file)
plot.admixture.correlations(t(broken.stick.admixtures), main="Broken stick admixture correlations")
d <- dev.off()


save.image(".Rdata")

q(status=0)

tumor.mins2 <- seq(from = 0.55, to = 0.6, by = 0.1)
admixtures2 <- llply(tumor.mins2,
                    .fun = function(tumor.min) {
                        set.seed(1234)
                        tumor.flat.model <- round(flat.model, digits=2)
                        tumor.flat.model["min", "tumor.fraction"] <- tumor.min
                        tumor.flat.model["max", "tumor.fraction"] <- tumor.min + 0.05
                        print(tumor.min)
                        sample.flat.unif.model(5, tumor.flat.model, min.prop = 0.01,
                                               constraint.col = "tumor.fraction")
                    })

save.image(".Rdata")

q(status = 0)

admixtures.tumor.low <- sample.flat.unif.model(20, flat.model, min.prop = 0.01,
                                               constraint.col = "tumor.fraction")

set.seed(1234)
tumor.mid.flat.model <- flat.model
tumor.mid.flat.model["min", "tumor.fraction"] <- 0.4
admixtures.tumor.mid <- sample.flat.unif.model(20, tumor.mid.flat.model, min.prop = 0.01,
                                               constraint.col = "tumor.fraction")

set.seed(1234)
tumor.high.flat.model <- flat.model
tumor.high.flat.model["min", "tumor.fraction"] <- 0.6
admixtures.tumor.high <- sample.flat.unif.model(20, tumor.high.flat.model, min.prop = 0.01,
                                                constraint.col = "tumor.fraction")

save.image(".Rdata")

model.names <- names(azizi.models)
##model.names <- c("coad-cytof")
for(model.name in model.names) {
    prefix <- paste0(model.name, "-admixture")

    file <- paste0(prefix, "-cell-type-correlations.png")
    png(file)
    plot.admixture.correlations(t(azizi.admixtures[[model.name]]), main=model.name)
    d <- dev.off()

    g.rand.freq <- plot.frequencies(azizi.admixtures[[model.name]])
    g.rand.freq <- g.rand.freq + geom_hline(yintercept = 0.01, linetype = "dashed")
    g.rand.freq <- g.rand.freq + ggtitle(paste0(model.name, " Admixtures"))
    
    file <- paste0(prefix, "-cell-type-proportions.png")
    png(file)
    print(g.rand.freq)
    d <- dev.off()
}

set.seed(1234)
n.admixtures <- 10
##brca.cytof.admixtures <- sample.hierarchical.model(n.admixtures, hierarchical.model, pops = final.pops, min.prop = 0.001)
##coad.cytof.admixtures <- sample.hierarchical.model(n.admixtures, models[["coad-cytof"]], pops = final.pops, min.prop = 0.01)
##admixtures[["coad-cytof"]] <- coad.cytof.admixtures
##cibersort.admixtures <- llply(cibersort.models, .parallel = TRUE,
##                    .fun = function(model) sample.hierarchical.model(n.admixtures, model, pops = final.pops, min.prop = 10^-3))

model.names <- names(cibersort.models)
## model.names <- c("brca-cibersort", "coad-cibersort")
cibersort.admixtures <- list()
##model.names <- c("coad-cytof")
## cibersort.min.prop <- 10^-3
cibersort.min.prop <- 10^-4
for(model.name in model.names) {
    prefix <- paste0(model.name, "-admixture")

    cibersort.admixtures[[model.name]] <-
        sample.hierarchical.model(n.admixtures, cibersort.models[[model.name]], pops = final.pops, min.prop = cibersort.min.prop)

    save.image(".Rdata")
    
    file <- paste0(prefix, "-cell-type-correlations.png")
    png(file)
    plot.admixture.correlations(t(cibersort.admixtures[[model.name]]), main=model.name)
    d <- dev.off()

    g.rand.freq <- plot.frequencies(cibersort.admixtures[[model.name]])
    g.rand.freq <- g.rand.freq + geom_hline(yintercept = 0.01, linetype = "dashed")
    g.rand.freq <- g.rand.freq + ggtitle(paste0(model.name, " Admixtures"))
    
    file <- paste0(prefix, "-cell-type-proportions.png")
    png(file)
    print(g.rand.freq)
    d <- dev.off()
}

## TODO
## - get min/max for
##   - textbook
