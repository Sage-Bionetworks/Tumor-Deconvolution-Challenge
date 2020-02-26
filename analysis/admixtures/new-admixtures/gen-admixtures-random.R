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

unif.broken.stick.proportion <- function(n, min.prop = 0) {
    sq <- seq(from = min.prop, to = 1 - min.prop, by = min.prop)
    sam <- sample(sq, size = n-1, replace = FALSE)
    sam <- c(0,sort(sam), 1)
    ret <- unlist(lapply(1:(length(sam)-1), function(i) sam[i+1]-sam[i]))
    if(sum(ret) != 1) { stop("unif.broken.stick.proportion sum is not 1\n") }
    ret
}

generate.random.uniform.admixtures <- function(populations, n, tumor.type, min.prop = 0) {

    ## Generate the random admixtures
    mat <- ldply(1:(n-1), .fun = function(i) unif.broken.stick.proportion(length(populations), min.prop = min.prop))
    print(dim(mat))
    ## Generate the admixture with only tumor content
    cancer.only.admixture <- rep(0, length(populations))
    cancer.only.admixture[populations == tumor.type] <- 1
    mat <- rbind(mat, cancer.only.admixture)
    mat <- t(mat)
    rownames(mat) <- populations
    colnames(mat) <- NULL
    mat
}

generate.random.dirichlet.admixtures <- function(proportions, n, tumor.type, min.prop = 0) {

    ## Generate the random admixtures
    mat <- rdirichlet(n=n-1, alpha=proportions*5)
    flag <- unlist(apply(mat, 1, function(row) any(row < min.prop)))
    while(any(flag)) {
        mat[flag,] <- rdirichlet(n=length(which(flag)), alpha=proportions*5)
        flag <- unlist(apply(mat, 1, function(row) any(row < min.prop)))
    }
    
    ## Generate the admixture with only tumor content
    cancer.only.admixture <- rep(0, length(proportions))
    cancer.only.admixture[names(proportions) == tumor.type] <- 1

    mat <- rbind(mat, cancer.only.admixture)
    mat <- t(mat)
    rownames(mat) <- names(proportions)
    colnames(mat) <- NULL
    mat
}

plot.admixture.correlations <- function(admixtures, ...) {

    fc <- as.matrix(cor(admixtures, method = "spearman"))
    
    mar <- c(1, 1, 1, 1)
    ##    corrplot(fc, method = "ellipse", type = "upper", order = "hclust", tl.cex = 0.6, mar = mar, ...)
    corrplot(fc, method = "ellipse", type = "upper", order = "original", tl.cex = 0.6, mar = mar, ...)
}

plot.frequencies <- function(admixtures) {
    m <- melt(admixtures)
    colnames(m) <- c("Population", "Admixture", "Proportion")
    g <- ggplot()
    m$Population <- factor(m$Population, levels = rownames(admixtures))
    ##    g <- g + geom_point(data = m, aes(x = Population, y = Proportion))
    ## g <- g + geom_beeswarm(data = m, aes(x = Population, y = log(Proportion)))
    ## g <- g + geom_violin(data = m, aes(x = Population, y = Proportion))
##    g <- g + geom_boxplot(data = m, aes(x = Population, y = log(Proportion)))
##    g <- g + geom_beeswarm(data = m, aes(x = Population, y = Proportion))    
    g <- g + geom_boxplot(data = m, aes(x = Population, y = Proportion))    
    g <- g + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    g
}

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

## Fit dirichlets to the tumor, leukocyte, non-leukocyte stromal fractions
## (excluding the few that are negative)
## Dirichlet parameters are k * p = most.likely.k * p
iatlas.leuko.stromal.dir.fits <-
    llply(iatlas.leuko.stromal.purity.tbls,
          .fun = function(df) {
              df <- na.omit(df)
              mat <- as.matrix(df)
##              colnames(mat) <- NULL
              fit <- fit.dirichlet(mat, type="mm")
              if(is.na(fit$weighted.k)) { return(NA) }
              alpha <- unname(fit$most.likely.k) * unname(fit$p)
              names(alpha) <- names(fit$p)
              alpha
          })


## Fit dirichlets to the tumor, leukocyte, non-leukocyte stromal fractions
## (excluding the few that are negative)
## Dirichlet parameters are k * p = most.likely.k * p
cibersort.leuko.dir.fits <-
    llply(iatlas.tbls,
          .fun = function(df) {
              df <- df[, old.cs.leuko.cols]
              colnames(df) <- new.cs.leuko.cols
              df <- na.omit(df)
              mat <- as.matrix(df)
##              colnames(mat) <- NULL
              fit <- fit.dirichlet(mat, type="mm")
              if(is.na(fit$weighted.k)) { return(NA) }
              alpha <- unname(fit$most.likely.k) * unname(fit$p)
              names(alpha) <- names(fit$p)
              alpha
          })

## Fit multivariate Polya (multinomial-Dirichlet) to CAF/Endothelial/Cancer in Tirosh melanoma

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
flag <- rowSums(tirosh.non.immune.cnt.df) > 0
tirosh.non.immune.polya.fits <- optim_polya(tirosh.non.immune.cnt.df[flag,])
tirosh.non.immune.polya.fits$par / sum(tirosh.non.immune.polya.fits$par)
colMeans(tirosh.non.immune.cnt.df / rowSums(tirosh.non.immune.cnt.df))
rdirichlet(1, tirosh.non.immune.polya.fits$par)

azizi.cnt.df <- load.azizi.cell.type.counts()
old.cols <- c("T.cell", "B.cell", "Neutrophil", "DC", "Monocyte", "Macrophage", "NK.cell", "Mast.cell")
## NB: We're going to treat eosinophils as noise.  Though Azizi doesn't assay them,
## let's call Mast cells Eosinophils.
new.cols <- c("T_cells", "B_cells", "Neutrophils", "Dendritic_cells",
              "Monocytes", "Macrophages", "NK_cells", "Eosinophils")
azizi.cnt.df <- azizi.cnt.df[, old.cols]
colnames(azizi.cnt.df) <- new.cols

azizi.polya.fits <- optim_polya(azizi.cnt.df)
azizi.polya.fits$par / sum(azizi.polya.fits$par)
colMeans(azizi.cnt.df / rowSums(azizi.cnt.df))
rdirichlet(1, azizi.polya.fits$par)

cytof.frac.df <- load.cytof.10k.cell.type.fractions()
cols <- c("Memory_CD4_T_cells", "Naive_CD4_T_cells",
          "Memory_CD8_T_cells", "Naive_CD8_T_cells", "Tregs")
tmp <- cytof.frac.df[, cols]
tmp <- tmp / rowSums(tmp)
tmp <- na.omit(tmp)
fit <- fit.dirichlet(as.matrix(tmp), type="mm")
alpha <- unname(fit$most.likely.k) * unname(fit$p)
names(alpha) <- names(fit$p)
cytof.t.cell.params <- alpha

cols <- c("Memory_CD8_T_cells", "Naive_CD8_T_cells")
tmp <- cytof.frac.df[, cols]
tmp <- tmp / rowSums(tmp)
tmp <- na.omit(tmp)
fit <- fit.dirichlet(as.matrix(tmp), type="mm")
alpha <- unname(fit$most.likely.k) * unname(fit$p)
names(alpha) <- names(fit$p)
cytof.cd8.t.cell.params <- alpha

cols <- c("Memory_B_cells", "Naive_B_cells")
tmp <- cytof.frac.df[, cols]
tmp <- tmp / rowSums(tmp)
tmp <- na.omit(tmp)
fit <- fit.dirichlet(as.matrix(tmp), type="mm")
alpha <- unname(fit$most.likely.k) * unname(fit$p)
names(alpha) <- names(fit$p)
cytof.b.cell.params <- alpha

## Compare iAtlas and Tirosh estimates of cancer vs non-leuko stromal
df <- data.frame(method = "iAtlas", tumor.type = iatlas.tbl$TCGA.Study,
                 tumor.vs.non.leuko.stromal.ratio = iatlas.tbl$tumor.fraction / iatlas.tbl$non.leukocyte.stromal.fraction)
ratio <- tirosh.non.immune.polya.fits$par[["Cancer"]] /
    (tirosh.non.immune.polya.fits$par[["Fibroblasts"]] +
     tirosh.non.immune.polya.fits$par[["Endothelial_cells"]])
df2 <- data.frame(method = "Tirosh", tumor.type = "Melanoma",
                  tumor.vs.non.leuko.stromal.ratio = ratio)
df <- rbind(df, df2)
df <- na.omit(df)
flag <- df$tumor.vs.non.leuko.stromal.ratio >= 0
df <- df[flag,]
## g <- ggplot() + geom_point(data = df, aes(x = tumor.type, y = tumor.vs.non.leuko.stromal.ratio))
## g <- g + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
##g <- ggplot() + geom_histogram(data = subset(df, tumor.type == "SKCM" & tumor.vs.non.leuko.stromal.ratio < 10),
##                             aes(x = tumor.vs.non.leuko.stromal.ratio))
g <- ggplot() + geom_density(data = subset(df, tumor.type == "SKCM" & tumor.vs.non.leuko.stromal.ratio < 10),
                             aes(x = tumor.vs.non.leuko.stromal.ratio))
g <- g + geom_vline(xintercept = ratio)

old.cols <- c("CAF", "Endo.")
new.cols <- c("Fibroblasts", "Endothelial_cells")
tirosh.non.immune.cnt.df <- tirosh.non.immune.cnt.orig.df[, old.cols]
colnames(tirosh.non.immune.cnt.df) <- new.cols
flag <- rowSums(tirosh.non.immune.cnt.df) > 0
tirosh.non.immune.polya.fits <- optim_polya(tirosh.non.immune.cnt.df[flag,])
tirosh.non.immune.polya.fits$par / sum(tirosh.non.immune.polya.fits$par)
colMeans(tirosh.non.immune.cnt.df[flag,] / rowSums(tirosh.non.immune.cnt.df[flag,]))
rdirichlet(1, tirosh.non.immune.polya.fits$par)

brca.purity.params <- iatlas.leuko.stromal.dir.fits[["BRCA"]]
names(brca.purity.params)[names(brca.purity.params) == "tumor.fraction"] <- "BRCA"

coad.purity.params <- iatlas.leuko.stromal.dir.fits[["COAD"]]
names(coad.purity.params)[names(coad.purity.params) == "tumor.fraction"] <- "COAD"

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
             children = c("BRCA", "leukocyte.fraction", "non.leukocyte.stromal.fraction"),
             dist = "dir", params = brca.purity.params),
        list(parent = "non.leukocyte.stromal.fraction",
             children = c("Fibroblasts", "Endothelial_cells"),
             dist = "dir", params = tirosh.non.immune.polya.fits$par),
        list(parent = "leukocyte.fraction",
             children = c("T_cells", "B_cells", "Neutrophils", "Dendritic_cells",
                          "Monocytes", "Macrophages", "NK_cells", "Eosinophils"),
             dist = "dir", params = azizi.polya.fits$par),
        list(parent = "T_cells",
             children = c("Memory_CD4_T_cells", "Naive_CD4_T_cells",
                          "Memory_CD8_T_cells", "Naive_CD8_T_cells", "Tregs"),
             dist = "dir", params = cytof.t.cell.params),
        list(parent = "B_cells",
             children = c("Memory_B_cells", "Naive_B_cells"),
             dist = "dir", params = cytof.b.cell.params)
    )

## "Cibersort" hierarchical model
## - iAtlas -> purity, leuko, non-leuko stromal
## - non-leuko stromal -> Tirosh -> endo, fibro

## - leuko -> cibersort -> CD4 memory / naive, CD8, Treg, B memory / naive, neutro, DC, mono, macro, NK, eos
## - CD8 -> 10k -> CD8 memory / naive

cibersort.hierarchical.model <-
    list(
        list(parent = "root",
             children = c("BRCA", "leukocyte.fraction", "non.leukocyte.stromal.fraction"),
             dist = "dir", params = brca.purity.params),
        list(parent = "non.leukocyte.stromal.fraction",
             children = c("Fibroblasts", "Endothelial_cells"),
             dist = "dir", params = tirosh.non.immune.polya.fits$par),
        list(parent = "leukocyte.fraction",
             children = c("Memory_B_cells", "Naive_B_cells", "Memory_CD4_T_cells", "Naive_CD4_T_cells",
                          "NK_cells", "Tregs", "Dendritic_cells", "Monocytes", "Macrophages",
                          "Neutrophils", "Eosinophils", "CD8_T_cells"),
             dist = "dir", params = cibersort.leuko.dir.fits[["BRCA"]]),
        list(parent = "CD8_T_cells",
             children = c("Memory_CD8_T_cells", "Naive_CD8_T_cells"),
             dist = "dir", params = cytof.cd8.t.cell.params)
    )

sample.dist <- function(n, dist.name, params) {
    ret <- NA
    if(dist.name == "dir") {
        ret <- rdirichlet(n, params)
    }
    return(ret)
}

sample.hierarchical.model_ <- function(n, hierarchical.model) {
    df <- data.frame(root = rep(1, n))
    l <- hierarchical.model
    while(length(l) > 0) {
        for(i in 1:length(l)) {
            if(!l[[i]]$parent %in% colnames(df)) { next }
            mat <- sample.dist(n, l[[i]]$dist, l[[i]]$params)
            mat <- mat * df[, l[[i]]$parent]
            colnames(mat) <- names(l[[i]]$params)
            df <- cbind(df, mat)
            l[[i]] <- NULL
            break
        }
    }
    df
}

sample.hierarchical.model <- function(n, hierarchical.model, pops, min.prop = 0.01) {
    df <- sample.hierarchical.model_(n, hierarchical.model)
    df <- df[, colnames(df) %in% pops]
    df <- df / rowSums(df)

    flag <- unlist(apply(df, 1, function(row) any(row < min.prop)))
    while(any(flag)) {
        print(length(which(flag)))
        tmp <- sample.hierarchical.model_(length(which(flag)), hierarchical.model)
        tmp <- tmp[, colnames(tmp) %in% pops]
        tmp <- tmp / rowSums(tmp)
        df[flag,] <- tmp
        flag <- unlist(apply(df, 1, function(row) any(row < min.prop)))
    }
    t(df)
}


azizi.models <- list("coad-azizi" = azizi.hierarchical.model, "brca-azizi" = azizi.hierarchical.model)
cibersort.models <- list("coad-cibersort" = cibersort.hierarchical.model, "brca-cibersort" = cibersort.hierarchical.model)
azizi.models[["coad-azizi"]][[1]] <- list(parent = "root",
                                          children = c("COAD", "leukocyte.fraction", "non.leukocyte.stromal.fraction"),
                                          dist = "dir", params = coad.purity.params)
cibersort.models[["coad-cibersort"]][[1]] <- list(parent = "root",
                                                  children = c("COAD", "leukocyte.fraction", "non.leukocyte.stromal.fraction"),
                                                  dist = "dir", params = coad.purity.params)
cibersort.models[["coad-cibersort"]][[3]] <-
    list(parent = "leukocyte.fraction",
         children = c("Memory_B_cells", "Naive_B_cells", "Memory_CD4_T_cells", "Naive_CD4_T_cells",
                      "NK_cells", "Tregs", "Dendritic_cells", "Monocytes", "Macrophages",
                      "Neutrophils", "Eosinophils", "CD8_T_cells"),
         dist = "dir", params = cibersort.leuko.dir.fits[["COAD"]])


final.pops <- c(populations, "tumor.fraction", "COAD", "READ", "BRCA")

set.seed(1234)
n.admixtures <- 10
admixtures <- list()
##brca.cytof.admixtures <- sample.hierarchical.model(n.admixtures, hierarchical.model, pops = final.pops, min.prop = 0.001)
##coad.cytof.admixtures <- sample.hierarchical.model(n.admixtures, models[["coad-cytof"]], pops = final.pops, min.prop = 0.01)
##admixtures[["coad-cytof"]] <- coad.cytof.admixtures
## azizi.min.prop <- 10^-2
azizi.min.prop <- 10^-3
azizi.admixtures <- llply(azizi.models, .parallel = TRUE,
                    .fun = function(model) sample.hierarchical.model(n.admixtures, model, pops = final.pops, min.prop = azizi.min.prop))

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

