usePackage <- function(p) 
{
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, repos = "http://cran.us.r-project.org", dep = TRUE)
  require(p, character.only = TRUE)
}

usePackage("pacman")

suppressPackageStartupMessages(p_load("GEOquery"))

home_dir <- "../../../../Tumor-Deconvolution-Challenge/"
cur_dir <- getwd()

setwd(home_dir)
source("scripts/utils.R")
setwd(cur_dir)

## One of the CIBERSORT training datasets
dataset <- "GSE22886"

gse.ids <- list(dataset)
names(gse.ids) <- gse.ids

gses <- llply(gse.ids, .parallel = TRUE, .fun = function(gse.id) getGEO(gse.id, GSEMatrix=TRUE))
gses <- unlist(gses)

## Get the GEO expression matrices, combining into one matrix
## First column is 'Gene'
expr.mats <- llply(gses, .fun = function(gse) exprs(gse) %>% as.data.frame)
expr.mat <- Reduce("cbind", expr.mats)
expr.mat <- drop.duplicate.columns(expr.mat) %>% rownames_to_column(var = "Gene")

probe.to.symbol.map <- get.probe.to.symbol.map(gses)

## Translate probes to symbols and to ensembl IDs
## Translate probes to symbols and to ensembl IDs
## compression.fun <- "choose.max.mad.row"
compression.fun <- "colMeans"
## compression.fun <- "choose.max.row"
symbol.compression.fun <- compression.fun
ensg.compression.fun <- compression.fun
expr.mat.symbol <- translate.genes(expr.mat, probe.to.symbol.map, fun = symbol.compression.fun)

