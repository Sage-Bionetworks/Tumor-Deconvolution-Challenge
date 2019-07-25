## Update the following:
source("dataset-setup.R")

suppressPackageStartupMessages(p_load(tidyverse))
suppressPackageStartupMessages(p_load(synapser))
suppressPackageStartupMessages(p_load(synapserutils))
suppressPackageStartupMessages(p_load(data.table))
suppressPackageStartupMessages(p_load(magrittr))
suppressPackageStartupMessages(p_load(Biobase))
suppressPackageStartupMessages(p_load(xlsx))
suppressPackageStartupMessages(p_load(reshape2))
suppressPackageStartupMessages(p_load(GEOquery))

## Begin configuration

file <- "create_challenge_files.R"

## End configuration

script_url <- paste0(url.base, "/", dataset, "/", file)
activity_name <- "create Challenge expression and ground truth files"

source("../../scripts/utils.R")
synLogin()

output.folder.synId <- create.folder(staging.folder.synId, dataset)

if(!grepl(dataset, pattern="GSE")) {
  stop("Was expecting an GSE dataset\n")
}

gse.ids <- list(dataset)
names(gse.ids) <- gse.ids

gses <- llply(gse.ids, .parallel = TRUE, .fun = function(gse.id) getGEO(gse.id, GSEMatrix=TRUE))
gses <- unlist(gses)

## Process ground truth

## Create a data frame holding the ground truth with columns sample, cell.type, and measured
## NB: this dataset includes cells from the HCT116 cell line--which is a CRC cell line
gt.df <- gses %>%
    get.geo.metadata.tbl %>%
    rownames_to_column(var = "sample") %>%
    dplyr::rename("hct116" = "hct116 mrna mass (ng):ch1") %>%
    dplyr::rename("monocytes" = "monocytes mrna mass (ng):ch1") %>%    
    dplyr::rename("neutrophils" = "neutrophils mrna mass (ng):ch1") %>%
    dplyr::rename("nk.cells" = "nk cells mrna mass (ng):ch1") %>%
    dplyr::rename("t.cells" = "t cells mrna mass (ng):ch1") %>%
    dplyr::rename("b.cells" = "b cells mrna mass (ng):ch1") %>%
    dplyr::select(c("sample", "monocytes", "neutrophils", "nk.cells", "t.cells", "b.cells")) %>%
    melt(id.vars = "sample") %>%
    set_colnames(c("sample", "cell.type", "measured")) %>%
    mutate(measured = as.numeric(measured))

## End processing ground truth

## Process GEO expression

## Get the GEO expression matrices, combining into one matrix
## First column is 'Gene'
expr.mats <- llply(gses, .fun = function(gse) exprs(gse) %>% as.data.frame)
if(length(expr.mats) != 1) { stop("Multiple GSEs: may be a batch effect\n") }
expr.mat <- Reduce("cbind", expr.mats)

expr.mat <- drop.duplicate.columns(expr.mat) %>% rownames_to_column(var = "Gene")

## May need to update the following get.geo.platform.name function
platform <- get.geo.platform.name(gses)
cancer.type <- "CRC"
data.processing <- unlist(get.geo.data.processing(gses))
normalization <- "RMA"
print(data.processing)

## scale <- get.log.or.linear.space(data.processing)
## Update the following based on data.processing:
scale <- "Log2"
## Update the following:
native.probe.type <- "Probe"

## Extract mappings from probe to gene symbol and Ensembl ID.
probe.to.symbol.map <- get.probe.to.symbol.map(gses)
probe.to.ensg.map <- get.probe.to.ensg.map(gses)

## End processing GEO expression

## Define the coarse- and fine-grained population mapping

fine.grained.definitions <-
  list(
##       "memory.B.cells" = c("Memory_B_cells"),
##       "naive.B.cells" = c("Naive_B_cells"),
##       "memory.CD4.T.cells" = c("Central memory", "Effector memory"),
##       "naive.CD4.T.cells" = c("Naive")
##       "regulatory.T.cells" = c("Tregs"),
##       "memory.CD8.T.cells" = c("memory.CD8.T.cells"), 
##       "naive.CD8.T.cells" = c("naive.CD8.T.cells"),
       "NK.cells" = c("nk.cells"),
       "neutrophils" = c("neutrophils"),
       "monocytes" = c("monocytes")
##       "myeloid.dendritic.cells" = c("MyeloidDC"),
##       "macrophages" = c("Macrophages.IHC")
##       "fibroblasts" = c("Fibroblasts"),
##       "endothelial.cells" = c("Endothelial")
      )

## Update the following:
coarse.grained.definitions <-
  list(
       "B.cells" = c("b.cells"),
##       "CD4.T.cells" = c("Central memory", "Effector memory", "Naive")
##       "CD8.T.cells" = c("CD8.IHC"),
       "NK.cells" = c("nk.cells"),
       "neutrophils" = c("neutrophils"),
       "monocytic.lineage" = c("monocytes")
##       "fibroblasts" = c("Fibroblasts"),
##       "endothelial.cells" = c("Endothelial")
      )
      
## End defining the coarse- and fine-grained population mapping

## May need to update the following get.geo.platform.name function
set.seed(1234)
obfuscated.dataset <- paste0("DS", sum(utf8ToInt(dataset)))

## Obfuscate the sample names
obfuscate.sample.names <- TRUE

ret <- subset.and.rename.samples(expr.mat, gt.df, obfuscate.sample.names)
expr.mat <- ret[["expr"]]
gt.df <- ret[["gt"]]
samples.map <- ret[["map"]]

## Spread ground truth into a matrix
gt.mat.raw <- gt.df %>%    
    spread(key = "cell.type", value = "measured") %>%
    dplyr::rename(sample.id = sample)

## Rename and combine ground truth columns into those required for
## coarse- and fine-grained challenges.
gt.df.coarse <- NA
if(length(coarse.grained.definitions) > 0) {
  gt.df.coarse <- map.and.format.populations(gt.mat.raw, coarse.grained.definitions)
}

gt.df.fine <- NA
if(length(fine.grained.definitions) > 0) {
  gt.df.fine <- map.and.format.populations(gt.mat.raw, fine.grained.definitions)
}

## Translate probes to symbols and to ensembl IDs
## compression.fun <- "choose.max.mad.row"
compression.fun <- "colMeans"
## compression.fun <- "choose.max.row"
symbol.compression.fun <- compression.fun
ensg.compression.fun <- compression.fun
expr.mat.symbol <- translate.genes(expr.mat, probe.to.symbol.map, fun = symbol.compression.fun)
expr.mat.ensg <- translate.genes(expr.mat, probe.to.ensg.map, fun = ensg.compression.fun)

expr.mats <- list("native" = expr.mat, "ensg" = expr.mat.ensg, "hugo" = expr.mat.symbol)
gt.mats <- list("fine" = gt.df.fine, "coarse" = gt.df.coarse)
mapping.mats <- list("symbol" = probe.to.symbol.map, "ensg" = probe.to.ensg.map)
ns <- list("n.coarse.pops" = length(unique(gt.df.coarse$cell.type)),
           "coarse.pops" = paste(sort(as.character(unique(gt.df.coarse$cell.type))), collapse=", "),
           "n.coarse" = nrow(gt.df.coarse),
	   "n.fine.pops" = length(unique(gt.df.fine$cell.type)),
           "fine.pops" = paste(sort(as.character(unique(gt.df.fine$cell.type))), collapse=", "),	   
           "n.fine" = nrow(gt.df.coarse),
           "n.samples" = ncol(expr.mat.symbol))

metadata <-
  list("dataset.name" = obfuscated.dataset,
       "orig.dataset.name" = dataset,
       "cancer.type" = cancer.type,
       "platform" = platform,
       "scale" = scale,
       "native.probe.type" = native.probe.type,
       "data.processing" = data.processing,
       "normalization" = normalization,
       "symbol.compression.function" = symbol.compression.fun,
       "ensg.compression.function" = ensg.compression.fun)

metadata <- c(metadata, ns)

identifier <- dataset
if(obfuscate.sample.names) {
  identifier <- obfuscated.dataset
}
## NB: the name of the metadata file uses (and _must_ use) the original dataset name
metadata.file.name <- paste0(dataset, "-metadata.tsv")
upload.data.and.metadata.to.synapse(identifier, expr.mats, gt.mats, mapping.mats, metadata,
                                    output.folder.synId, metadata.file.name,
                                    executed = script_url, used = NULL, sample.mapping = samples.map)

