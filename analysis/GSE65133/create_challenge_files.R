## Update the following:
source("dataset-setup.R")

suppressPackageStartupMessages(p_load(tidyverse))
suppressPackageStartupMessages(p_load(synapser))
suppressPackageStartupMessages(p_load(synapserutils))
suppressPackageStartupMessages(p_load(data.table))
suppressPackageStartupMessages(p_load(magrittr))
suppressPackageStartupMessages(p_load(Biobase))
suppressPackageStartupMessages(p_load(ImmuneSpaceR))
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
metadata.file.name <- paste0(dataset, "-metadata.tsv")

if(!grepl(dataset, pattern="GSE")) {
  stop("Was expecting an GSE dataset\n")
}

gse.ids <- list(dataset)
names(gse.ids) <- gse.ids

gses <- llply(gse.ids, .parallel = TRUE, .fun = function(gse.id) getGEO(gse.id, GSEMatrix=TRUE))
gses <- unlist(gses)

## Process ground truth

## Create a data frame holding the ground truth with columns sample, cell.type, and measured
gt.df <- gses %>%
  get.geo.metadata.tbl %>%
  rownames_to_column("sample") %>% 
  as.data.frame() %>% 
  dplyr::select(sample, title, `flow cytometry cell subset proportions:ch1`) %>% 
  set_colnames(c("sample", "id", "cell_types")) %>% 
  filter(cell_types != "NA") %>% 
  separate(cell_types, sep = "; ", into = as.character(1:20), fill = "right") %>%
  gather(key = "key", value = "value", -c(sample, id)) %>% 
  dplyr::select(-key) %>%
  dplyr::select(-id) %>%   
  drop_na %>% 
  separate(value, sep = " = ", into = c("cell.type", "measured")) %>% 
  mutate(cell.type = str_replace_all(cell.type, " ", "_")) %>% 
  mutate(cell.type = str_replace_all(cell.type, "ï", "i")) %>% 
  mutate(measured = str_remove_all(measured, "%")) %>%
  mutate(measured = as.numeric(measured))

## End processing ground truth

## Process GEO expression

## Get the GEO expression matrices, combining into one matrix
## First column is 'Gene'
expr.mats <- llply(gses, .fun = function(gse) exprs(gse) %>% as.data.frame)
if(length(expr.mats) != 1) { stop("Multiple GSEs: may be a batch effect\n") }
expr.mat <- Reduce("cbind", expr.mats)

suppressPackageStartupMessages(p_load(illuminaHumanv4.db))

query_df <- 
    AnnotationDbi::select(
        illuminaHumanv4.db, 
        keys=keys(illuminaHumanv4.db,keytype="PROBEID"),
        columns=c("SYMBOL"), 
        keytype="PROBEID") %>% 
    as_tibble() %>% 
    set_colnames(c("Illum", "Hugo")) %>% 
    drop_na()


log_expr_df <- 
    expr.mat %>%
    matrix_to_df("Illum") %>% 
    inner_join(query_df) %>% 
    dplyr::select(Hugo, everything()) %>% 
    dplyr::select(-Illum) %>% 
    group_by(Hugo) %>% 
    summarise_all(max) %>% 
    drop_na() %>% 
    filter(Hugo != "") %>%
    dplyr::rename(Gene = Hugo)

expr.mat <- drop.duplicate.columns(expr.mat) %>% rownames_to_column(var = "Gene")

## May need to update the following get.geo.platform.name function
platform <- get.geo.platform.name(gses)
cancer.type <- NA
data.processing <- unlist(get.geo.data.processing(gses))
print(data.processing)
normalization <- "normexp"

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
       "memory.B.cells" = c("Memory_B_cells"),
       "naive.B.cells" = c("Naive_B_cells"),
       "memory.CD4.T.cells" = c("Resting_memory_CD4_T_cells", "Activated_memory_CD4_T_cells"),
       "naive.CD4.T.cells" = c("Naive_CD4_T_cells"),
##       "regulatory.T.cells" = c("Tregs"),
##       "memory.CD8.T.cells" = c("memory.CD8.T.cells"), 
##       "naive.CD8.T.cells" = c("naive.CD8.T.cells"),
       "NK.cells" = c("NK_cells"),
##       "neutrophils" = c("Neutrophils"),
       "monocytes" = c("Monocytes")
##       "myeloid.dendritic.cells" = c("MyeloidDC"),
##       "macrophages" = c("X"),
##       "fibroblasts" = c("Fibroblasts"),
##       "endothelial.cells" = c("Endothelial")
      )

## Update the following:
coarse.grained.definitions <-
  list("B.cells" = c("Naive_B_cells", "Memory_B_cells"),
       "CD4.T.cells" = c("Naive_CD4_T_cells", "Resting_memory_CD4_T_cells", "Activated_memory_CD4_T_cells"),
       "CD8.T.cells" = c("CD8_T_cells"),
       "NK.cells" = c("NK_cells")
##       "neutrophils" = c("Neutrophils"),
##       "monocytic.lineage" = c("Monocytes"),
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
## Translate probes to symbols and to ensembl IDs
## compression.fun <- "choose.max.mad.row"
compression.fun <- "colMeans"
## compression.fun <- "choose.max.row"
symbol.compression.fun <- compression.fun
ensg.compression.fun <- compression.fun
expr.mat.symbol <- translate.genes(expr.mat, probe.to.symbol.map, fun = symbol.compression.fun)
expr.mat.ensg <- translate.genes(expr.mat, probe.to.ensg.map, fun = ensg.compression.fun)
## expr.mat.symbol <- log_expr_df 

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
