## Update the following:
source("dataset-setup.R")

suppressPackageStartupMessages(p_load(tidyverse))
suppressPackageStartupMessages(p_load(synapser))
suppressPackageStartupMessages(p_load(synapserutils))
suppressPackageStartupMessages(p_load(data.table))
suppressPackageStartupMessages(p_load(magrittr))
suppressPackageStartupMessages(p_load(xlsx))
suppressPackageStartupMessages(p_load(reshape2))

## Begin configuration

file <- "create_challenge_files.R"

expr_id         <- "syn18144910"
ground_truth_id <- "syn18144911"

## End configuration

script_url <- paste0(url.base, "/", dataset, "/", file)
used          <- str_c(expr_id, ground_truth_id, sep = ";")
activity_name <- "create Challenge expression and ground truth files"

source("../../scripts/utils.R")
synLogin()

## Store data and metadata in this folder: staging/<dataset>
output.folder.synId <- create.folder(staging.folder.synId, dataset)

metadata.file.name <- paste0(dataset, "-metadata.tsv")

fpkm_df <- expr_id %>% 
    create_df_from_synapse_id() %>% 
    dplyr::rename(Gene = Gene.Symbol) %>% 
    gather(key = "sample", value = "expr", -"Gene")

ground_truth_df <- ground_truth_id %>% 
    download_from_synapse() %>% 
    read.xlsx(sheetIndex = 1) %>% 
    dplyr::rename(cell_type = `NA.`) %>%
    gather(key = "sample", value = "value", -"cell_type")

total_df <- ground_truth_df %>%
    group_by(sample) %>%
    dplyr::summarise(total = sum(value))

ground_truth_df <- ground_truth_df %>%
    inner_join(total_df) %>%
    mutate(value = value / total) %>%
    select(-total) %>%
    dplyr::rename(cell.type = cell_type) %>%
    dplyr::rename(measured = value)

samples_in_common <- intersect(fpkm_df$sample, ground_truth_df$sample)

## Obfuscate the sample names
obfuscate.sample.names <- TRUE
samples.map <- data.frame(from = as.character(samples_in_common), to = as.character(samples_in_common))
if(obfuscate.sample.names) {
  new_samples_in_common <- paste0("S", 1:length(samples_in_common))
  map <- data.frame(sample = samples_in_common, new.sample = new_samples_in_common)
  ground_truth_df <- merge(ground_truth_df, map) %>%
    select(-sample) %>%
    dplyr::rename(sample.id = new.sample)
  fpkm_df <- merge(fpkm_df, map) %>%
    select(-sample) %>%
    dplyr::rename(sample = new.sample)
  samples.map <- data.frame(from = as.character(samples_in_common), to = as.character(new_samples_in_common))    
}

gt.mat.raw <- ground_truth_df %>%
    spread(key = "cell.type", value = "measured")

linear_expr_df <- fpkm_df %>%
    spread(key = "sample", value = "expr") 

log_expr_df <- fpkm_df %>%
    mutate(expr = log2(expr + 1)) %>% 
    spread(key = "sample", value = "expr")

expr.mat <- linear_expr_df

fine.grained.definitions <-
  list(
##       "memory.B.cells" = c("memory.B.cells"),
##       "naive.B.cells" = c("naive.B.cells"),
##       "memory.CD4.T.cells" = c("memory.CD4.T.cells"),
##       "naive.CD4.T.cells" = c("naive.CD4.T.cells"),
##       "regulatory.T.cells" = c("Tregs"),
##       "memory.CD8.T.cells" = c("memory.CD8.T.cells"), 
##       "naive.CD8.T.cells" = c("naive.CD8.T.cells"),
       "NK.cells" = c("NK Cells"),
       "neutrophils" = c("Neutrophils"),
       "monocytes" = c("Monocytes"),
       "myeloid.dendritic.cells" = c("MyeloidDC"),
##       "macrophages" = c("X"),
       "fibroblasts" = c("Fibroblasts"),
       "endothelial.cells" = c("Endothelial")
      )

## Update the following:
coarse.grained.definitions <-
  list("B.cells" = c("B Cells"),
       "CD4.T.cells" = c("CD4"),
       "CD8.T.cells" = c("CD8"),
       "NK.cells" = c("NK Cells"),
       "neutrophils" = c("Neutrophils"),
       "monocytic.lineage" = c("Monocytes", "MyeloidDC"),
       "fibroblasts" = c("Fibroblasts"),
       "endothelial.cells" = c("Endothelial")
      )
      
gt.mat.fine <- combine.columns(gt.mat.raw, fine.grained.definitions)
rownames(gt.mat.fine) <- gt.mat.raw$sample.id
gt.mat.coarse <- combine.columns(gt.mat.raw, coarse.grained.definitions)
rownames(gt.mat.coarse) <- gt.mat.raw$sample.id

## May need to update the following get.geo.platform.name function
set.seed(1234)
obfuscated.dataset <- paste0("DS", sum(utf8ToInt(dataset)))

gt.df.fine <- melt(as.matrix(gt.mat.fine)) %>%
  dplyr::rename(sample.id = Var1) %>%
  dplyr::rename(cell.type = Var2) %>%
  dplyr::rename(measured = value) %>%
  add_column(dataset.name = obfuscated.dataset) %>%
  select(dataset.name, sample.id, cell.type, measured)

gt.df.coarse <- melt(as.matrix(gt.mat.coarse)) %>%
  dplyr::rename(sample.id = Var1) %>%
  dplyr::rename(cell.type = Var2) %>%
  dplyr::rename(measured = value) %>%
  add_column(dataset.name = obfuscated.dataset) %>%
  select(dataset.name, sample.id, cell.type, measured)

## Extract mappings from probe to gene symbol and Ensembl ID.
probe.to.symbol.map <- unique(data.frame(from = as.character(expr.mat$Gene), to = as.character(expr.mat$Gene)))
probe.to.ensg.map <- get.symbol.to.ensg.map(expr.mat$Gene)

## Translate symbols to ensg
expr.mat.symbol <- expr.mat

## Translate probes to symbols and to ensembl IDs
## compression.fun <- "choose.max.mad.row"
compression.fun <- "colMeans"
## compression.fun <- "choose.max.row"
symbol.compression.fun <- "identity"
ensg.compression.fun <- compression.fun

tmp.map <- subset(probe.to.ensg.map, from %in% expr.mat$Gene)
tmp.map <- tmp.map[!duplicated(tmp.map$to), ]
expr.mat.ensg <- expr.mat %>%
  column_to_rownames(var = "Gene") %>%
  aggregate_rows(., tmp.map, fun = ensg.compression.fun, parallel = TRUE) %>%
  rownames_to_column(var = "Gene")

expr.mats <- list("native" = expr.mat, "ensg" = expr.mat.ensg, "hugo" = expr.mat.symbol)
gt.mats <- list("fine" = gt.df.fine, "coarse" = gt.df.coarse)
mapping.mats <- list("symbol" = probe.to.symbol.map, "ensg" = probe.to.ensg.map)

cancer.type <- "CRC"
platform <- "RNA-seq"
scale <- "Linear"
native.probe.type <- "Hugo"
data.processing <- "FPKM"
normalization <- "FPKM"


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

identifier <- dataset
if(obfuscate.sample.names) {
  identifier <- obfuscated.dataset
}
## NB: the name of the metadata file uses (and _must_ use) the original dataset name
metadata.file.name <- paste0(dataset, "-metadata.tsv")
upload.data.and.metadata.to.synapse(identifier, expr.mats, gt.mats, mapping.mats, metadata,
                                    output.folder.synId, metadata.file.name,
                                    executed = script_url, used = NULL, sample.mapping = samples.map)
