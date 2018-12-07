library(synapser)
library(data.table)
library(magrittr)
library(tidyverse)

home_dir <- "../../../Tumor-Deconvolution-Challenge/"
tmp_dir <- tempdir()

setwd(home_dir)
source("scripts/utils.R")
setwd(tmp_dir)
synLogin()

## Folder: leaderboard_datasets/GSE64385/
dataset_upload_id <- "syn17088595"
## Folder: leaderboard_datasets/GSE64385/pre-processed
preprocessed_upload_id <- "syn17088596"
## Folder: leaderboard_datasets/GSE64385/raw
raw_upload_id <- "syn17088597"

## Begin download/processing from GEO

suppressPackageStartupMessages(library("foreach"))
suppressPackageStartupMessages(library("parallel"))

num.cores <- detectCores()
if(!is.na(num.cores) && (num.cores > 1)) {
  suppressPackageStartupMessages(library("doMC"))
  cat(paste("Registering ", num.cores-1, " cores.\n", sep=""))
  registerDoMC(cores=(num.cores-1))
}

library(GEOquery)
gse <- getGEO("GSE64385", GSEMatrix=TRUE)

annotations <- pData(phenoData(gse[[1]]))
annotations$sample <- rownames(annotations)

anno_df <- annotations %>% 
    as_data_frame %>%
    dplyr::rename("hct116.mrna.mass" = "hct116 mrna mass (ng):ch1") %>%
    dplyr::rename("monocytes.mrna.mass" = "monocytes mrna mass (ng):ch1") %>%    
    dplyr::rename("neutrophils.mrna.mass" = "neutrophils mrna mass (ng):ch1") %>%
    dplyr::rename("nk.cells.mrna.mass" = "nk cells mrna mass (ng):ch1") %>%
    dplyr::rename("t.cells.mrna.mass" = "t cells mrna mass (ng):ch1") %>%
    dplyr::rename("b.cells.mrna.mass" = "b cells mrna mass (ng):ch1") %>%    
    dplyr::mutate(hct116.mrna.mass = as.numeric(hct116.mrna.mass)) %>%
    dplyr::mutate(monocytes.mrna.mass = as.numeric(monocytes.mrna.mass)) %>%
    dplyr::mutate(neutrophils.mrna.mass = as.numeric(neutrophils.mrna.mass)) %>%
    dplyr::mutate(nk.cells.mrna.mass = as.numeric(nk.cells.mrna.mass)) %>%
    dplyr::mutate(t.cells.mrna.mass = as.numeric(t.cells.mrna.mass)) %>%
    dplyr::mutate(b.cells.mrna.mass = as.numeric(b.cells.mrna.mass)) %>%    
    dplyr::mutate(hct116.mrna.percent = hct116.mrna.mass /
                  (hct116.mrna.mass + monocytes.mrna.mass + neutrophils.mrna.mass +
		   nk.cells.mrna.mass + t.cells.mrna.mass + b.cells.mrna.mass) ) %>%
    dplyr::mutate(monocytes.mrna.percent = monocytes.mrna.mass /
                  (hct116.mrna.mass + monocytes.mrna.mass + neutrophils.mrna.mass +
		   nk.cells.mrna.mass + t.cells.mrna.mass + b.cells.mrna.mass) ) %>%
    dplyr::mutate(neutrophils.mrna.percent = neutrophils.mrna.mass /
                  (hct116.mrna.mass + monocytes.mrna.mass + neutrophils.mrna.mass +
		   nk.cells.mrna.mass + t.cells.mrna.mass + b.cells.mrna.mass) ) %>%
    dplyr::mutate(nk.cells.mrna.percent = nk.cells.mrna.mass /
                  (hct116.mrna.mass + monocytes.mrna.mass + neutrophils.mrna.mass +
		   nk.cells.mrna.mass + t.cells.mrna.mass + b.cells.mrna.mass) ) %>%
    dplyr::mutate(t.cells.mrna.percent = t.cells.mrna.mass /
                  (hct116.mrna.mass + monocytes.mrna.mass + neutrophils.mrna.mass +
		   nk.cells.mrna.mass + t.cells.mrna.mass + b.cells.mrna.mass) ) %>%
    dplyr::mutate(b.cells.mrna.percent = b.cells.mrna.mass /
                  (hct116.mrna.mass + monocytes.mrna.mass + neutrophils.mrna.mass +
		   nk.cells.mrna.mass + t.cells.mrna.mass + b.cells.mrna.mass) ) %>%
    select(sample,
           "hct116.mrna.mass", "monocytes.mrna.mass", "neutrophils.mrna.mass",
           "nk.cells.mrna.mass", "t.cells.mrna.mass", "b.cells.mrna.mass",
           "hct116.mrna.percent", "monocytes.mrna.percent", "neutrophils.mrna.percent",
           "nk.cells.mrna.percent", "t.cells.mrna.percent", "b.cells.mrna.percent")

activity_obj <- Activity(
    name = "download-annotations",
    description = "download GEO annotations",
    used = list("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64385"),
    executed = list("https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/GSE64385/create_processed_tables.R")
)

write_tsv(anno_df, "annotation.tsv")
upload_file_to_synapse("annotation.tsv", preprocessed_upload_id, activity_obj = activity_obj)

expr <- as.data.frame(exprs(gse[[1]]))

activity_obj <- Activity(
    name = "download-expression",
    description = "download raw GEO files",
    used = list("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64385"),
    executed = list("https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/GSE64385/create_processed_tables.R")
)

write_tsv(expr, "GSE64385-expr-probes.tsv")
upload_file_to_synapse("GSE64385-expr-probes.tsv", raw_upload_id, activity_obj = activity_obj)

## Translate the probe-based expression to gene-based expression
gpl <- getGEO(gse[[1]]@annotation, destdir=".")
mapping <- Table(gpl)[, c("ID", "Gene Symbol")]
colnames(mapping) <- c("from", "to")
if(!all(rownames(expr) %in% mapping$from)) {
    cat("Some probes not in mapping\n")
    table(rownames(expr) %in% mapping$from)
    stop("Stopping")
} else {
    cat("All probes in mapping\n")
}

compressed <- expr %>% compressGenes(mapping, fun = max) 
expr_symbols <- compressed %>% matrix_to_df("Hugo")

## NB: according to the annotations, these data are output by RMA, which is in log space.
write_tsv(expr_symbols, "log_expression_microarray.tsv")

## Get the Synapse ID of the raw file we saved above
children <- synGetChildren(raw_upload_id)
l <- as.list(children)
df <- do.call(rbind.data.frame, l)

raw_expr_id <- as.character(df$id[df$name == "GSE64385-expr-probes.tsv"])

activity_obj <- Activity(
    name = "create",
    description = "process GEO expression files",
    used = list(raw_expr_id),
    executed = list("https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/GSE64385/create_processed_tables.R")
)
upload_file_to_synapse("log_expression_microarray.tsv", preprocessed_upload_id, activity_obj = activity_obj)

## Process the ground truth file
gt_df <- anno_df

write_tsv(gt_df, "ground_truth.tsv")

activity_obj <- Activity(
    name = "create",
    description = "standardize format of raw ground truth file",
    used = list("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64385"),
    executed = list("https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/GSE64385/create_processed_tables.R")
)

upload_file_to_synapse("ground_truth.tsv", preprocessed_upload_id, activity_obj = activity_obj)

