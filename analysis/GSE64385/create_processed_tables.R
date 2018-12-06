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

## ground truth: IHC_MCP.TXT (sent by Aurelien)
gt_id   <- "syn17014397"

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

stop("stop")

anno_df <- annotations %>% 
    as_data_frame %>%
    dplyr::rename("age" = "age.at.diagnosis (year):ch1") %>%
    dplyr::rename("gender" = "Sex:ch1") %>%
    select(sample, age, gender)

activity_obj <- Activity(
    name = "download-annotations",
    description = "download GEO annotations",
    used = NULL,
    executed = list("https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/GSE39582/create_processed_tables.R")
)

write_tsv(anno_df, "annotation.tsv")
upload_file_to_synapse("annotation.tsv", preprocessed_upload_id, activity_obj = activity_obj)

expr <- as.data.frame(exprs(gse[[1]]))

activity_obj <- Activity(
    name = "download-expression",
    description = "download raw GEO files",
    used = NULL,
    executed = list("https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/GSE39582/create_processed_tables.R")
)

write_tsv(expr, "GSE39582-expr-probes.tsv")
upload_file_to_synapse("GSE39582-expr-probes.tsv", raw_upload_id, activity_obj = activity_obj)

## Translate the probe-based expression to gene-based expression
gpl <- getGEO(gse[[1]]@annotation, destdir=".")
mapping <- Table(gpl)[, c("ID", "Gene Symbol")]
colnames(mapping) <- c("from", "to")
if(!all(rownames(ex) %in% mapping$from)) {
    cat("Some probes not in mapping\n")
    table(rownames(ex) %in% mapping$from)
    stop("Stopping")
} else {
    cat("All probes in mapping\n")
}

library(plyr)
## Translate/compress genes from one name space (e.g., probe ids) to another (e.g., symbols)
## Take the max probe for each gene as that gene's expression
compressGenes <- function(e, mapping, from.col = "from", to.col = "to")
{
  e$to    <- mapping[match(rownames(e), mapping[, from.col]), to.col]
  e           <- e[!is.na(e$to),]
  e           <- ddply(.data = e, .variables = "to", .fun = function(x){apply(x[,-ncol(x)],2,max)},.parallel = T)
  rownames(e) <- e$to
  e           <- e[,-1]
  return(e)
}

expr_symbols <- expr %>% compressGenes(mapping) %>% matrix_to_df("Hugo")

write_tsv(expr_symbols, "expression_microarray.tsv")

## Get the Synapse ID of the raw file we saved above
children <- synGetChildren(raw_upload_id)
l <- as.list(children)
df <- do.call(rbind.data.frame, l)

raw_expr_id <- as.character(df$id[df$name == "GSE39582-expr-probes.tsv"])

activity_obj <- Activity(
    name = "create",
    description = "process GEO expression files",
    used = list(raw_expr_id),
    executed = list("https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/GSE39582/create_processed_tables.R")
)
upload_file_to_synapse("expression_microarray.tsv", preprocessed_upload_id, activity_obj = activity_obj)

## Process the ground truth file
gt_df <- gt_id %>% 
    create_df_from_synapse_id %>%
    dplyr::rename(sample = CEL.ID)

write_tsv(gt_df, "ground_truth.tsv")

activity_obj <- Activity(
    name = "create",
    description = "standardize format of raw ground truth file provided by Aurelien",
    used = list(gt_id),
    executed = list("https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/GSE39582/create_processed_tables.R")
)

upload_file_to_synapse("ground_truth.tsv", preprocessed_upload_id, activity_obj = activity_obj)

