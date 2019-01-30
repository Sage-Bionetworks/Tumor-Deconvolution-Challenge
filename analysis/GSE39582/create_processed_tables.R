library(synapser)
library(synapserutils)
library(data.table)
library(magrittr)
library(tidyverse)

home_dir <- "/home/aelamb/repos/Tumor-Deconvolution-Challenge/"
home_dir <- "../../../Tumor-Deconvolution-Challenge/"
tmp_dir  <- "/home/aelamb/tmp/tumor_deconvolution/GSE59654/"
tmp_dir <- tempdir()

setwd(home_dir)
source("scripts/utils.R")
setwd(tmp_dir)
synLogin()

## Folder: leaderboard_datasets/GSE39582/
dataset_upload_id <- "syn17014395"
## Folder: leaderboard_datasets/GSE39582/pre-processed
preprocessed_upload_id <- "syn17014396"
## Folder: leaderboard_datasets/GSE39582/raw
raw_upload_id <- "syn17014401"

## ground truth: IHC_MCP.TXT (sent by Aurelien)
gt_id   <- "syn17014397"

dataset       <- "GSE39582"
script_url    <- paste0("https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/", dataset, "/create_processed_tables.R")
activity_name <- "process_files_from_GEO"

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
gse <- getGEO("GSE39582", GSEMatrix=TRUE)

annotations <- pData(phenoData(gse[[1]]))
annotations$sample <- rownames(annotations)

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

## This is RMA-normalized data (according to annotations), which is in log2 space
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
if(!all(rownames(expr) %in% mapping$from)) {
    cat("Some probes not in mapping\n")
    table(rownames(expr) %in% mapping$from)
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

log2_expr_symbols <- expr %>% compressGenes(mapping) %>% matrix_to_df("Hugo")
linear_expr_symbols <- log2_expr_symbols %>%
    df_to_matrix("Hugo") %>%
    raise_to_power(x=2, power=.) %>%
    matrix_to_df("Hugo")

write_tsv(log2_expr_symbols, "expression_log.tsv")
write_tsv(linear_expr_symbols, "expression_linear.tsv")

## Get the Synapse ID of the raw file we saved above
children <- synGetChildren(raw_upload_id)
l <- as.list(children)
df <- do.call(rbind.data.frame, l)

raw_expr_id <- as.character(df$id[df$name == "GSE39582-expr-probes.tsv"])

microarray_type <- NA
if(Meta(gpl)$title == "[HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array") {
  microarray_type <- "Affymetrix HG-U133 Plus 2.0"
} else {
  stop(paste0("Unknown array type", Meta(gpl)$title))
}

manifest_df1 <- tibble(
    path = c("expression_log.tsv",
             "expression_linear.tsv"),
    parent = preprocessed_upload_id,
    used = raw_expr_id,
    executed = script_url,
    activityName = "process GEO expression files",
    dataset = dataset,
    file_type = "expression",
    expression_type = "microarray",
    microarray_type = microarray_type,
    expression_space = c("log2", "linear")
)

write_tsv(manifest_df1, "manifest1.tsv")

syncToSynapse("manifest1.tsv")

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

