library(tidyverse)
library(synapser)
library(synapserutils)
library(data.table)
library(magrittr)
library(GEOquery)
library(AnnotationDbi)
library(hugene11sttranscriptcluster.db)

ground_truth_id <- "syn17089680"

upload_id     <- "syn17026189"
gt_upload_id  <- "syn17026188"

source("../../scripts/utils.R")
synLogin()

dataset       <- "GSE20300"
github.path   <- "https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/"
dataset.path  <- paste0(github.path, dataset)
script_url    <- paste0(dataset.path, "/create_processed_tables.R")
used          <- NA
activity_name <- "process GEO data into usable tables"

library(xlsx)

ground_truth_df <- ground_truth_id %>% 
    download_from_synapse() %>% 
    xlsx::read.xlsx(1)

gse_object <- getGEO(dataset, GSEMatrix = TRUE)

series_df <- gse_object %>% 
    extract2(1) %>% 
    phenoData() %>% 
    pData() %>%
    rownames_to_column("sample") %>%
    as_data_frame %>% 
    arrange(`transplant state:ch1`)

translation_df <- 
    bind_cols(ground_truth_df, series_df) %>% 
    dplyr::select(sample, title, Sample.ID)

ground_truth_df <- ground_truth_df %>% 
    left_join(translation_df) %>% 
    dplyr::select(-c(Sample.ID, Patient.Group, Total, title)) %>% 
    dplyr::select(sample, everything())

annotation_df <- series_df %>% 
    dplyr::select(sample, `transplant state:ch1`) %>% 
    dplyr::rename(transplant_state = `transplant state:ch1`)

expr <-
    gse_object %>%
    extract2(1) %>%
    assayDataElement('exprs') %>%
    as.data.frame()

gpl <- getGEO(gse_object[[1]]@annotation, destdir=".")
mapping <- Table(gpl)[, c("ID", "Gene Symbol")]
colnames(mapping) <- c("from", "to")
if(!all(rownames(expr) %in% mapping$from)) {
    cat("Some probes not in mapping\n")
    table(rownames(expr) %in% mapping$from)
    stop("Stopping")
} else {
    cat("All probes in mapping\n")
}

microarray_type <- NA
if(Meta(gpl)$title == "[HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array") {
  microarray_type <- "Affymetrix HG-U133 Plus 2.0"
} else if(Meta(gpl)$title == "[HuGene-1_0-st] Affymetrix Human Gene 1.0 ST Array [transcript (gene) version]") {
  microarray_type <- "Affymetrix Human Gene 1.0 ST Array"
} else {
  stop(paste0("Unknown array type ", Meta(gpl)$title))
}

## Confirm that these data are gcRMA-normalized data--i.e., in log2 space
data.processing <-
    gse_object %>%
    extract2(1) %>%
    phenoData() %>%
    pData() %>%
    dplyr::pull(data_processing) %>%
    unique() %>%
    as.character()

cat(paste0(dataset, " data processing: ", data.processing))
if(data.processing %in% c("gcRMA", "RMA")) {
  cat(paste0("Based on ", data.processing, " processing, assuming data are in log2 space\n"))
} else {
  stop(paste0("Data were not processed via RMA, but with ", data.processing, ".  Can not confirm in log2 space\n"))
}

expr_df <- expr %>%
    compressGenes(mapping, fun = max) %>%
    matrix_to_df("Hugo") %>%
    gather(key = "sample", value = "expr", - Hugo) %>%
    filter(Hugo != "")

samples_in_common <- 
    purrr::reduce(list(expr_df$sample, ground_truth_df$sample, annotation_df$sample, translation_df$sample),
           intersect)

expr <- expr[, colnames(expr) %in% samples_in_common]

log_expr_df <- expr_df %>% 
    filter(sample %in% samples_in_common) %>% 
    spread(key = "sample", value = "expr")

linear_expr_df <- log_expr_df %>%
    df_to_matrix("Hugo") %>%
    raise_to_power(x=2, power=.) %>%
    matrix_to_df("Hugo")

ground_truth_df <- ground_truth_df %>%
    filter(sample %in% samples_in_common)

annotation_df <- annotation_df %>%
    filter(sample %in% samples_in_common)

unity <- 100
eps <- unity * 0.01

flag <- abs(rowSums(ground_truth_df[, !(colnames(ground_truth_df) == "sample")]) - unity) > eps
if(any(flag)) {
    cat("The following rows in the ground truth do not sum to one:\n")
    print(ground_truth_df[flag,,drop=F])
}

ground_truth_df <- ground_truth_df %>%
    df_to_matrix("sample") %>%
    divide_by(rowSums(.)) %>%
    matrix_to_df("sample")

expression_manifest_df <- tibble(
    path = c("expression_log.tsv", "expression_linear.tsv", "expression_log_probe.tsv"),
    parent = upload_id,
    executed = script_url,
    activityName = activity_name,
    dataset = dataset,
    used = used,
    file_type = "expression",
    expression_type = "microarray",
    microarray_type = microarray_type,
    probe_space = c("gene", "gene", "probe"),
    expression_space = c("log2", "linear", "log2")
)

annotation_manifest_df <- tibble(
    path = "annotation.tsv",
    parent = upload_id,
    executed = script_url,
    activityName = activity_name,
    dataset = dataset,
    used = used,
    file_type = "annotations",
    annotations = str_c(colnames(annotation_df)[-1], collapse = ";")
)

ground_truth_manifest_df <- tibble(
    path = "ground_truth.tsv",
    parent = gt_upload_id,
    executed = script_url,
    activityName = activity_name,
    dataset = dataset,
    used = used,
    file_type = "ground truth",
    unit = "fraction",
    cell_types = str_c(colnames(ground_truth_df)[-1], collapse = ";")
)

write_tsv(expr, "expression_log_probe.tsv")
write_tsv(log_expr_df, "expression_log.tsv")
write_tsv(linear_expr_df, "expression_linear.tsv")
write_tsv(annotation_df, "annotation.tsv")
write_tsv(ground_truth_df, "ground_truth.tsv")

write_tsv(expression_manifest_df, "expression_manifest.tsv")
write_tsv(annotation_manifest_df, "annotation_manifest.tsv")
write_tsv(ground_truth_manifest_df, "ground_truth_manifest.tsv")

syncToSynapse("expression_manifest.tsv")
syncToSynapse("annotation_manifest.tsv")
syncToSynapse("ground_truth_manifest.tsv")

