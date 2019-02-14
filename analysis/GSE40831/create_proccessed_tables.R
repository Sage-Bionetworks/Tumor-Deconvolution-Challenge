library(tidyverse)
library(synapser)
library(synapserutils)
library(data.table)
library(magrittr)
library(GEOquery)
library(AnnotationDbi)
library(hgu133plus2.db)
library(xlsx)

upload_id     <- "syn18139812"
gt_id         <- "syn18139845"

source("../../scripts/utils.R")
synLogin()

script_url    <- "https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/GSE40831/create_processed_tables.R"
dataset       <- "GSE40831"
activity_name <- "process GEO files into usable tables"
used          <- gt_id





mononucleated_ground_truth_df <- gt_id %>% 
    download_from_synapse() %>% 
    read.xlsx(sheetIndex = 1, startRow = 2, endRow = 13, stringsAsFactors = F) %>% 
    inset("sample", value = "mononucleated cell")

lin_ground_truth_df <- gt_id %>% 
    download_from_synapse() %>% 
    read.xlsx(sheetIndex = 1, startRow = 16, endRow = 27, stringsAsFactors = F) %>% 
    inset("sample", value = "lin- cell")

ground_truth_df <-
    bind_rows(mononucleated_ground_truth_df, lin_ground_truth_df) %>% 
    as_tibble() %>% 
    dplyr::select(sample, Reference.populations, Flow.cytometry) %>% 
    set_colnames(c("sample", "cell_type", "percent")) %>% 
    mutate(percent = as.double(percent)) %>% 
    filter(!is.na(percent)) %>% 
    mutate(sample = c("lin_minus_cell", "mononucleated_cell"))

    

gse_object <- getGEO(dataset, GSEMatrix = TRUE)

annotation_df <- gse_object %>% 
    extract2(2) %>% 
    phenoData() %>% 
    pData() %>%
    rownames_to_column("sample") %>%
    as_tibble() %>% 
    dplyr::select(sample, `cell type:ch1`) %>% 
    dplyr::rename(cell_type = `cell type:ch1`)



query_df <-
    AnnotationDbi::select(
        hgu133plus2.db,
        keys=keys(hgu133plus2.db,keytype="PROBEID"),
        columns=c("SYMBOL"),
        keytype="PROBEID") %>%
    as_tibble() %>%
    set_colnames(c("AFFY", "Hugo")) %>%
    drop_na()

expr_df <- gse_object$`GSE40831-GPL571_series_matrix.txt.gz`@assayData %>%
    assayDataElement('exprs') %>%
    matrix_to_df("AFFY") %>%
    inner_join(query_df) %>%
    dplyr::select(Hugo, everything()) %>%
    dplyr::select(-AFFY) %>%
    group_by(Hugo) %>%
    summarise_all(mean) %>%
    drop_na() %>%
    filter(Hugo != "") %>%
    gather(key = "sample", value = "expr", -Hugo) %>% 
    left_join(annotation_df) %>% 
    group_by(cell_type, Hugo) %>% 
    dplyr::summarise(expr = mean(expr)) %>% 
    dplyr::rename(sample = cell_type) 

linear_expr_df <- expr_df %>%
    mutate(expr = 2^expr) %>% 
    spread(key = "sample", value = "expr") %>% 
    set_colnames(c("Hugo","lin_minus_cell", "mononucleated_cell"))

log_expr_df <- expr_df %>%
    spread(key = "sample", value = "expr") %>% 
    set_colnames(c("Hugo","lin_minus_cell", "mononucleated_cell"))

ground_truth_df <- ground_truth_df %>%
    spread(key = "cell_type", value = "percent")


expression_manifest_df <- tibble(
    path = c("expression_log.tsv", "expression_linear.tsv"),
    parent = upload_id,
    executed = script_url,
    activityName = activity_name,
    dataset = dataset,
    used = used,
    file_type = "expression",
    expression_type = "microarray",
    microarray_type = "Affymetrix HG-U133 Plus 2.0",
    expression_space = c("log2", "linear")
)

ground_truth_manifest_df <- tibble(
    path = "ground_truth.tsv",
    parent = upload_id,
    executed = script_url,
    activityName = activity_name,
    dataset = dataset,
    used = used,
    file_type = "ground truth",
    unit = "percent",
    cell_types = str_c(colnames(ground_truth_df)[-1], collapse = ";")
)

write_tsv(log_expr_df, "expression_log.tsv")
write_tsv(linear_expr_df, "expression_linear.tsv")
write_tsv(ground_truth_df, "ground_truth.tsv")

write_tsv(expression_manifest_df, "expression_manifest.tsv")
write_tsv(ground_truth_manifest_df, "ground_truth_manifest.tsv")

syncToSynapse("expression_manifest.tsv")
syncToSynapse("ground_truth_manifest.tsv")
