library(tidyverse)
library(synapser)
library(synapserutils)
library(data.table)
library(magrittr)
library(GEOquery)
library(AnnotationDbi)
library(hgu133plus2.db)

upload_id     <- "syn17091411"
gt_upload_id  <- "syn17091410"
dataset       <- "GSE11057"
script_url    <- "https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/GSE11057/create_processed_tables.R"
activity_name <- "process_files_from_GEO"

source("../../scripts/utils.R")
synLogin()


gse_object <- getGEO(dataset, GSEMatrix = TRUE)

geo_df <- gse_object %>% 
    extract2(1) %>% 
    phenoData() %>% 
    pData() %>%
    rownames_to_column("sample") %>%
    as_tibble %>% 
    dplyr::select(sample, title, `fractional abundance:ch1`) %>% 
    set_colnames(c("sample", "title", "fraction")) %>% 
    separate(col = title, into = c("cell_type", "donor"), sep = ", donor ")

translation_df <- geo_df %>% 
    filter(cell_type == "PBMCs") %>% 
    dplyr::select("sample", "donor")

ground_truth_df <- geo_df %>% 
    dplyr::select(-sample) %>% 
    spread(key = "cell_type", value = "fraction") %>% 
    drop_na() %>% 
    inner_join(translation_df, .) %>% 
    dplyr::select(- c(donor, PBMCs))

query_df <-
    AnnotationDbi::select(
        hgu133plus2.db,
        keys=keys(hgu133plus2.db,keytype="PROBEID"),
        columns=c("SYMBOL"),
        keytype="PROBEID") %>%
    as_tibble() %>%
    set_colnames(c("Probe", "Hugo")) %>%
    drop_na()

expr_df <- gse_object$GSE11057_series_matrix.txt.gz@assayData %>% 
    assayDataElement('exprs') %>%
    matrix_to_df("Probe") %>%
    inner_join(query_df) %>%
    dplyr::select(Hugo, everything()) %>%
    dplyr::select(-Probe) %>%
    group_by(Hugo) %>%
    summarise_all(mean) %>%
    drop_na() %>%
    filter(Hugo != "") %>% 
    gather(key = "sample", value = "expr", -Hugo) %>% 
    mutate(sample = str_remove_all(sample, "X"))

samples_in_common <- intersect(expr_df$sample, ground_truth_df$sample)

linear_expr_df <- expr_df %>% 
    filter(sample %in% samples_in_common) %>% 
    spread(key = "sample", value = "expr")

log_expr_df <- expr_df %>% 
    filter(sample %in% samples_in_common) %>% 
    mutate(expr = log2(expr +1)) %>% 
    spread(key = "sample", value = "expr")

ground_truth_df <- ground_truth_df %>% 
    filter(sample %in% samples_in_common)



manifest_df1 <- tibble(
    path = c("expression_log.tsv",
             "expression_linear.tsv"),
    parent = upload_id,
    executed = script_url,
    activityName = activity_name,
    dataset = dataset,
    file_type = "expression",
    expression_type = "microarray",
    microarray_type = "Affymetrix HG-U133 Plus 2.0"
    expression_space = c("log2", "linear")
)

manifest_df3 <- tibble(
    path = "ground_truth.tsv",
    parent = gt_upload_id,
    executed = script_url,
    activityName = activity_name,
    dataset = dataset,
    file_type = "ground truth",
    unit = "fraction", 
    cell_types = str_c(colnames(ground_truth_df)[-1], collapse = ";")
)


write_tsv(manifest_df1, "manifest1.tsv")
write_tsv(manifest_df3, "manifest3.tsv")

write_tsv(log_expr_df, "expression_log.tsv")
write_tsv(linear_expr_df, "expression_linear.tsv")
write_tsv(ground_truth_df, "ground_truth.tsv")

syncToSynapse("manifest1.tsv")
syncToSynapse("manifest3.tsv")

