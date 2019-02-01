library(tidyverse)
library(synapser)
library(synapserutils)
library(data.table)
library(magrittr)
library(GEOquery)
library(AnnotationDbi)
library(hgu133plus2.db)
library(preprocessCore)

upload_id    <- "syn15664986"
gt_upload_id <- "syn15664985"

expr_gse65133_id <- "syn15667753"

source("../../scripts/utils.R")
synLogin()


dataset       <- "GSE65135"
script_url    <- "https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/GSE65135/create_processed_tables.R"
activity_name <- "Process files from Geo."

gse <- getGEO(dataset, GSEMatrix = TRUE)

m_expr_gse65133_v <- expr_gse65133_id %>% 
    create_df_from_synapse_id %>% 
    df_to_matrix("Hugo") %>% 
    rowMeans() 

series_df <- 
    pData(phenoData(gse[[1]])) %>% 
    rownames_to_column("sample") %>% 
    as_data_frame %>% 
    dplyr::select(sample, title, `flow cytometry cell subset proportions:ch1`) %>% 
    set_colnames(c("sample", "id", "cell_types")) %>% 
    filter(cell_types != "NA") %>% 
    separate(cell_types, sep = "; ", into = as.character(1:20), fill = "right") %>%
    gather(key = "key", value = "value", -c(sample, id)) %>% 
    dplyr::select(-key) %>% 
    drop_na %>% 
    separate(value, sep = " = ", into = c("cell_type", "percent")) %>% 
    mutate(cell_type = str_replace_all(cell_type, " ", "_")) %>% 
    mutate(cell_type = str_replace_all(cell_type, "ï", "i")) %>% 
    mutate(percent = str_remove_all(percent, "%")) %>% 
    spread(key = "cell_type", value = "percent")


anno_df <-         dplyr::select(series_df, sample, id)
ground_truth_df <- dplyr::select(series_df, -id)


query_df <- 
    AnnotationDbi::select(
        hgu133plus2.db, 
        keys=keys(hgu133plus2.db, keytype="PROBEID"),
        columns=c("SYMBOL"), 
        keytype="PROBEID") %>% 
    as_data_frame() %>% 
    set_colnames(c("Affy", "Hugo")) %>% 
    drop_na()

expr_df <- 
    assayDataElement(gse$GSE65135_series_matrix.txt.gz@assayData, 'exprs') %>% 
    matrix_to_df("Affy") %>% 
    dplyr::select(c("Affy", anno_df$sample)) %>% 
    inner_join(query_df) %>% 
    dplyr::select(Hugo, everything()) %>% 
    dplyr::select(-Affy) %>% 
    group_by(Hugo) %>% 
    summarise_all(max) 

log_expr_df <- expr_df %>% 
    gather(key = "sample", value = "expr", -Hugo) %>% 
    mutate(expr = log2(expr)) %>% 
    spread(key = "sample", value = "expr")

norm_log_expr_df <- expr_df %>% 
    df_to_matrix("Hugo") %>% 
    normalize.quantiles.use.target(m_expr_gse65133_v) %>% 
    set_rownames(expr_df$Hugo) %>% 
    matrix_to_df("Hugo") %>% 
    set_colnames(colnames(expr_df))


write_tsv(log_expr_df, "expression_log.tsv")
write_tsv(expr_df, "expression_linear.tsv")
write_tsv(anno_df, "annotation.tsv")
write_tsv(ground_truth_df, "ground_truth.tsv")

expression_manifest_df <- tibble(
    path = c("expression_log.tsv", "expression_linear.tsv"),
    parent = upload_id,
    executed = script_url,
    activityName = activity_name,
    dataset = dataset,
    file_type = "expression",
    expression_type = "microarray", 
    microarray_type = "Affymetrix HG-U133 Plus 2.0",
    expression_space = c("log2", "linear")
)

annotation_manifest_df <- tibble(
    path = "annotation.tsv",
    parent = upload_id,
    executed = script_url,
    activityName = activity_name,
    dataset = dataset,
    file_type = "annotations",
    annotations = str_c(colnames(anno_df)[-1], collapse = ";")
)

ground_truth_manifest_df <- tibble(
    path = "ground_truth.tsv",
    parent = gt_upload_id,
    executed = script_url,
    activityName = activity_name,
    dataset = dataset,
    file_type = "ground truth",
    unit = "percent",
    cell_types = str_c(colnames(ground_truth_df)[-1], collapse = ";")
)

write_tsv(expression_manifest_df, "expression_manifest.tsv")
write_tsv(annotation_manifest_df, "annotation_manifest.tsv")
write_tsv(ground_truth_manifest_df, "ground_truth_manifest.tsv")

syncToSynapse("expression_manifest.tsv")
syncToSynapse("annotation_manifest.tsv")
syncToSynapse("ground_truth_manifest.tsv")


write_tsv(norm_log_expr_df, "normalized_log_expression_affy.tsv")

activity_obj2 <- Activity(
    name = "create",
    description = "process GEO data into usable tables, normalize with GSE65133",
    used = list(expr_gse65133_id),
    executed = list("https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/GSE65135/create_processed_tables.R")
)

upload_file_to_synapse("normalized_log_expression_affy.tsv", upload_id, activity_obj = activity_obj2)

