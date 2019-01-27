library(tidyverse)
library(synapser)
library(synapserutils)
library(data.table)
library(magrittr)
library(GEOquery)
library(AnnotationDbi)
library(hugene11sttranscriptcluster.db)

upload_id     <- "syn17091857"
gt_upload_id  <- "syn17091856"

source("../../scripts/utils.R")
synLogin()

script_url <- "https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/GSE77344/create_processed_tables.R"
dataset <- "GSE77344"
activity_name <- "process GEO files into usable tables"


gse_object <- getGEO(dataset, GSEMatrix = TRUE)

geo_df <- gse_object %>% 
    extract2(1) %>% 
    phenoData() %>% 
    pData() %>%
    rownames_to_column("sample") %>%
    as_data_frame %>%
    dplyr::select(
        "sample", 
        "age" = "age:ch1", 
        "gender" = "Sex:ch1",
        "cd4_t_lymphocyte" = "cell prop. cd4 t lymphocyte:ch1",
        "cd8_t_lymphocyte" = "cell prop. cd8 t lymphocyte:ch1", 
        "granulocyte" = "cell prop. granulocyte:ch1",
        "monocyte" = "cell prop. monocyte:ch1",
        "nk_lymphocyte" = "cell prop. nk lymphocyte:ch1") 


query_df <-
    AnnotationDbi::select(
        hugene11sttranscriptcluster.db,
        keys=keys(hugene11sttranscriptcluster.db,keytype="PROBEID"),
        columns=c("SYMBOL"),
        keytype="PROBEID") %>%
    as_data_frame() %>%
    set_colnames(c("Probe", "Hugo")) %>%
    drop_na()

## NB: these data are in log2 space
expr_df <- gse_object$GSE77344_series_matrix.txt.gz@assayData %>% 
    assayDataElement('exprs') %>%
    matrix_to_df("Probe") %>%
    inner_join(query_df) %>%
    dplyr::select(Hugo, everything()) %>%
    dplyr::select(-Probe) %>%
    drop_na() %>%
    filter(Hugo != "") %>% 
    group_by(Hugo) %>%
    summarise_all(mean) %>%
    gather(key = "sample", value = "expr", -Hugo) %>% 
    mutate(sample = str_remove_all(sample, "X"))

samples_in_common <- intersect(expr_df$sample, geo_df$sample)

expr_df <- expr_df %>% 
    filter(sample %in% samples_in_common) %>% 
    spread(key = "sample", value = "expr")

linear_expr_df <- expr_df %>% 
    df_to_matrix("Hugo") %>% 
    raise_to_power(x=2, power=.) %>% 
    matrix_to_df("Hugo")

ground_truth_df <- geo_df %>%
    dplyr::select(-c(age, gender)) %>% 
    filter(sample %in% samples_in_common)

annotation_df <- geo_df %>% 
    dplyr::select(sample, age, gender) %>% 
    filter(sample %in% samples_in_common)

write_tsv(expr_df, "expression_log.tsv")
write_tsv(linear_expr_df, "expression_linear.tsv")
write_tsv(ground_truth_df, "ground_truth.tsv")
write_tsv(annotation_df, "annotation.tsv")

expression_manifest_df <- tibble(
    path = c("expression_log.tsv", "expression_linear.tsv"),
    parent = upload_id,
    executed = script_url,
    activityName = activity_name,
    dataset = dataset,
    file_type = "expression",
    expression_type = "microarray", 
    microarray_type = "Affymetrix Human Gene 1.1 ST Array",
    expression_space = c("log2", "linear")
)

annotation_manifest_df <- tibble(
    path = "annotation.tsv",
    parent = upload_id,
    executed = script_url,
    activityName = activity_name,
    dataset = dataset,
    file_type = "annotations",
    annotations = str_c(colnames(annotation_df)[-1], collapse = ";")
)

ground_truth_manifest_df <- tibble(
    path = "ground_truth.tsv",
    parent = gt_upload_id,
    executed = script_url,
    activityName = activity_name,
    dataset = dataset,
    file_type = "ground truth",
    unit = "fraction",
    cell_types = str_c(colnames(ground_truth_df)[-1], collapse = ";")
)


write_tsv(expression_manifest_df, "expression_manifest.tsv")
write_tsv(annotation_manifest_df, "annotation_manifest.tsv")
write_tsv(ground_truth_manifest_df, "ground_truth_manifest.tsv")

syncToSynapse("expression_manifest.tsv")
syncToSynapse("annotation_manifest.tsv")
syncToSynapse("ground_truth_manifest.tsv")

