library(tidyverse)
library(synapser)
library(data.table)
library(magrittr)
library(GEOquery)
library(AnnotationDbi)
library(hugene11sttranscriptcluster.db)

upload_id     <- "syn17091857"
gt_upload_id  <- "syn17091856"

source("../../scripts/utils.R")
synLogin()


gse_object <- getGEO("GSE77344", GSEMatrix = TRUE)

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
  
activity_obj <- Activity(
    name = "create",
    description = "process GEO data into usable tables",
    used = list(),
    executed = list("https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/GSE77344/create_processed_tables.R")
)

write_tsv(expr_df, "expression_affy_log.tsv")
write_tsv(linear_expr_df, "expression_affy_linear.tsv")
write_tsv(ground_truth_df, "ground_truth.tsv")
write_tsv(annotation_df, "annotation.tsv")

upload_file_to_synapse("expression_affy_log.tsv", upload_id, activity_obj = activity_obj)
upload_file_to_synapse("expression_affy_linear.tsv", upload_id, activity_obj = activity_obj)
upload_file_to_synapse("annotation.tsv", upload_id, activity_obj = activity_obj)
upload_file_to_synapse("ground_truth.tsv", gt_upload_id, activity_obj = activity_obj)

