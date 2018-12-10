library(tidyverse)
library(synapser)
library(data.table)
library(magrittr)
library(GEOquery)
library(AnnotationDbi)
library(hugene11sttranscriptcluster.db)

upload_id     <- "syn17091886"
gt_upload_id  <- "syn17091885"

source("../../scripts/utils.R")
synLogin()


gse_object <- getGEO("GSE40240", GSEMatrix = TRUE)

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
        "race" = "ethnicity:ch1",
        "basophils" =  "relative.basophils:ch1",                
        "eosinophils" = "relative.eosinophils:ch1",
        "lymphocytes" = "relative.lymphocytes:ch1",         
        "monocytes" = "relative.monocytes:ch1",
        "neutrophils" = "relative.neutrophils:ch1") 


query_df <-
    AnnotationDbi::select(
        hugene11sttranscriptcluster.db,
        keys=keys(hugene11sttranscriptcluster.db,keytype="PROBEID"),
        columns=c("SYMBOL"),
        keytype="PROBEID") %>%
    as_data_frame() %>%
    set_colnames(c("Probe", "Hugo")) %>%
    drop_na()

expr_df <- gse_object$GSE40240_series_matrix.txt.gz@assayData %>% 
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

ground_truth_df <- geo_df %>%
    dplyr::select(-c(age, gender, race)) %>% 
    filter(sample %in% samples_in_common)

annotation_df <- geo_df %>% 
    dplyr::select(sample, age, gender, race) %>% 
    filter(sample %in% samples_in_common)
  
activity_obj <- Activity(
    name = "create",
    description = "process GEO data into usable tables",
    used = list(),
    executed = list("https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/GSE40240/create_processed_tables.R")
)

write_tsv(expr_df, "expression_affy.tsv")
write_tsv(ground_truth_df, "ground_truth.tsv")
write_tsv(annotation_df, "annotation.tsv")

upload_file_to_synapse("expression_affy.tsv", upload_id, activity_obj = activity_obj)
upload_file_to_synapse("annotation.tsv", upload_id, activity_obj = activity_obj)
upload_file_to_synapse("ground_truth.tsv", gt_upload_id, activity_obj = activity_obj)

