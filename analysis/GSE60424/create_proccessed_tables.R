library(tidyverse)
library(synapser)
library(data.table)
library(magrittr)

hugo_id          <- "syn11536071"
counts_id        <- "syn17091968"
series_matrix_id <- "syn17091969"
upload_id        <- "syn17091967"
gt_upload_id     <- "syn17091966"

source("../../scripts/utils.R")
synLogin()

hugo_df <-  create_df_from_synapse_id(hugo_id)

counts_df <- counts_id %>% 
    create_df_from_synapse_id(unzip = T) %>% 
    left_join(hugo_df, by = c("genenames" = "ensembl_gene_id")) %>% 
    select(hgnc_symbol, everything()) %>% 
    select(-genenames) %>% 
    drop_na() %>% 
    group_by(hgnc_symbol) %>% 
    summarise_all(sum) %>% 
    filter(!hgnc_symbol == "") %>% 
    ungroup %>% 
    dplyr::rename(Hugo = hgnc_symbol) %>% 
    gather(key = "sample", value = "counts", - Hugo)

series_df <- series_matrix_id %>% 
    create_df_from_synapse_id(unzip = T, skip = 38, nrow = 50) %>% 
    .[c(9, 10, 11, 13, 15, 17, 19), -1] %>% 
    mutate(column = c(
        "age",
        "cell_count",
        "cell_type",
        "disease_status",
        "gender",
        "race",
        "smoker")) %>% 
    transpose_df("column", "sample") %>% 
    gather(key = "key", value = "value", -sample) %>% 
    mutate(value = str_match(value, "^[:print:]+: ([:print:]+$)")[,2]) %>% 
    mutate(value = ifelse(value == "--", NA, value)) %>% 
    spread(key = "key", value = "value") %>% 
    filter(!disease_status == "ALS")

celltype_df <- series_df %>% 
    filter(!cell_type == "Whole Blood") %>% 
    dplyr::rename(cell_type_sample = sample)

ground_truth_df <- series_df %>% 
    filter(cell_type == "Whole Blood") %>% 
    select(-c(cell_type, cell_count)) %>% 
    inner_join(celltype_df) %>% 
    select(sample, cell_type, cell_count) %>% 
    complete(sample, cell_type, fill = list(cell_count = 0)) %>% 
    spread(key = "cell_type", value = "cell_count")

annotation_df <- series_df %>% 
    select(-c(cell_type, cell_count))


samples_in_common <- 
    reduce(list(counts_df$sample, ground_truth_df$sample, annotation_df$sample), intersect)

annotation_df <- annotation_df %>% 
    filter(sample %in% samples_in_common)

ground_truth_df <- ground_truth_df %>%
    filter(sample %in% samples_in_common) %>% 
    gather(key = "cell_type", value = "count", -sample) %>% 
    mutate(count = as.numeric(count)) %>% 
    split(.$sample) %>% 
    map(mutate, total = sum(count)) %>% 
    bind_rows %>% 
    mutate(fraction = count / total) %>% 
    select(-c(count, total)) %>% 
    spread(key = "cell_type", value = "fraction")
    


counts_df <- counts_df %>%
    filter(sample %in% samples_in_common) %>%
    spread(key = "sample", value = "counts") 

cpm_m <- counts_df %>% 
    df_to_matrix("Hugo") %>% 
    calculate_cpm()

cpm_df <- cpm_m %>% 
    matrix_to_df("Hugo")
    
log_cpm_df <- cpm_m %>% 
    add(1) %>% 
    log10 %>% 
    matrix_to_df("Hugo")


activity_obj <- Activity(
    name = "create",
    description = "process GEO data into usable tables",
    used = list(hugo_id, counts_id, series_matrix_id),
    executed = list("https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/GSE60424/create_processed_tables.R")
)

write_tsv(cpm_df, "cpm.tsv")
write_tsv(log_cpm_df, "log_cpm.tsv")
write_tsv(ground_truth_df, "ground_truth.tsv")
write_tsv(annotation_df, "annotation.tsv")

upload_file_to_synapse("cpm.tsv", upload_id, activity_obj = activity_obj)
upload_file_to_synapse("log_cpm.tsv", upload_id, activity_obj = activity_obj)
upload_file_to_synapse("annotation.tsv", upload_id, activity_obj = activity_obj)
upload_file_to_synapse("ground_truth.tsv", gt_upload_id, activity_obj = activity_obj)

