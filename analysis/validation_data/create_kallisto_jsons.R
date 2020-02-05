library(magrittr)

synapser::synLogin()
devtools::source_url("https://raw.githubusercontent.com/Sage-Bionetworks/synapse_tidy_utils/master/utils.R")

parameter_list <- list(
    "kallisto_index_file" = list(
        "path" = "Homo_sapiens.GRCh38.cdna.all.idx",
        "class" = "File"),
    "synapse_config" = list(
        "path" = ".synapseConfig",
        "class" = "File"),
    "destination_id" = "syn21570813"
)

uploaded_samples <- "select sample from syn21571153" %>%
    query_synapse_table() %>%
    dplyr::pull(sample)

fastq_tbl <- 
    "SELECT id, pair, sample, size, type FROM syn21570575" %>%
    query_synapse_table() %>%  
    dplyr::filter(!sample %in% uploaded_samples) %>% 
    dplyr::group_by(sample) %>% 
    dplyr::mutate(size = sum(size)) %>% 
    dplyr::ungroup() %>% 
    tidyr::pivot_wider(names_from = pair, values_from = id)
   

df_to_json <- function(dfs, output_file_names, json_file_names){
    dfs %>% 
        purrr::map(as.list) %>% 
        purrr::map(c, parameter_list) %>% 
        purrr::map2(output_file_names, ~c(.x, list("output_file_name" = .y))) %>% 
        purrr::map(RJSONIO::toJSON) %>% 
        purrr::walk2(json_file_names, writeLines)
}

# Purified Cell

fastq_tbl %>%
    dplyr::filter(type == "Purified Cell") %>%
    dplyr::mutate(cum_size = cumsum(size)) %>%
    dplyr::filter(cum_size < 2.0e11) %>%
    dplyr::select(fastq1_ids = R1, fastq2_ids = R2, sample_names = sample) %>% 
    as.list() %>%
    c(parameter_list) %>%
    RJSONIO::toJSON() %>% 
    writeLines("PC.json")


# Biological Mix 

fastq_tbl %>%
    dplyr::filter(type == "Biological Mix") %>%
    dplyr::mutate(cum_size = cumsum(size)) %>%
    dplyr::select(fastq1_ids = R1, fastq2_ids = R2, sample_names = sample) %>% 
    as.list() %>%
    c(parameter_list) %>%
    RJSONIO::toJSON() %>% 
    writeLines("BM.json")


# Random Mix

fastq_tbl %>%
    dplyr::filter(type == "Random Mix") %>%
    dplyr::mutate(cum_size = cumsum(size)) %>%
    dplyr::select(fastq1_ids = R1, fastq2_ids = R2, sample_names = sample) %>% 
    as.list() %>%
    c(parameter_list) %>%
    RJSONIO::toJSON() %>% 
    writeLines("RM.json")
