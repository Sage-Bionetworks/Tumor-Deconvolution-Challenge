library(magrittr)

synapser::synLogin()
devtools::source_url("https://raw.githubusercontent.com/Sage-Bionetworks/synapse_tidy_utils/master/utils.R")

parameter_list <- list(
    "index_file" = list(
        "path" = "Homo_sapiens.GRCh38.cdna.all.idx",
        "class" = "File"),
    "synapse_config" = list(
        "path" = ".synapseConfig",
        "class" = "File"),
    "destination_id"= "syn22364870"
)

uploaded_samples <- "syn22364870" %>% 
    create_entity_tbl() %>%
    dplyr::pull("name") %>% 
    stringr::str_remove_all(".tsv")

"syn22364869" %>% 
    create_entity_tbl() %>% 
    dplyr::select("id", "name") %>% 
    dplyr::mutate("name" = stringr::str_remove_all(name, ".fastq.gz")) %>% 
    tidyr::separate("name", sep = "_", into = c("sample", "pair")) %>% 
    tidyr::pivot_wider(names_from = "pair", values_from = "id") %>% 
    dplyr::filter(!sample %in% uploaded_samples) %>% 
    dplyr::select("fastq1_ids" = "R1", "fastq2_ids" = "R2", "sample_name" = "sample") %>% 
    as.list() %>% 
    c(parameter_list) %>%
    RJSONIO::toJSON() %>%
    writeLines("kallisto.json")
