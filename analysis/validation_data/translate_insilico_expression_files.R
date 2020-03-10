library(magrittr)

synapser::synLogin()
devtools::source_url(
    "https://raw.githubusercontent.com/Sage-Bionetworks/synapse_tidy_utils/master/utils.R"
)

translation_tbl <- synapse_file_to_tbl("syn21574276")
counts_transcripts_tbl <- "syn21574262" %>% 
    synapse_file_to_tbl() %>% 
    dplyr::mutate(
        transcript = stringr::str_remove_all(transcript, ".[:digit:]+$")
    ) 
    
counts_ensg_tbl <-   
    dplyr::inner_join(
        translation_tbl,
        counts_transcripts_tbl,
        by = c("ensembl_transcript_id" = "transcript")
    ) %>% 
    dplyr::select(-c(ensembl_transcript_id, hgnc_symbol)) %>% 
    dplyr::group_by(ensembl_gene_id) %>% 
    dplyr::summarise_all(sum) %>% 
    dplyr::ungroup() %>% 
    dplyr::rename(Gene = ensembl_gene_id)

readr::write_csv(counts_ensg_tbl, "ensg_counts.csv")   
synapser::File("ensg_counts.csv", parent = "syn21571479") %>% 
    synapser::synStore()
file.remove("ensg_counts.csv")

counts_symbol_tbl <- counts_ensg_tbl %>% 
    dplyr::inner_join(
        translation_tbl,
        .,
        by = c("ensembl_gene_id" = "Gene")
    ) %>% 
    dplyr::select(-c(ensembl_transcript_id, ensembl_gene_id)) %>% 
    tidyr::drop_na() %>% 
    dplyr::group_by(hgnc_symbol) %>% 
    dplyr::summarise_all(mean) %>% 
    dplyr::ungroup() %>% 
    dplyr::rename(Gene = hgnc_symbol)

readr::write_csv(counts_symbol_tbl, "symbol_counts.csv")   
synapser::File("symbol_counts.csv", parent = "syn21571479") %>% 
    synapser::synStore()
file.remove("symbol_counts.csv")

tpm_transcripts_tbl <- "syn21574261" %>% 
    synapse_file_to_tbl() %>% 
    dplyr::mutate(
        transcript = stringr::str_remove_all(transcript, ".[:digit:]+$")
    ) 

tpm_ensg_tbl <-   
    dplyr::inner_join(
        translation_tbl,
        tpm_transcripts_tbl,
        by = c("ensembl_transcript_id" = "transcript")
    ) %>% 
    dplyr::select(-c(ensembl_transcript_id, hgnc_symbol)) %>% 
    dplyr::group_by(ensembl_gene_id) %>% 
    dplyr::summarise_all(sum) %>% 
    dplyr::ungroup() %>% 
    dplyr::rename(Gene = ensembl_gene_id)

readr::write_csv(tpm_ensg_tbl, "ensg_tpm.csv")   
synapser::File("ensg_tpm.csv", parent = "syn21571479") %>% 
    synapser::synStore()
file.remove("ensg_tpm.csv")

tpm_symbol_tbl <- tpm_ensg_tbl %>% 
    dplyr::inner_join(
        translation_tbl,
        .,
        by = c("ensembl_gene_id" = "Gene")
    ) %>% 
    dplyr::select(-c(ensembl_transcript_id, ensembl_gene_id)) %>% 
    tidyr::drop_na() %>% 
    dplyr::group_by(hgnc_symbol) %>% 
    dplyr::summarise_all(mean) %>% 
    dplyr::ungroup() %>% 
    dplyr::rename(Gene = hgnc_symbol)

readr::write_csv(tpm_symbol_tbl, "symbol_tpm.csv")   
synapser::File("symbol_tpm.csv", parent = "syn21571479") %>% 
    synapser::synStore()
file.remove("symbol_tpm.csv")

