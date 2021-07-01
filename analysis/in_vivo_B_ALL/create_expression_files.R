library(magrittr)

synapser::synLogin()
devtools::source_url(
    "https://raw.githubusercontent.com/Sage-Bionetworks/synapse_tidy_utils/master/utils.R"
)

kallisto_tbl <- "syn22364870" %>% 
    create_entity_tbl(.) %>% 
    dplyr::mutate("sample" = stringr::str_remove_all(.data$name, ".tsv")) %>% 
    dplyr::select("id", "sample") %>% 
    dplyr::mutate(.tbl = purrr::map(id, synapse_file_to_tbl)) 

tpm_tbl <- kallisto_tbl %>% 
    dplyr::mutate(.tbl = purrr::map(
        .tbl, 
        ~dplyr::select(.x, transcript = target_id, tpm
    ))) %>% 
    tidyr::unnest(cols = .tbl) %>% 
    dplyr::select(-id) %>% 
    tidyr::pivot_wider(names_from = sample, values_from = tpm)

counts_tbl <- kallisto_tbl %>% 
    dplyr::mutate(.tbl = purrr::map(.tbl, ~dplyr::select(
        .x, 
        transcript = target_id, est_counts
    ))) %>% 
    tidyr::unnest(cols = .tbl) %>% 
    dplyr::select(-id) %>% 
    tidyr::pivot_wider(names_from = sample, values_from = est_counts)

readr::write_tsv(tpm_tbl, "transcript_tpms.tsv")  
readr::write_tsv(counts_tbl, "transcript_est_counts.tsv")  
synapser::File("transcript_tpms.tsv", parent = "syn22492020") %>% 
    synapser::synStore()
synapser::File("transcript_est_counts.tsv", parent = "syn22492020") %>% 
    synapser::synStore()
file.remove("transcript_tpms.tsv")
file.remove("transcript_est_counts.tsv")
