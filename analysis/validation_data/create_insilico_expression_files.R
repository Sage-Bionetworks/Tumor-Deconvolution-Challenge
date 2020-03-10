library(magrittr)

synapser::synLogin()
devtools::source_url(
    "https://raw.githubusercontent.com/Sage-Bionetworks/synapse_tidy_utils/master/utils.R"
)

kallisto_tbl <- "SELECT id, sample FROM syn21571153" %>% 
    query_synapse_table() %>% 
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
synapser::File("transcript_tpms.tsv", parent = "syn21624952") %>% 
    synapser::synStore()
synapser::File("transcript_est_counts.tsv", parent = "syn21624952") %>% 
    synapser::synStore()
file.remove("transcript_tpms.tsv")
file.remove("transcript_est_counts.tsv")
