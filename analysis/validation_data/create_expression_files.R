library(magrittr)

synapser::synLogin()
devtools::source_url("https://raw.githubusercontent.com/Sage-Bionetworks/synapse_tidy_utils/master/utils.R")

kallisto_tbl <- "SELECT id, sample FROM syn21571153" %>% 
    query_synapse_table()

expression_tbl <- kallisto_tbl %>% 
    dplyr::mutate(.tbl = purrr::map(id, synapse_file_to_tbl)) %>% 
    dplyr::mutate(.tbl = purrr::map(.tbl, ~dplyr::select(.x, transcript = target_id, tpm))) %>% 
    tidyr::unnest(cols = .tbl) %>% 
    dplyr::select(-id) %>% 
    tidyr::pivot_wider(names_from = sample, values_from = tpm)

readr::write_tsv(expression_tbl, "transcript_expression.tsv")    
synapser::File("transcript_expression.tsv", parent = "syn21571479") %>% 
    synapser::synStore()
file.remove("transcript_expression.tsv")
