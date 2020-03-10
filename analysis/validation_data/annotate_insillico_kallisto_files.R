library(magrittr)

synapser::synLogin()
devtools::source_url(
    "https://raw.githubusercontent.com/Sage-Bionetworks/synapse_tidy_utils/master/utils.R"
)
    
fastq_tbl <- query_synapse_table("SELECT sample, type FROM syn21570575")

insillico_kallisto_tbl <- query_synapse_table(
    "SELECT id, name FROM syn21609795"
)

insillico_annotations_tbl <- insillico_kallisto_tbl %>% 
    dplyr::mutate(sample = stringr::str_remove_all(name, ".tsv")) %>% 
    dplyr::select(-name) %>% 
    dplyr::inner_join(fastq_tbl) %>% 
    dplyr::rename(entity = id) %>% 
    dplyr::distinct() %>% 
    tidyr::nest(annotations = -entity)

result <- purrr::pmap(insillico_annotations_tbl, synapser::synSetAnnotations)






    
