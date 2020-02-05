library(magrittr)

synapser::synLogin()
devtools::source_url("https://raw.githubusercontent.com/Sage-Bionetworks/synapse_tidy_utils/master/utils.R")
    
fastq_tbl <- "SELECT sample, celltype, replicate, type FROM syn21570575" %>% 
    query_synapse_table()

kallisto_tbl <- "SELECT id, name FROM syn21571153" %>% 
    query_synapse_table()

annotations_tbl <- kallisto_tbl %>% 
    dplyr::mutate(sample = stringr::str_remove_all(name, ".tsv")) %>% 
    dplyr::select(-name) %>% 
    dplyr::inner_join(fastq_tbl) %>% 
    dplyr::rename(entity = id) %>% 
    tidyr::nest(annotations = -entity)

purrr::pmap(annotations_tbl, synapser::synSetAnnotations)






    
