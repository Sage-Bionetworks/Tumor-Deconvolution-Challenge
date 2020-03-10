library(tidyverse)
library(synapser)

devtools::source_url(
    "https://raw.githubusercontent.com/Sage-Bionetworks/synapse_tidy_utils/master/utils.R"
)

synapser::synLogin()

"SELECT id FROM syn21570575 WHERE type = 'Purified Cell'" %>% 
    query_synapse_table(.) %>% 
    dplyr::pull(id) %>% 
    paste("synapse get", ., "--downloadLocation fastqs" ) %>% 
    writeLines("synapse_download_cell.sh")

"SELECT id FROM syn21570575 WHERE type = 'Biological Mix'" %>% 
    query_synapse_table(.) %>% 
    dplyr::pull(id) %>% 
    paste("synapse get", ., "--downloadLocation bms" ) %>% 
    writeLines("synapse_download_bm.sh")

"SELECT id FROM syn21570575 WHERE type = 'Random Mix'" %>% 
    query_synapse_table(.) %>% 
    dplyr::pull(id) %>% 
    paste("synapse get", ., "--downloadLocation rms" ) %>% 
    writeLines("synapse_download_rm.sh")
