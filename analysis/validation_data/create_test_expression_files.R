library(magrittr)

synapser::synLogin()
devtools::source_url("https://raw.githubusercontent.com/Sage-Bionetworks/synapse_tidy_utils/master/utils.R")

dataset_file_id <- "syn21590362"
expr_file_id    <- "syn21576632"
coarse_gs_id    <- "syn21590364"
fine_gs_id      <- "syn21590365"
upload_dir_id   <- "syn21576356"

dataset_tbl   <- synapse_file_to_tbl(dataset_file_id)
coarse_gs_tbl <- synapse_file_to_tbl(coarse_gs_id, delim = ",")
fine_gs_tbl   <- synapse_file_to_tbl(fine_gs_id, delim = ",")
expr_tbl      <- synapse_file_to_tbl(expr_file_id, delim = ",")

dataset_list <- dataset_tbl %>% 
    dplyr::select(dataset, id) %>% 
    dplyr::group_by(dataset) %>% 
    dplyr::summarize(samples = list(id)) %>% 
    dplyr::ungroup() %>% 
    dplyr::arrange(dataset) %>% 
    tibble::deframe()
        

expr_tbl %>% 
    dplyr::select(c("Gene", dataset_list$DS1)) %>% 
    readr::write_csv("symbol_tpm_DS1.csv")

"symbol_tpm_DS1.csv" %>% 
    synapser::File(parent = upload_dir_id) %>% 
    synapser::synStore()

expr_tbl %>% 
    dplyr::select(c("Gene", dataset_list$DS2)) %>% 
    readr::write_csv("symbol_tpm_DS2.csv")

"symbol_tpm_DS2.csv" %>% 
    synapser::File(parent = upload_dir_id) %>% 
    synapser::synStore()
    
expr_tbl %>% 
    dplyr::select(c("Gene", dataset_list$DS3)) %>% 
    readr::write_csv("symbol_tpm_DS3.csv")

"symbol_tpm_DS3.csv" %>% 
    synapser::File(parent = upload_dir_id) %>% 
    synapser::synStore()

expr_tbl %>% 
    dplyr::select(c("Gene", dataset_list$DS4)) %>% 
    readr::write_csv("symbol_tpm_DS4.csv")

"symbol_tpm_DS4.csv" %>% 
    synapser::File(parent = upload_dir_id) %>% 
    synapser::synStore()

