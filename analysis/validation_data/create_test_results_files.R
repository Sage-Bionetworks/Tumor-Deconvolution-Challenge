library(magrittr)

synapser::synLogin()
devtools::source_url("https://raw.githubusercontent.com/Sage-Bionetworks/synapse_tidy_utils/master/utils.R")

prediction_files_tbl <- "syn21576641" %>% 
    create_entity_tbl() %>% 
    dplyr::select(name, id) %>% 
    dplyr::mutate(name = stringr::str_remove_all(name, ".csv")) %>% 
    tidyr::separate(name, into = c("method", "subchallenge"), sep = "_") %>% 
    dplyr::mutate(tbl = purrr::map(id, synapse_file_to_tbl, delim = ",")) %>% 
    dplyr::select(-id) %>% 
    tidyr::unnest(cols = c(tbl))

coarse_tbl <- "syn21590364" %>% 
    synapse_file_to_tbl(delim = ",") %>% 
    dplyr::full_join(
        dplyr::filter(prediction_files_tbl, subchallenge == "coarse"),
        .
    )

fine_tbl <- "syn21590365" %>% 
    synapse_file_to_tbl(delim = ",") %>% 
    dplyr::full_join(
        dplyr::filter(prediction_files_tbl, subchallenge == "fine"),
        .
    )

tbl <- dplyr::bind_rows(coarse_tbl, fine_tbl)

# expr_tbl %>% 
#     dplyr::select(c("Gene", dataset_list$DS4)) %>% 
#     readr::write_csv("symbol_tpm_DS4.csv")
# 
# "symbol_tpm_DS4.csv" %>% 
#     synapser::File(parent = upload_dir_id) %>% 
#     synapser::synStore()

