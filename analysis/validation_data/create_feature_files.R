library(magrittr)

synapser::synLogin()
devtools::source_url("https://raw.githubusercontent.com/Sage-Bionetworks/synapse_tidy_utils/master/utils.R")

destination_dir_id    <- "syn20825179"
destination_file_name <- "validation_features.csv"
input_file_id         <- "syn21821122"
input_dir_id          <- "syn21821096"

file_to_synapse <- function(file, id){
    file %>% 
        synapser::File(parent = id) %>% 
        synapser::synStore()
    file.remove(file)
}

insilico_expr_files <- input_file_id %>% 
    synapse_file_to_tbl(delim = ",") %>% 
    dplyr::select("dataset.name", dplyr::contains("expr")) %>% 
    tidyr::pivot_longer(., -"dataset.name", names_to = "x", values_to = "name") %>% 
    dplyr::select("dataset" = "dataset.name", "name") %>% 
    dplyr::distinct()

input_dir_id  %>% 
    create_entity_tbl(.) %>% 
    dplyr::select("name", "id") %>% 
    dplyr::left_join(insilico_expr_files) %>% 
    dplyr::mutate("type" = dplyr::if_else(
        stringr::str_detect(.data$name, "ensg"),
        "ensg",
        "hugo"
    )) %>% 
    dplyr::mutate(
        "tbl" = purrr::map(.data$id, synapse_file_to_tbl, delim = ",")
    ) %>% 
    tidyr::unnest(cols = c(.data$tbl)) %>% 
    dplyr::select("name","dataset", "type", "Gene") %>% 
    dplyr::distinct() %>% 
    readr::write_csv(destination_file_name)

file_to_synapse(destination_file_name, destination_dir_id)
