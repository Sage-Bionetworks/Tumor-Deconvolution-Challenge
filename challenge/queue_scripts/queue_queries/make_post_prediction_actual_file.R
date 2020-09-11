library(tidyverse)
library(synapser)
library(magrittr)

synapser::synLogin()

devtools::source_url("https://raw.githubusercontent.com/Sage-Bionetworks/synapse_tidy_utils/master/utils.R")

coarse_submissions <- "select id, dockerrepositoryname from syn22274598" %>% 
    query_synapse_table() %>% 
    dplyr::select(id, repo_name = dockerrepositoryname)

fine_submissions <- "select id, dockerrepositoryname from syn22276620" %>% 
    query_synapse_table() %>% 
    dplyr::select(id, repo_name = dockerrepositoryname)

prediction_file_df <- "syn19518404" %>%
    create_entity_tbl() %>%
    dplyr::select(id, name) %>%
    dplyr::mutate(objectId = stringr::str_remove_all(name, ".csv$")) %>%
    dplyr::select(-name)

fine_gold_standard <- synapse_file_to_tbl("syn21820376",  delim = ",")

fine_predictions <- prediction_file_df %>% 
    dplyr::right_join(fine_submissions, by = c("objectId" = "id")) %>% 
    dplyr::mutate(data = purrr::map(id, synapse_file_to_tbl, delim = ",")) %>%
    dplyr::select(-id) %>% 
    tidyr::unnest(cols = c(data)) %>% 
    dplyr::full_join(fine_gold_standard, by = c("dataset.name", "sample.id", "cell.type")) %>% 
    dplyr::mutate(subchallenge = "fine")

coarse_gold_standard <- synapse_file_to_tbl("syn22267267",  delim = ",")

coarse_predictions <- prediction_file_df %>% 
    dplyr::right_join(coarse_submissions, by = c("objectId" = "id")) %>% 
    dplyr::mutate(data = purrr::map(id, synapse_file_to_tbl, delim = ",")) %>%
    dplyr::select(-id) %>% 
    tidyr::unnest(cols = c(data)) %>% 
    dplyr::full_join(coarse_gold_standard, by = c("dataset.name", "sample.id", "cell.type")) %>% 
    dplyr::mutate(subchallenge = "coarse")

readr::write_csv(
    dplyr::bind_rows(fine_predictions, coarse_predictions),
    "post_predictions.csv"
)

"post_predictions.csv" %>%
    synapser::File("syn16925027") %>%
    synapser::synStore()
