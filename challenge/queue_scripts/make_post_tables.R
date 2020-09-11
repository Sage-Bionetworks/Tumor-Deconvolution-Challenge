library(tidyverse)

synapser::synLogin()

devtools::source_url(
    "https://raw.githubusercontent.com/Sage-Bionetworks/synapse_tidy_utils/master/utils.R"
)

project_id <- "syn16925027"

create_table_from_query_result <- function(query_result){
    df <- query_result %>%
        dplyr::select(
            objectId, 
            repositoryName, 
            contains("rounded")
        ) %>% 
        dplyr::rename(repo_name = repositoryName) %>% 
        tidyr::gather(
            "metric", 
            "metric_value", 
            -c(objectId, repo_name)
        ) %>% 
        tidyr::separate(
            metric, 
            sep = "_", 
            into = c("dataset", "celltype", "metric", "rounded"),
            fill = "left"
        ) %>%  
        dplyr::filter(!is.na(rounded)) %>%
        dplyr::select(-rounded) %>% 
        dplyr::mutate(
            dataset = dplyr::case_when(
                is.na(dataset) & is.na(celltype)  ~ "Grand mean",
                is.na(dataset)                    ~ "Celltype mean",
                T                                 ~ dataset
            )
        ) %>% 
        dplyr::mutate(
            celltype = dplyr::case_when(
                dataset ==  "Grand mean" ~ "Grand mean",
                T                        ~ celltype
            )
        )
}

clear_synapse_tbl <- function(id){
    current <- synapser::synTableQuery(paste0("SELECT * FROM ", id))
    synapser::synDelete(current) 
}


# tbls <- purrr::map(
#     table_tbl$query_result_file,
#     create_table_from_query_result
# )
# 
# purrr::map2(table_tbl$name, tbls, ~synapser::synBuildTable(.x, project_id, .y)) %>%
#     purrr::walk(synapser::synStore)

coarse_tbl <-  "challenge/queue_scripts/queue_queries/coarse_post.csv" %>% 
    readr::read_csv(.) %>% 
    create_table_from_query_result() %>% 
    synapser::synBuildTable("Post Challenge Coarse Results", project_id, .) %>% 
    synapser::synStore()


fine_tbl <-  "challenge/queue_scripts/queue_queries/" %>% 
    list.files(full.names = T) %>% 
    purrr::keep(., stringr::str_detect(., "fine_post")) %>% 
    purrr::map(readr::read_csv) %>% 
    dplyr::bind_rows() %>% 
    dplyr::distinct() %>% 
    create_table_from_query_result() %>% 
    synapser::synBuildTable("Post Challenge Fine Results", project_id, .) %>% 
    synapser::synStore()


coarse_post.csv %>% 
    readr::read_csv(.) %>% 
    create_table_from_query_result() %>% 
    synapser::synBuildTable("Post Challenge Coarse Results", project_id, .) %>% 
    synapser::synStore()
    















