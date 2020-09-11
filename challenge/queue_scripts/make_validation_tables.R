library(tidyverse)

synapser::synLogin()

devtools::source_url(
    "https://raw.githubusercontent.com/Sage-Bionetworks/synapse_tidy_utils/master/utils.R"
)

query_file_dir <- "challenge/queue_scripts/queue_queries/"

project_id <- "syn16925027"
# project_id <- "syn15589870"

get_submitter_name <- function(id){
    name <- try(synapser::synGetUserProfile(id)$userName, silent = T)
    if (class(name) == "try-error") {
        name <- synapser::synGetTeam(id)$name
    }
    return(name)
}

create_table_from_query_result <- function(query_result){
    df <- query_result %>%
        readr::read_csv(.) %>% 
        dplyr::select(
            objectId, 
            createdOn, 
            submitterId, 
            repositoryName, 
            contains("rounded")
        ) %>% 
        dplyr::mutate(repo_name = as.character(ifelse(
            repositoryName %in% names(baseline_methods),
            baseline_methods[repositoryName],
            basename(repositoryName)
        ))) %>% 
        dplyr::select(-repositoryName) %>% 
        dplyr::mutate(submitter = purrr::map_chr(submitterId, get_submitter_name)) %>% 
        dplyr::group_by(submitterId) %>%
        dplyr::arrange(desc(createdOn)) %>% 
        dplyr::mutate(n = 1:dplyr::n()) %>% 
        dplyr::mutate(is_latest = dplyr::if_else(n == 1, T, F)) %>% 
        dplyr::ungroup() %>% 
        dplyr::select(-c(n, createdOn)) %>% 
        tidyr::gather(
            "metric", 
            "metric_value", 
            -c(objectId, is_latest, submitter, submitterId, repo_name)
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

baseline_methods <- list(
    "docker.synapse.org/syn16925027/cibersort_course" = "baseline_method1",
    "docker.synapse.org/syn16925027/mcpcounter_course" = "baseline_method2",
    "docker.synapse.org/syn16925027/quantiseq_course" = "baseline_method3",
    "docker.synapse.org/syn16925027/xcell_course" = "baseline_method4",
    "docker.synapse.org/syn16925027/epic_course2" = "baseline_method5",
    "docker.synapse.org/syn16925027/timer_coarse" = "baseline_method6",
    "docker.synapse.org/syn16925027/cibersort_fine" = "baseline_method1",
    "docker.synapse.org/syn16925027/mcpcounter_fine" = "baseline_method2",
    "docker.synapse.org/syn16925027/quantiseq_fine" = "baseline_method3",
    "docker.synapse.org/syn16925027/xcell_fine" = "baseline_method4",
    "docker.synapse.org/syn16925027/epic_fine" = "baseline_method5",
    "docker.synapse.org/syn16925027/timer_fine" = "baseline_method6",
    "docker.synapse.org/syn20505859/epic" = "CelEst"
)

table_tbl <- dplyr::tibble(
    query_result_file = c(
        paste0(query_file_dir, "coarse_final.csv"),
        paste0(query_file_dir, "fine_final.csv")
    ),
    name = c(
        "Validation Coarse Results",
        "Validation Fine Results"
    )
)

tbls <- purrr::map(
    table_tbl$query_result_file,
    create_table_from_query_result
)

purrr::map2(table_tbl$name, tbls, ~synapser::synBuildTable(.x, project_id, .y)) %>%
    purrr::walk(synapser::synStore)















