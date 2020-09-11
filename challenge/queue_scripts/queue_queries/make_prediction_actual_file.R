library(tidyverse)
library(synapser)
library(magrittr)

synapser::synLogin()

devtools::source_url("https://raw.githubusercontent.com/Sage-Bionetworks/synapse_tidy_utils/master/utils.R")

query_result_files <- list(
    "coarse_final.csv",
    c("fine_final1.csv", "fine_final2.csv", "fine_final3.csv", "fine_final4.csv", "fine_final5.csv") 
)

gs_ids <- c(
    "syn21820375",
    "syn21820376"
)

subchallenge <- c(
    "coarse",
    "fine"
)

methods <- list(
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

coarse_df <- "coarse_final.csv" %>% 
    readr::read_csv(.) %>% 
    dplyr::distinct() %>% 
    dplyr::mutate(round = "final", subchallenge = "coarse") 
    
fine_df <-  
    c(
        "fine_final1.csv",
        "fine_final2.csv", 
        "fine_final3.csv", 
        "fine_final4.csv", 
        "fine_final5.csv"
    ) %>% 
    purrr::map(readr::read_csv) %>%
    dplyr::bind_rows() %>% 
    dplyr::distinct() %>% 
    dplyr::mutate(round = "final", subchallenge = "fine") 
    

submission_df <- list(coarse_df, fine_df) %>% 
    dplyr::bind_rows() %>%
    dplyr::filter(status == "ACCEPTED") %>% 
    dplyr::select(round, subchallenge, objectId, createdOn, submitterId, repositoryName) %>%
    dplyr::group_by(submitterId, round, subchallenge) %>%
    dplyr::arrange(desc(createdOn)) %>%
    dplyr::mutate(n = 1:dplyr::n()) %>%
    dplyr::mutate(is_latest = dplyr::if_else(n == 1, T, F)) %>%
    dplyr::ungroup() %>% 
    dplyr::select(-c(n, createdOn))

prediction_file_df <- "syn19518404" %>%
    create_entity_tbl() %>%
    dplyr::select(id, name) %>%
    dplyr::mutate(objectId = as.numeric(stringr::str_remove_all(name, ".csv$"))) %>%
    dplyr::select(-name)

cibersortx_coarse_df <- "coarse_cibersortx.csv" %>% 
    readr::read_csv() %>% 
    dplyr::mutate("subchallenge" = "coarse")

cibersortx_fine_df <- "fine_cibersortx.csv" %>% 
    readr::read_csv() %>% 
    dplyr::mutate("subchallenge" = "fine")

cibersort_df <- 
    dplyr::bind_rows(cibersortx_coarse_df, cibersortx_fine_df) %>% 
    dplyr::mutate(
        "round" = "final",
        "objectId" = 0L,
        "repo_name" = "baseline_method7",
        "is_latest" = F,
        "submitterId" = 3360851L
    )

prediction_df <- submission_df %>%
    dplyr::left_join(prediction_file_df) %>%
    dplyr::mutate(data = purrr::map(id, synapse_file_to_tbl, delim = ",")) %>%
    unnest(cols = c(data)) %>%
    dplyr::mutate(repo_name = as.character(ifelse(
        repositoryName %in% names(methods),
        methods[repositoryName],
        basename(repositoryName)
    ))) %>%
    dplyr::select(
        round,
        subchallenge,
        dataset.name,
        sample.id,
        submitterId,
        is_latest,
        objectId,
        repo_name,
        cell.type,
        prediction) %>% 
    dplyr::bind_rows(cibersort_df, .)




gs_df <- gs_ids %>%
    purrr::map(synapse_file_to_tbl, delim = ",") %>%
    purrr::map2(subchallenge, ~dplyr::mutate(.x, subchallenge = .y)) %>%
    dplyr::bind_rows() %>% 
    dplyr::mutate(round = "final")

df <- prediction_df %>%
    dplyr::left_join(gs_df)

readr::write_csv(df, "final_predictions.csv")

"final_predictions.csv" %>%
    synapser::File("syn16925027") %>%
    synapser::synStore()
