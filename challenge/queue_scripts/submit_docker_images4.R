synapser::synLogin()

devtools::source_url(
    "https://raw.githubusercontent.com/Sage-Bionetworks/synapse_tidy_utils/master/utils.R"
)

source("challenge/queue_scripts/submission_functions.R")

all_images <- "syn16925027" %>% 
    create_entity_tbl() %>% 
    dplyr::filter(type == "org.sagebionetworks.repo.model.docker.DockerRepository") %>% 
    dplyr::select(id) %>% 
    dplyr::mutate(name = purrr::map_chr(
        id, 
        ~synapser::synGet(.x, downloadFile = F)$properties$repositoryName   
    )) 

team_images <- all_images %>% 
    dplyr::mutate(
        "submission_id" = stringr::str_match(
            .data$name, 
            "[:print:]+/([:print:]+)$"
        )[,2],
        "submission_id" = as.integer(.data$submission_id)
    ) %>% 
    tidyr ::drop_na()

coarse_submitted_images <- 
    query_synapse_table("SELECT id FROM syn22280982 where status = 'ACCEPTED' AND createdBy <> 3360851") %>% 
    dplyr::pull(id)

fine_submitted_images <- 
    query_synapse_table("SELECT id FROM syn22287785 where status = 'ACCEPTED' AND createdBy <> 3360851") %>% 
    dplyr::pull(id)

# team_images %>%
#     dplyr::filter(.data$submission_id %in% coarse_submitted_images) %>%
#     dplyr::pull(id) %>%
#     purrr::walk(make_submission, 9614655L)
# 
# team_images %>%
#     dplyr::filter(.data$submission_id %in% fine_submitted_images) %>%
#     dplyr::pull(id) %>%
#     purrr::walk(make_submission, 9614656L)
# 
# all_images %>%
#     dplyr::filter(stringr::str_detect(name, "co[au]rse")) %>%
#     dplyr::pull(id) %>%
#     purrr::walk(make_submission, 9614655L)
# 
# all_images %>%
#     dplyr::filter(stringr::str_detect(name, "fine")) %>%
#     dplyr::pull(id) %>%
#     purrr::walk(make_submission, 9614656L)
