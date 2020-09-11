synapser::synLogin()

devtools::source_url(
    "https://raw.githubusercontent.com/Sage-Bionetworks/synapse_tidy_utils/master/utils.R"
)

source("submission_functions.R")

tbl <- "syn16925027" %>% 
    create_entity_tbl() %>% 
    dplyr::filter(type == "org.sagebionetworks.repo.model.docker.DockerRepository") %>% 
    dplyr::select(id) %>% 
    dplyr::mutate(name = purrr::map_chr(
        id, 
        ~synapser::synGet(.x, downloadFile = F)$properties$repositoryName   
    ))

tbl %>% 
    dplyr::filter(stringr::str_detect(name, "co[au]rse")) %>% 
    dplyr::pull(id) %>% 
    purrr::walk(make_submission, 9614582L)

tbl %>% 
    dplyr::filter(stringr::str_detect(name, "fine")) %>% 
    dplyr::pull(id) %>% 
    purrr::walk(make_submission, 9614583L)
