synapser::synLogin()

devtools::source_url(
    "https://raw.githubusercontent.com/Sage-Bionetworks/synapse_tidy_utils/master/utils.R"
)

source("challenge/queue_scripts/submission_functions.R")


fine_images <-  query_synapse_table("SELECT id, status FROM syn23517766")
dplyr::count(fine_images, status)

coarse_images <- query_synapse_table("SELECT id, status FROM syn23517749")
dplyr::count(coarse_images, status)


fine_images %>%
    dplyr::filter(status == "INVALID") %>%
    dplyr::pull(id) %>%
    purrr::walk(reset_submission)

coarse_images %>%
    dplyr::filter(status == "INVALID") %>%
    dplyr::pull(id) %>%
    purrr::walk(reset_submission)

