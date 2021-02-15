synapser::synLogin()

devtools::source_url(
    "https://raw.githubusercontent.com/Sage-Bionetworks/synapse_tidy_utils/master/utils.R"
)

source("challenge/queue_scripts/submission_functions.R")


fine_images <-  query_synapse_table("SELECT id, status FROM syn22396079")
dplyr::count(fine_images, status)

coarse_images <- query_synapse_table("SELECT id, status FROM syn22396060")
dplyr::count(coarse_images, status)


# fine_images %>%
#     dplyr::filter(status == "INVALID") %>%
#     dplyr::pull(id) %>%
#     purrr::walk(reset_submission)
# 
# coarse_images %>%
#     dplyr::filter(status == "INVALID") %>%
#     dplyr::pull(id) %>%
#     purrr::walk(reset_submission)

# fine_images %>%
#     dplyr::pull(id) %>%
#     purrr::walk(reset_submission, status_string = "INVALID")
# 
# coarse_images %>%
#     dplyr::pull(id) %>%
#     purrr::walk(reset_submission, status_string = "INVALID")
