old_tables  <- list(
    "coarse1" = "syn21071329",
    "fine1"   = "syn20782304",
    "coarse2" = "syn21071330",
    "fine2"   = "syn21046130",
    "coarse3" = "syn21425960",
    "fine3"   = "syn21411988"
)

tables <- c("syn21822682", "syn21822681")

synapser::synLogin()

devtools::source_url(
    "https://raw.githubusercontent.com/Sage-Bionetworks/synapse_tidy_utils/master/utils.R"
)

old_submissions <- old_tables %>% 
    paste0("select submitterId, submitter from ", .) %>% 
    purrr::map(query_synapse_table) %>% 
    dplyr::bind_rows() %>% 
    dplyr::distinct()


ids <- tables %>% 
    paste0("select submitterId from ", .) %>% 
    purrr::map(query_synapse_table) %>% 
    dplyr::bind_rows() %>% 
    dplyr::distinct() %>% 
    dplyr::pull(submitterId)

send_message <- function(submitterId, submitter){
    synapser::synSendMessage(
        list(as.character(submitterId)), 
        "Tumor Deconvolution Final Round",
        paste0(
            "Dear ", submitter,
            "\n\n",
            "The tumor deconvolution DREAM challenge is set to end on Friday, ",
            "May 22nd at 1PM Pacific US time.",
            "\n",
            "We noticed that you previously made a successful submission to an ",
            "earlier round, but have not submitted to the final validation round.",
            "\n\n",
            "Please let us know if you would still like to submit and an extension ",
            "would help you do so. We would likely consider extensions of a few ",
            "days up to two weeks. Please let us know your preference.",
            "\n\n",
            "If you do not plan to submit, may be include your last successful ",
            "submission in the validation round by executing your ",
            "previously-submitted Docker container against the validation data?",
            "\n\n",
            "Thank you,",
            "\n\n",
            "Tumor Deconvolution Challenge Administrators."
        )
    )
}

old_submissions %>% 
    dplyr::filter(!submitterId %in% ids) %>%
    purrr::pmap(send_message)


    