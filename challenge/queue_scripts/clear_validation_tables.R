synapser::synLogin()

devtools::source_url(
    "https://raw.githubusercontent.com/Sage-Bionetworks/synapse_tidy_utils/master/utils.R"
)

tables  <- list(
    # coarse_tbl_id = "syn21822681",
    fine_tbl_id   = "syn21822682"
)

clear_synapse_tbl <- function(id){
    current <- synapser::synTableQuery(paste0("SELECT * FROM ", id))
    synapser::synDelete(current) 
}

purrr::walk(tables, clear_synapse_tbl)
