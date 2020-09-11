delete_submission <- function(id){
    synapser::synDelete(synapser::synGetSubmission(id))
}

reset_submission <- function(id){
    status <- synapser::synGetSubmissionStatus(id)
    status$status <- "RECEIVED"
    synapser::synStore(status)
}

make_submission <- function(entity_id, evaluation_id){
    synapser::synSubmit(
        synapser::synGetEvaluation(evaluation_id),
        synapser::synGet(entity_id)
    )
}