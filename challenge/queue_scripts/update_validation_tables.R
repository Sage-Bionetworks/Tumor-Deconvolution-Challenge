synapser::synLogin()

devtools::source_url(
    "https://raw.githubusercontent.com/Sage-Bionetworks/synapse_tidy_utils/master/utils.R"
)

query_file_dir <- "challenge/queue_scripts/queue_queries/"

tables  <- list(
    "coarse1" = "syn21071329",
    "fine1"   = "syn20782304",
    "coarse2" = "syn21071330",
    "fine2"   = "syn21046130",
    "coarse3" = "syn21425960",
    "fine3"   = "syn21411988"
    # "coarsef" = "syn21822681",
    # "finef"   = "syn21822682",
)

files <- list(
    "coarse1" = paste0(query_file_dir, "coarse1.csv"),
    "fine1"   = paste0(query_file_dir, "fine1.csv"),
    "coarse2" = paste0(query_file_dir, "coarse2.csv"),
    "fine2"   = paste0(query_file_dir, "fine2.csv"),
    "coarse3" = paste0(query_file_dir, "coarse3.csv"),
    "fine3"   = paste0(query_file_dir, "fine3.csv"),
    "coarsef" = paste0(query_file_dir, "coarse_final.csv"),
    "finef"   = paste0(query_file_dir, "fine_final.csv")
)

get_submitter_name <- function(id){
    name <- try(synapser::synGetUserProfile(id)$userName, silent = T)
    if (class(name) == "try-error") {
        name <- synapser::synGetTeam(id)$name
    }
    return(name)
}

create_table_from_query_result <- function(query_result){
    query_result %>%
        readr::read_csv(.) %>% 
        dplyr::filter(status == "ACCEPTED") %>% 
        dplyr::select(
            objectId, 
            createdOn, 
            submitterId, 
            repositoryName, 
            dplyr::contains("rounded")
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

upload_new_results_to_synapse <- function(result_file, table_id){
    synapse_table <- table_id %>%
        paste0("select objectId from ", .) %>%
        query_synapse_table()
    results_tbl <- result_file %>%
        create_table_from_query_result()
    new_results_tbl <- results_tbl %>%
        dplyr::filter(!objectId %in% synapse_table$objectId)
    old_result_ids <- synapse_table %>%
        dplyr::filter(!objectId %in% results_tbl$objectId) %>%
        dplyr::filter(objectId > 0) %>%
        dplyr::pull(objectId)
    if (nrow(new_results_tbl) > 0) {
        new_results_tbl %>%
            synapser::Table(table_id, .) %>%
            synapser::synStore()
    }
    if (length(old_result_ids) > 0) {
        old_result_ids %>%
            paste0(collapse = ", ") %>%
            paste0("select * from ", table_id, " where objectId in (", ., ")") %>%
            synapser::synTableQuery() %>%
            synapser::synDelete()
    }
}

# purrr::walk2(
#     files,
#     tables, 
#     upload_new_results_to_synapse
# )






