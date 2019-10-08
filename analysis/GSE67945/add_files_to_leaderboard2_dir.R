library(magrittr)
synapser::synLogin()
devtools::source_url("https://raw.githubusercontent.com/Sage-Bionetworks/synapse_tidy_utils/master/utils.R")

synapser::synLogin()

upload_file_to_synapse <- function(path, destination_id){
    path %>% 
        synapser::File(destination_id) %>% 
        synapser::synStore()
    
    file.remove(path)
}

leaderboard_dir_id <- "syn20555541"
gs_dir_id          <- "syn19987461"
# copy files ----
file_ids_to_copy <- c(
    "syn20745158",
    "syn20745172",
    "syn20745160",
    "syn20745179"
)

purrr::walk(
    file_ids_to_copy, 
    synapserutils::copy, 
    leaderboard_dir_id, 
    updateExisting = T
)

# input.csv ----
input_file_id <- "syn20564947"
ds500_file_id <- "syn20745114"

ds500_df <- synapse_file_to_tbl(ds500_file_id, delim = ",") 

input_file_id %>% 
    synapse_file_to_tbl(delim = ",") %>% 
    dplyr::filter(!dataset.name == "DS500") %>% 
    dplyr::bind_rows(ds500_df) %>% 
    readr::write_csv("input.csv")

upload_file_to_synapse("input.csv", leaderboard_dir_id)

# coarse gt ----
coarse_gt_id    <- "syn20564964"
ds500_coarse_id <- "syn20745131"

ds500_df <- synapse_file_to_tbl(ds500_coarse_id, delim = ",") 

coarse_gt_id  %>% 
    synapse_file_to_tbl(delim = ",") %>% 
    dplyr::filter(!dataset.name == "DS500") %>% 
    dplyr::bind_rows(ds500_df) %>% 
    readr::write_csv("lb_coarse_r2.csv")

upload_file_to_synapse("lb_coarse_r2.csv", gs_dir_id)

# fine gt ----
fine_gt_id    <- "syn20564966"
ds500_fine_id <- "syn20745128"

ds500_df <- synapse_file_to_tbl(ds500_fine_id, delim = ",") 

x <- fine_gt_id  %>% 
    synapse_file_to_tbl(delim = ",") %>% 
    dplyr::filter(!dataset.name == "DS500") %>% 
    dplyr::filter(!dataset.name == "DS483") %>% 
    dplyr::bind_rows(ds500_df) %>% 
    readr::write_csv("lb_fine_r2.csv")

upload_file_to_synapse("lb_fine_r2.csv", gs_dir_id)





    
