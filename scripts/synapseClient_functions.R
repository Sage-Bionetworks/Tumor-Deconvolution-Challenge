##require(synapseClient)
require(synapser)
require(stringr)
library(dplyr)
require(data.table)
require(magrittr)

create_df_from_synapse_id <- function(syn_id, location = NULL, unzip = F, ...){
    path <- download_from_synapse(syn_id, location)
    if(unzip) path <- stringr::str_c("zcat ", path)
    path %>% 
        data.table::fread(...) %>% 
        dplyr::as_data_frame() 
}

download_from_synapse <- function(syn_id, location = NULL){
    path <- synGet(syn_id, downloadLocation = location)$path
    return(path)
}

upload_file_to_synapse <- function(
    path, synapse_id, 
    annotation_list = NULL, 
    activity_obj = NULL, 
    ret = "entity"){
    
    entity <- File(
        path = path, 
        parentId = synapse_id, 
        annotations = annotation_list)
    entity <- synStore(entity, activity = activity_obj)
    if(ret == "entity") return(entity)
    if(ret == "syn_id") return(entity$properties$id)
}

get_file_df_from_synapse_dir_id <- function(syn_id){
    str_c('select id, name from file where parentId=="', syn_id, '"') %>% 
        synQuery()
}