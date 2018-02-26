
require(dplyr)
require(tidyr)
require(purrr)

create_annotation_df <- function(synapse_id = "syn11898281"){
    df <- synapse_id %>% 
        create_df_from_synapse_id("./") %>% 
        dplyr::select(title, characteristics_ch1) %>% 
        tidyr::separate(characteristics_ch1, c("string", "person"), ": ") %>% 
        dplyr::mutate(person = str_replace_all(person, ";", "")) %>% 
        dplyr::select(title, person) %>% 
        purrr::set_names(c("sample", "person"))
}
