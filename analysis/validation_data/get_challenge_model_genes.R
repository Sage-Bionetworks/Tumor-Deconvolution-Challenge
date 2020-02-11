EPIC::TRef$refProfiles %>% 
    data.frame() %>% 
    tibble::rownames_to_column("Hugo") %>% 
    dplyr::as_tibble() %>% 
    tidyr::pivot_longer(- "Hugo", names_to = "cell_type") %>% 
    dplyr::select(-value) %>% 
    dplyr::mutate(Method = "Epic", from = "tref")


EPIC::BRef$refProfiles %>% 
    data.frame() %>% 
    tibble::rownames_to_column("Hugo") %>% 
    dplyr::as_tibble() %>% 
    tidyr::pivot_longer(- "Hugo", names_to = "cell_type") %>% 
    dplyr::select(-value) %>% 
    dplyr::mutate(Method = "Epic", from = "bref")

xCell::xCell.data$genes %>% 
    dplyr::tibble("Hugo" = .) %>% 
    dplyr::mutate(Method = "xCell", from = "xCell.data$genes")
