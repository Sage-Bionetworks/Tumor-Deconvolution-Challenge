
###############################################################################
## This code uses timer to calculate cell type predictions for
## the coarse-grained sub-Challenge.
###############################################################################


library(magrittr)
library(immunedeconv)

## Read in the round and sub-Challenge-specific input file 
## listing each of the datasets
print(list.files())
print(getwd())
input_df <- readr::read_csv("input/input.csv")

## Extract the names of each dataset
dataset_names <- input_df$dataset.name

## Extract the names of the expression files that use 
## Hugo symbols as gene identifiers
expression_files  <- input_df$hugo.expr.file

## Form the paths of the expression files
expression_paths <- paste0("input/", expression_files)

cancer_types <- input_df %>% 
    dplyr::mutate(
        cancer.type = ifelse(cancer.type == "BRCA", "brca", "coad")
    ) %>% 
    magrittr::use_series(cancer.type)
    




# timer to challenge cell type names
# translation_df <- tibble::tribble(
#     ~cell.type, ~timer.cell.type,
#     "B.cells", "Bcells",
#     "CD4.T.cells", "CD4_Tcells",
#     "CD8.T.cells", "CD8_Tcells",
#     "NK.cells", "NKcells",
#     "neutrophils", "NKcells",
#     "monocytic.lineage", "Macrophages",
#     "fibroblasts", "CAFs",
#     "endothelial.cells", "Endothelial"
# )

## Assumes that expression_path points to a CSV whose gene identifiers
## are HUGO symbols.
do_timer <- function(expression_path, dataset_name, cancer_type){
    x <- expression_path %>%
        readr::read_csv() %>%
        tidyr::drop_na() %>% 
        as.data.frame() %>%
        tibble::column_to_rownames(., colnames(.)[[1]]) %>%
        as.matrix() %>%
        immunedeconv::deconvolute(
            "timer",
            indications = rep(cancer_type, ncol(.))
        ) %>% 
        as.data.frame() %>%
        tibble::rownames_to_column("xcell.cell.type") %>%
        dplyr::as_tibble() %>%
        tidyr::gather(
            key = sample.id,
            value = prediction,
            -xcell.cell.type) %>%
        dplyr::mutate(dataset.name = dataset_name)
}

## Run EPIC on each of the expression files
result_dfs <- purrr::pmap(
    list(expression_paths, dataset_names, cancer_types), 
    do_timer) 

## Combine all results into one dataframe
combined_result_df <- dplyr::bind_rows(result_dfs)

## Translate cell type names as output from MCP-Counter to those
## required for the coarse-grained sub-Challenge.
combined_result_df <- combined_result_df %>% 
    dplyr::inner_join(translation_df) %>% 
    dplyr::select(dataset.name, sample.id, cell.type, prediction)

##### MCPcounter example code above ########

## Create the directory the output will go into
dir.create("output")

## Write result into output directory
readr::write_csv(combined_result_df, "output/predictions.csv")

print(list.files("output"))
    
