
###############################################################################
## This code uses Epic to calculate cell type predictions for
## the coarse-grained sub-Challenge.
###############################################################################

library(EPIC)
library(readr)
library(tibble)
library(magrittr)
library(dplyr)
library(purrr)
library(tibble)
library(tidyr)


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


# EPIC to challenge cell type names
translation_df <- tibble::tribble(
    ~cell.type, ~epic.cell.type,
    "B.cells", "Bcells",
    "CD4.T.cells", "CD4_Tcells",
    "CD8.T.cells", "CD8_Tcells",
    "NK.cells", "NKcells",
    "neutrophils", "NKcells",
    "monocytic.lineage", "Macrophages",
    "fibroblasts", "CAFs",
    "endothelial.cells", "Endothelial"
)

## Assumes that expression_path points to a CSV whose gene identifiers
## are HUGO symbols.
do_epic <- function(expression_path, dataset_name){
    
    expression_path %>% 
        readr::read_csv() %>% 
        as.data.frame() %>% 
        tibble::column_to_rownames(., colnames(.)[[1]]) %>% 
        as.matrix() %>% 
        EPIC::EPIC() %>% 
        magrittr::use_series(mRNAProportions) %>% 
        as.data.frame() %>% 
        tibble::rownames_to_column("sample.id") %>% 
        dplyr::as_tibble() %>% 
        tidyr::gather(
            key = epic.cell.type, 
            value = prediction, 
            -sample.id) %>% 
        dplyr::mutate(dataset.name = dataset_name)
}

## Run EPIC on each of the expression files
result_dfs <- purrr::map2(expression_paths, dataset_names, do_epic) 

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
    
