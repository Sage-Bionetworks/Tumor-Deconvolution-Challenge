
###############################################################################
## This code uses Cibersort to calculate cell type predictions for
## the coarse-grained sub-Challenge.
###############################################################################


library(dplyr)
library(readr)
library(purrr)
library(magrittr)
library(tibble)
library(tidyr)
library(e1071)
library(parallel)
library(preprocessCore)

print("libraries loaded")

source("CIBERSORT.R")

print("Cibersort loaded")

print(list.files())
print(getwd())


## Read in the round and sub-Challenge-specific input file 
## listing each of the datasets
input_df <- readr::read_csv("input/input.csv")
print(input_df)

## Extract the names of each dataset
dataset_names <- input_df$dataset.name

## Extract the names of the expression files that use 
## Hugo symbols as gene identifiers
expression_files  <- input_df$hugo.expr.file

## Form the paths of the expression files
expression_paths <- paste0("input/", expression_files)

# Cibersort translation table
translation_df <- tibble::tribble(
    ~cell.type, ~cibersort.cell.type,
    "B.cells", "B cells naive",
    "B.cells", "B cells memory",
    "CD4.T.cells", "T cells CD4 naive", 
    "CD4.T.cells", "T cells CD4 memory resting", 
    "CD4.T.cells", "T cells CD4 memory activated",
    "CD4.T.cells", "T cells regulatory (Tregs)", 
    "CD4.T.cells", "T cells follicular helper",
    "CD8.T.cells", "T cells CD8",
    "NK.cells", "NK cells resting", 
    "NK.cells", "NK cells activated",
    "neutrophils", "Neutrophils",
    "monocytic.lineage", "Monocytes",
    "monocytic.lineage", "Macrophages M0",
    "monocytic.lineage", "Macrophages M1",
    "monocytic.lineage", "Macrophages M2",
    "monocytic.lineage", "Dendritic cells resting",
    "monocytic.lineage", "Dendritic cells activated",
    "fibroblasts", "non.immune",
    "endothelial.cells", "non.immune"
)

## Execute Cibersort against a dataset.
## Assumes that expression_path points to a CSV whose gene identifiers
## are HUGO symbols.
do_cibersort <- function(expression_path, dataset_name){
    print(dataset_name)
    
    # files must be in tsv format with no missing values
    expression_path %>% 
        readr::read_csv() %>% 
        tidyr::drop_na() %>% 
        readr::write_tsv("expr.tsv")
    
    result_matrix <- CIBERSORT(
        "LM22.tsv", 
        "expr.tsv", 
        absolute = TRUE,
        abs_method = "sig.score"
    )
    
    # remove expresion file
    file.remove("expr.tsv")
    
    # Convert the result matrix back to a dataframe, remove non cell type, 
    # columns, create non-immune columns
    result_df <- result_matrix %>% 
        as.data.frame() %>% 
        tibble::rownames_to_column("sample.id") %>% 
        dplyr::as_tibble() %>% 
        dplyr::select(-c(
            "P-value", 
            "Correlation", 
            "RMSE", 
            dplyr::starts_with("Absolute score")
        )) %>% 
        dplyr::mutate(non.immune = 1 - rowSums(.[-1])) 
        
    
    # Stack the predictions into one column
    result_df <- tidyr::gather(
        result_df,
        key = cibersort.cell.type, 
        value = prediction, 
        -sample.id) 
    
    # Add dataset column
    result_df <- dplyr::mutate(result_df, dataset.name = dataset_name)
}

## Run Cibersort on each of the expression files
result_dfs <- purrr::map2(expression_paths, dataset_names, do_cibersort) 

## Combine all results into one dataframe
combined_result_df <- dplyr::bind_rows(result_dfs)

## Translate cell type names as output from MCP-Counter to those
## required for the coarse-grained sub-Challenge.
combined_result_df <- combined_result_df %>% 
    dplyr::inner_join(translation_df) %>% 
    dplyr::select(dataset.name, sample.id, cell.type, prediction) %>% 
    dplyr::group_by(dataset.name, sample.id, cell.type) %>% 
    dplyr::summarise_all(mean)

## Create the directory the output will go into
dir.create("output")

## Write result into output directory
readr::write_csv(combined_result_df, "output/predictions.csv")

print(list.files("output"))
    
