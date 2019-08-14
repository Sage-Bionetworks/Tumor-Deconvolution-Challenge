
###############################################################################
## This code uses quantiseq to calculate cell type predictions for
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


# quantiseq to challenge cell type names
translation_df <- tibble::tribble(
    ~cell.type, ~quantiseq.cell.type,
    "memory.B.cells", "Memory B-cells",
    "naive.B.cells", "naive B-cells",
    "memory.CD4.T.cells", "CD4+ memory T-cells",
    "naive.CD4.T.cells", "CD4+ naive T-cells",
    "regulatory.T.cells", "T cell regulatory (Tregs)",
    "memory.CD8.T.cells", "CD8+ Tem",
    "naive.CD8.T.cells", "CD8+ naive T-cells",
    "NK.cells", "NK cell",
    "neutrophils", "Neutrophil",
    "monocytes", "Monocyte",
    "myeloid.dendritic.cells", "Dendritic cell",
    "macrophages", "Macrophage M1", 
    "macrophages", "Macrophage M2",
    "fibroblasts", "Fibroblasts",
    "endothelial.cells", "Endothelial cells"
)

## Assumes that expression_path points to a CSV whose gene identifiers
## are HUGO symbols.
do_quantiseq <- function(expression_path, dataset_name){
    df <- expression_path %>%
        readr::read_csv() %>%
        tidyr::drop_na() %>% 
        as.data.frame() %>%
        tibble::column_to_rownames(., colnames(.)[[1]]) %>%
        as.matrix() %>%
        immunedeconv::deconvolute("quantiseq") %>% 
        dplyr::as_tibble() %>%
        dplyr::rename(quantiseq.cell.type = cell_type) %>% 
        tidyr::gather(
            key = sample.id,
            value = prediction,
            -quantiseq.cell.type) 
    
    dplyr::tibble(
            "quantiseq.cell.type" = c(
                "Memory B-cells",
                "naive B-cells",
                "CD4+ memory T-cells",
                "CD4+ naive T-cells",
                "CD8+ Tem",
                "CD8+ naive T-cells",
                "Fibroblasts", 
                "Endothelial cells"),
            "prediction" = 0
        ) %>% 
        merge(df$sample.id) %>%
        dplyr::distinct() %>% 
        dplyr::rename(sample.id = y) %>% 
        dplyr::bind_rows(df) %>% 
        dplyr::mutate(dataset.name = dataset_name)
    
    
}

## Run Quantiseq on each of the expression files
result_dfs <- purrr::map2(expression_paths, dataset_names, do_quantiseq) 

## Combine all results into one dataframe
combined_result_df <- dplyr::bind_rows(result_dfs)

## Translate cell type names as output from Quantiseq to those
## required for the coarse-grained sub-Challenge.
combined_result_df2 <- combined_result_df %>% 
    dplyr::inner_join(translation_df) %>% 
    dplyr::select(dataset.name, sample.id, cell.type, prediction) %>% 
    dplyr::group_by(dataset.name, sample.id, cell.type) %>% 
    dplyr::summarise(prediction = sum(prediction)) %>% 
    dplyr::ungroup()


## Create the directory the output will go into
dir.create("output")

## Write result into output directory
readr::write_csv(combined_result_df2, "output/predictions.csv")

print(list.files("output"))