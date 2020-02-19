
###############################################################################
## This code uses xcell to calculate cell type predictions for
## the fine-grained sub-Challenge.
###############################################################################


library(magrittr)
library(xCell)

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


# xcell to challenge cell type names
translation_df <- tibble::tribble(
    ~cell.type, ~xcell.cell.type,
    "memory.B.cells", "Memory B-cells",
    "naive.B.cells", "naive B-cells",
    "memory.CD4.T.cells", "CD4+ memory T-cells",
    "naive.CD4.T.cells", "CD4+ naive T-cells",
    "regulatory.T.cells", "Tregs",
    "memory.CD8.T.cells", "CD8+ Tem",
    "naive.CD8.T.cells", "CD8+ naive T-cells",
    "NK.cells", "NK cells",
    "neutrophils", "Neutrophils",
    "monocytes", "Monocytes",
    "myeloid.dendritic.cells", "DC",
    "macrophages", "Macrophages",
    "fibroblasts", "Fibroblasts",
    "endothelial.cells", "Endothelial cells"
)

## get the scale methods from the input file
scales <- input_df$scale

## get the scale methods from the input file
normalizations <- input_df$normalization

allowed_normalization_methods <- c(
    "CPM", "MAS5", "gcRMA", "RMA", "RMA+quantile normalization+FARMS", 
    "average", "TMM", "RMA+quantile normalization", "normexp", "TPM"
)

## Get cancer types 
cancer_types <- input_df$cancer.type

## Get platforms
platforms <- input_df$platform


## Assumes that expression_path points to a CSV whose gene identifiers
## are HUGO symbols.
do_xcell <- function(
    expression_path, 
    dataset_name,
    scale, 
    normalization, 
    cancer_type,
    platform
){
    if (platform %in% c("Illumina HiSeq 2000", "Illumina NovaSeq")) {
        rnaseq <- TRUE
    } else if (platform %in% c(
        "Affymetrix HG-U133 Plus 2.0", "Affymetrix Human Gene 1.0 ST",
        "Illumina HumanHT-12 V4.0", "Affymetrix Human Gene 1.1 ST",
        "Affymetrix Human Gene PrimeView")) {
        rnaseq <- FALSE
    } else {
        stop("platform not allowed")
    }
    
    # normalization must be one of these methods
    if (!normalization %in% allowed_normalization_methods) {
        stop("non-accepted normalization method")
    }
    
    expression_matrix <- expression_path %>%
        readr::read_csv() %>%
        as.data.frame() %>%
        tibble::column_to_rownames(., colnames(.)[[1]]) %>%
        as.matrix()
    
    if (scale == "Linear") {
        expression_matrix <- expression_matrix
    } else if (scale == "Log2") {
        expression_matrix <- 2^expression_matrix
    } else if (scale == "Log10") {
        expression_matrix <- 10^expression_matrix
    } else {
        stop("non-accepted scale method")
    }
    
    expression_matrix %>% 
        xCell::xCellAnalysis(
            rnaseq = rnaseq,
            cell.types.use = translation_df$xcell.cell.type
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

## Run Xcell on each of the expression files
result_dfs <- purrr::pmap(
    list(expression_paths, dataset_names, scales, 
         normalizations, cancer_types, platforms),
    do_xcell
)

## Combine all results into one dataframe
combined_result_df <- dplyr::bind_rows(result_dfs)

## Translate cell type names as output from Xcell to those
## required for the coarse-grained sub-Challenge.
combined_result_df <- combined_result_df %>% 
    dplyr::inner_join(translation_df) %>% 
    dplyr::select(dataset.name, sample.id, cell.type, prediction)

## Create the directory the output will go into
dir.create("output")

## Write result into output directory
readr::write_csv(combined_result_df, "output/predictions.csv")

print(list.files("output"))
    
