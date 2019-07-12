suppressPackageStartupMessages(p_load(MCPcounter))
suppressPackageStartupMessages(p_load(dplyr))
suppressPackageStartupMessages(p_load(readr))
suppressPackageStartupMessages(p_load(purrr))
suppressPackageStartupMessages(p_load(magrittr))
suppressPackageStartupMessages(p_load(tibble))
suppressPackageStartupMessages(p_load(tidyr))

## Execute MCP-Counter against a dataset.
## Assumes that expression_path points to a CSV whose gene identifiers
## are HUGO symbols.
do_mcpcounter <- function(expression_matrix, dataset_name){
    
    # We are using the HUGO version of the expression file, so this needs to
    # indicate that here. probests and genes are the dataframes created above.
    result_matrix <- MCPcounter::MCPcounter.estimate(
        expression_matrix,
        featuresType = 'HUGO_symbols')

    # Convert the result matrix back to a dataframe
    result_df <- result_matrix %>% 
        as.data.frame() %>% 
        tibble::rownames_to_column("mcpcounter.cell.type") %>% 
        dplyr::as_tibble()
    
    # Stack the predictions into one column
    result_df <- tidyr::gather(
        result_df,
        key = "sample.id", 
        value = "prediction", 
        -mcpcounter.cell.type) 
    
    # Add dataset column
    result_df <- dplyr::mutate(result_df, dataset.name = dataset_name)
}

run.coarse.grained.mcpcounter.model <- function(expression_matrices, dataset_names) {
  ## Create a table that translates the cell types output by
  ## MCP-Counter ('mcpcounter.cell.type' column) to the cell types
  ## required of the course-grained sub-Challenge ('cell.type' column).
  ## Note that MCP-Counter does not predict CD4 T cells.
  translation_df <- tibble::tribble(
      ~cell.type, ~mcpcounter.cell.type,
      "B.cells", "B lineage",
##      "CD4.T.cells", "T cells",
      "CD8.T.cells", "CD8 T cells",
      "NK.cells", "NK cells",
      "neutrophils", "Neutrophils",
      "monocytic.lineage", "Monocytic lineage",
      "fibroblasts", "Fibroblasts",
      "endothelial.cells", "Endothelial cells"
  )

  ## Run MCP-Counter on each of the expression files
  result_dfs <- purrr::map2(expression_matrices, dataset_names, do_mcpcounter) 

  ## Combine all results into one dataframe
  combined_result_df <- dplyr::bind_rows(result_dfs)

  ## Translate cell type names as output from MCP-Counter to those
  ## required for the coarse-grained sub-Challenge.
  combined_result_df <- combined_result_df %>% 
      dplyr::inner_join(translation_df) %>% 
      dplyr::select(dataset.name, sample.id, cell.type, prediction)
}