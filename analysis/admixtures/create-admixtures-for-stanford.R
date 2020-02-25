library(openxlsx)
file <- "Copy\ of\ Cells_CellLines-DeconvolutionProject\ for\ Bioanalyzer.xlsx"

tbl <- read.xlsx(file, sheet=1)

tbl$`Standard Description`
