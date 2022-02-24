suppressPackageStartupMessages(library("pacman"))
suppressPackageStartupMessages(p_load("synapser"))
suppressPackageStartupMessages(p_load("openxlsx"))


synLogin()


rename.samples <- function(mat) {
    lst <- list(
        "Naive_B_cells" = "B_naive",
        "Macrophages" = "Macro",
        "Dendritic_cells" = "DC",
        "Endothelial_cells" = "Endo",
        "Memory_CD8_T_cells" = "CD8T_mem",
        "Naive_CD8_T_cells" = "CD8T_naive",
        "Memory_CD4_T_cells" = "CD4T_mem",
        "Naive_CD4_T_cells" = "CD4T_naive",
        "Monocytes" = "Mono",
        "Neutrophils" = "Neutro",
        "NK_cells" = "NK")
    for(entry in names(lst)) {
      colnames(mat) <- gsub(colnames(mat), pattern = entry, replacement = lst[[entry]])
    }
    mat
}

## Our original admixture specification includes the vendor for each sample
synId <- "syn21577258"
obj <- synGet(synId, downloadFile = TRUE)
vendors1 <- read.xlsx(obj$path, sheet = 1)[3,,drop=F]
vendors2 <- read.xlsx(obj$path, sheet = 2)[3,,drop=F]
vendors1 <- vendors1[,-1,drop=F]
vendors2 <- vendors2[,-1,drop=F]
exclude.cols <- c("Tregs", "Endothelial_cells", "breast", "CRC", "Fibroblasts")
colnames(vendors1) <-
    unlist(lapply(colnames(vendors1),
                  function(str) ifelse(str %in% exclude.cols,
                                       str,
                                       paste0(str, "_1"))))
colnames(vendors2) <-
    unlist(lapply(colnames(vendors2),
                  function(str) ifelse(str %in% exclude.cols,
                                       str,
                                       paste0(str, "_2"))))
vendors <- cbind(vendors1, vendors2[, !(colnames(vendors2) %in% exclude.cols), drop=F])
vendors <- t(rename.samples(vendors))
rownames(vendors)[grepl(rownames(vendors), pattern="breast")] <- "Breast"
colnames(vendors) <- c("vendor")
vendors <- data.frame(sample = rownames(vendors), vendors)

file <- "sample-vendor-assignment.tsv"
write.table(vendors, file = file, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

syn.folder <- "syn21557719"
f <- File(file, parent = syn.folder)
synStore(f)