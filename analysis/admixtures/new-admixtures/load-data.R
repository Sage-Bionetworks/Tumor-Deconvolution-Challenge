library(plyr)
library(tidyverse)

load.iAtlas.supp <- function(cache.dir = ".") {
    library(openxlsx)
    iatlas.file <- paste0(cache.dir, "/iAtlas-sup1.xlsx")
    iatlas.tbl <- read.xlsx(iatlas.file)
    iatlas.tbl
}

## NB: Li has very few counts
load.li.cell.type.counts <- function(cache.dir = ".") {
    ## Li CRC (GSE81861)
    ## synId <- "syn11898281"
    synId <- "syn11898215"  ## this is the expression file
    li.tbl <- read.table(synGet(synId)$path, sep="\t", header=TRUE, comment.char="", quote="\"")
    li.tbl <- li.tbl[, c("sample_id", "cell_type")]
    colnames(li.tbl) <- c("patient", "cell_type")
    li.tbl$patient <- gsub(li.tbl$patient, pattern="_Colorectal_Tumor", replacement="")
    li.tbl <- na.omit(li.tbl)

    tmp <- ddply(li.tbl, .variables = c("patient", "cell_type"),
                 .fun = function(df) data.frame(cnt = nrow(df)))
    tmp <- spread_(tmp[, c("patient", "cell_type", "cnt")], key="cell_type", value="cnt", fill = NA)

    rownames(tmp) <- tmp$patient
    tmp <- tmp[, !(colnames(tmp) == "patient")]
    tmp
}

load.azizi.cell.type.counts <- function(cache.dir = ".") {
    library(openxlsx)
    file <- paste0(cache.dir, "/azizi-fig1d-immune-proportions.xlsx")
    tbl <- read.xlsx(file, startRow = 4)
    cols <- colnames(tbl)
    cols <- cols[!(cols %in% c("Sample", "Total"))]
    tbl[,cols] <- tbl[,cols] / tbl$Total
    tbl <- tbl[, !(colnames(tbl) %in% c("Total"))]

    ## tbl now holds the population fractions
    
    ## Read in the supplemental table that has the total number of cells per sample
    file <- paste0(cache.dir, "/221994-1.xlsx")
    cnt.tbl <- read.xlsx(file, startRow = 3)
    cnt.tbl <- subset(cnt.tbl, sample == "TUMOR")
    cnt.tbl <- ddply(cnt.tbl, .variables = "patient", .fun = function(df) data.frame(n_cells = sum(df$n_cells)))
    tbl <- merge(tbl, cnt.tbl, all.x = TRUE, by.x = "Sample", by.y = "patient")
    tbl[,cols] <- tbl[,cols] * tbl$n_cells
    rownames(tbl) <- tbl$Sample
    tbl <- tbl[, cols]
    tbl <- round(tbl)
    tbl
}

load.cytof.10k.cell.type.fractions <- function(cache.dir = ".") {
  cytof.10k.file <- paste0(cache.dir, "/10KImmunomes.CyTOF_\ PBMC.2018-07-09.csv")
  cytof.tbl <- read.table(cytof.10k.file, sep=",", header=TRUE, as.is=TRUE)
  cytof.tbl <- cytof.tbl[!duplicated(cytof.tbl$subject_accession, fromLast = TRUE) &
                         !duplicated(cytof.tbl$subject_accession, fromLast = FALSE), ]

  cytof.tbl$Memory_CD4_T_cells <- as.numeric(cytof.tbl$Central_Memory_CD4_T_cells) + as.numeric(cytof.tbl$Effector_Memory_CD4_T_cells)
  cytof.tbl$Memory_CD8_T_cells <- as.numeric(cytof.tbl$Central_Memory_CD8_T_cells) + as.numeric(cytof.tbl$Effector_Memory_CD8_T_cells)
  cols <- c("Monocytes", "Tregs", "Naive_CD4_T_cells", "Memory_CD4_T_cells",
            "Naive_CD8_T_cells", "Memory_CD8_T_cells", "NK_cells",
            "Naive_B_cells", "Memory_B_cells")
  cytof.mat <- as.matrix(cytof.tbl[, cols])
  cytof.mat
}

## Can not compare immune (CD45+) and non-immune (CD45-), as these
## were isolated separately.
## In fact, there was further isolation into CD45-/CD90-/CD31- (deplete fibroblasts and endothelial and enrich for malignant)
## and CD45+/CD3+ (to enrich for T cells).
## Hence, really not use to quantitate relative abundance.
load.puram.cell.type.counts <- function(cache.dir = ".") {
    ## Puram Head and Neck (GSE103322)
    ## Get the annotations from the Puram Head and Neck data
    synId <- "syn11990378"
    puram.tbl <- read.table(synGet(synId)$path, sep="\t", header=TRUE)

    ## Exclude myocytes
    puram.tbl <- puram.tbl[!(puram.tbl$cell_type == "myocyte"),]

    ## Clean up some Puram ids
    puram.tbl$sample <- unlist(lapply(puram.tbl$sample, function(x) gsub(x, pattern="^HNSCC_(.+)$", replacement="HNSCC\\1")))
    puram.tbl$sample <- unlist(lapply(puram.tbl$sample, function(x) gsub(x, pattern="^HNSCC(.+)$", replacement="HN\\1")))
    puram.tbl$patient <- unlist(lapply(puram.tbl$sample, function(x) gsub(x, pattern="^([^_]+)_.+$", replacement="\\1")))
    flag <- !grepl(x=puram.tbl$patient, pattern="combo")
    puram.tbl <- puram.tbl[flag, ]
    flag <- grepl(x=puram.tbl$cell_source, pattern="Lymph")
    puram.tbl$patient[flag] <- paste0(puram.tbl$patient[flag], "-LN")
    flag <- grepl(x=puram.tbl$cell_source, pattern="Primary")
    puram.tbl$patient[flag] <- paste0(puram.tbl$patient[flag], "-P")

    tmp <- ddply(puram.tbl, .variables = c("patient", "cell_type"),
                 .fun = function(df) data.frame(cnt = nrow(df)))
    tmp <- spread_(tmp[, c("patient", "cell_type", "cnt")], key="cell_type", value="cnt", fill = NA)

    rownames(tmp) <- tmp$patient
    tmp <- tmp[, !(colnames(tmp) == "patient")]
    tmp
}

## Can not compare immune (CD45+) and non-immune (CD45-), as these
## were isolated separately.
## Hence, can not compare relative abundance of immune vs non-immune.
## In four tumors (Mel58, 67, 72 and 74), we sequenced primarily the immune infiltrates (CD45+ cells) and there were only zero or one malignant cells by this definition.
## This seems also to be true for CY75, which doesn't have any tumor cells.
load.tirosh.cell.type.counts <- function(cache.dir = ".") {
    ## Tirosh Melanoma (GSE72056)
    synId <- "syn12119636"
    tirosh.tbl <- read.table(synGet(synId)$path, sep="\t", header=TRUE)

    tmp <- ddply(tirosh.tbl, .variables = c("patient", "cell_type"),
                 .fun = function(df) data.frame(cnt = nrow(df)))
    tmp <- spread_(tmp[, c("patient", "cell_type", "cnt")], key="cell_type", value="cnt", fill = NA)

    rownames(tmp) <- tmp$patient
    tmp <- tmp[, !(colnames(tmp) == "patient")]
    tmp
}


