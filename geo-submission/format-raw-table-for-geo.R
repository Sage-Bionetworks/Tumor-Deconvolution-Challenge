library(pacman)
p_load(plyr)
p_load(synapser)
p_load(xlsx)
synLogin()

# Get the sample vendor assignment
vendor.info <- read.table(synGet("syn27245024", downloadFile=TRUE)$path, sep="\t", header=TRUE, as.is=TRUE)

# Get the supplemental table with cell type / vendor info
supp.table.synId <- "syn27336433"
supp.tbl <- read.xlsx(synGet(supp.table.synId, downloadFile=TRUE)$path, 1, startRow=2)

supp.tbl <- supp.tbl[, c("Cell.Type", "Vendor", "Vendor.ID", "tissue", "cell.line", "disease.status", "Extract.Protocol")]
colnames(supp.tbl) <- c("characteristics.cell.type", "vendor", "characteristics.vendor.id", "characteristics.tissue", "characteristics.cell.line", "characteristics.disease.status", "characteristics.extract.protocol")
for(col in colnames(supp.tbl)) { supp.tbl[,col] <- as.character(supp.tbl[,col]) }

tbl <- read.table("fastq-table.tsv", sep="\t", header=TRUE, as.is=TRUE)


rename.samples <- function(vec) {
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
      vec <- gsub(vec, pattern = entry, replacement = lst[[entry]])
    }
    vec
}

# tbl$sample <- rename.samples(tbl$sample)

raw.tbl <- tbl

paired.tbl <-
  ddply(tbl, .variables = c("sample"),
        .fun = function(df) {
                 stopifnot(nrow(df) == 2)
                 files <- df$name
                 names(files) <- df$pair
                 r1 <- files[["R1"]]
                 r2 <- files[["R2"]]
                 data.frame(raw.file1 = r1, raw.file2 = r2, type = df[1,"type"])
               })

paired.tbl$Sample.name <- rename.samples(paired.tbl$sample)
# paired.tbl$Sample.name <- paired.tbl$sample

paired.tbl <- merge(paired.tbl, vendor.info, by.x = c("Sample.name"), by.y = c("sample"), all.x=TRUE)

colnames(paired.tbl)[colnames(paired.tbl) == "sample"] <- "title"
paired.tbl$organism <- "Homo sapiens"
paired.tbl$molcule <- "Total RNA"
paired.tbl$processed.data.file1 <- "ensg_counts.csv"
paired.tbl$processed.data.file2 <- "symbol_counts.csv"
paired.tbl$processed.data.file3 <- "transcript_counts.tsv"
paired.tbl$processed.data.file4 <- "ensg_tpm.csv"
paired.tbl$processed.data.file5 <- "symbol_tpm.csv"
paired.tbl$processed.data.file6 <- "transcript_tpms.tsv"

flag <- !is.na(paired.tbl$vendor) & grepl(paired.tbl$vendor, pattern="HIMC")
paired.tbl[flag,"vendor"] <- "Stanford HIMC"

flag <- paired.tbl$type == "Purified Cell"
paired.tbl$characteristics.cell.type <- ""
paired.tbl[flag,"characteristics.cell.type"] <- paired.tbl[flag,"title"]
paired.tbl[flag,"characteristics.cell.type"] <- gsub(paired.tbl[flag,"characteristics.cell.type"], pattern="_1", replacement="")
paired.tbl[flag,"characteristics.cell.type"] <- gsub(paired.tbl[flag,"characteristics.cell.type"], pattern="_2", replacement="")
paired.tbl[flag,"characteristics.cell.type"] <- gsub(paired.tbl[flag,"characteristics.cell.type"], pattern="_", replacement=" ")


paired.tbl <- merge(paired.tbl, supp.tbl, all.x = TRUE, by = c("characteristics.cell.type", "vendor"))

flag <- paired.tbl$type == "Purified Cell"

print(paired.tbl[flag, c("Sample.name","title", "vendor", "characteristics.cell.type", "characteristics.tissue")])

colnames(paired.tbl)[colnames(paired.tbl)=="vendor"] <- "characteristics.vendor"
paired.tbl$source.name <- ""
paired.tbl[flag, "source.name"] <- paste0(paired.tbl$characteristics.cell.type[flag], ", ", paired.tbl$characteristics.disease.status[flag])
paired.tbl$description <- paired.tbl$type
paired.tbl$molecule <- "total RNA"

cols <- c("Sample.name", "title", "source.name", "organism", "characteristics.tissue", "characteristics.cell.type", "characteristics.disease.status",
          "characteristics.cell.line", "characteristics.vendor", "characteristics.vendor.id", "characteristics.extract.protocol", "molecule", 
          "description", "processed.data.file1", "processed.data.file2", "processed.data.file3", "processed.data.file4", "processed.data.file5", 
          "processed.data.file6", "raw.file1", "raw.file2")

paired.tbl[is.na(paired.tbl)] <- ""

paired.tbl <- paired.tbl[,cols]

flag <- paired.tbl$description == "Purified Cell"

mixtures <- paired.tbl[!flag,]
purified <- paired.tbl[flag,]

mixtures <- mixtures[order(mixtures$Sample.name, decreasing=FALSE),]
purified <- purified[order(purified$Sample.name, decreasing=FALSE),]
print(head(purified))
print(head(mixtures))

paired.tbl <- rbind(purified, mixtures)

xlsx.file <- "fastq-files-for-geo.xlsx"
write.xlsx(paired.tbl, file=xlsx.file)
synStore(File(xlsx.file, parent = "syn23395211"))

rownames(raw.tbl) <- raw.tbl$name
raw.tbl$file.type <- "fastq"
paired.tbl$raw.file1 <- as.character(paired.tbl$raw.file1)
paired.tbl$raw.file2 <- as.character(paired.tbl$raw.file2)
ordered.fastqs <- 
  ldply(1:nrow(paired.tbl), .fun = function(i) { data.frame(file = c(paired.tbl[i,"raw.file1"], paired.tbl[i, "raw.file2"])) })
ordered.fastqs$file <- as.character(ordered.fastqs$file)

raw.tbl <- raw.tbl[ordered.fastqs$file,c("name", "file.type", "md5", "sample")]
sample.order <- rename.samples(as.character(raw.tbl$sample))
sample.order <- sample.order[!duplicated(sample.order)]
raw.tbl <- raw.tbl[, c("name", "file.type", "md5")]
xlsx.file <- "fastq-file-md5s-for-geo.xlsx"
write.xlsx(raw.tbl, file=xlsx.file)
synStore(File(xlsx.file, parent = "syn23395211"))

rownames(paired.tbl) <- paired.tbl$Sample.name
paired.tbl <- paired.tbl[sample.order,]

xlsx.file <- "fastq-file-pairs-for-geo.xlsx"
write.xlsx(paired.tbl[, c("raw.file1", "raw.file2")], file=xlsx.file)
synStore(File(xlsx.file, parent = "syn23395211"))

