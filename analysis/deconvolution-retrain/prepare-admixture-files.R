suppressPackageStartupMessages(library(pacman))
suppressPackageStartupMessages(p_load(synapser))
suppressPackageStartupMessages(p_load(plyr))
suppressPackageStartupMessages(p_load(dplyr))

synLogin()

# Get the input file for the challenge.
# Note that this is the input using the same admixtures for both fine- and coarse-grained,
# i.e., the 'input-translated-from-fine.csv'
# Originally, we had different admixtures for fine- and for coarse-grained.
input.synId <- "syn22267272"

# Read in the input file
input.obj <- synGet(input.synId, downloadFile=TRUE)
input.tbl <- read.table(input.obj$path, header=TRUE, sep=",")

# Get the folder of the input file, which will hold the data files listed within it.
data.folder.synId <- input.obj$properties$parentId
children <- synGetChildren(data.folder.synId)
l <- as.list(children)
df <- do.call(rbind.data.frame, l)

# Download the data files of interest (e.g., symbol or ensembl id-based, TPM or count-based)
# "hugo.expr.file"                        
# "hugo.expr.est.counts.file"             
# "ensg.expr.file"                        
# "ensg.expr.est.counts.file"       
files <- input.tbl[, c("dataset.name", "hugo.expr.file")]
colnames(files) <- c("dataset.name", "file")

files <- merge(files, df[, c("name", "id")], by.x = c("file"), by.y = c("name"))
for(col in colnames(files)) { files[, col] <- as.character(files[, col]) }

for(i in 1:nrow(files)) {

  # Download the data file
  mat <- read.table(synGet(files[i, "id"], downloadFile=TRUE)$path, sep=",", header=TRUE)

  # Make the Gene column the rowname
  mat$Gene <- as.character(mat$Gene)
  rownames(mat) <- mat$Gene
  mat <- mat[, !(colnames(mat) == "Gene")]
  
  # Save the file
  ofile <- paste0(files[i, "dataset.name"], "_symbol_tpm_gene_as_row.csv")
  write.table(file=ofile, mat, row.names=TRUE, col.names=TRUE, quote=FALSE, sep=",")
}

cat("Done\n")