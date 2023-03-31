suppressPackageStartupMessages(library(pacman))
suppressPackageStartupMessages(p_load(synapser))
suppressPackageStartupMessages(p_load(plyr))
suppressPackageStartupMessages(p_load(dplyr))

set.seed(1234)

synLogin()

# Original coarse-grained ground truth
# coarse.gt.synId <- "syn21820375"
# New coarse-grained ground truth based on fine-grained
coarse.gt.synId <- "syn22267267"
fine.gt.synId <- "syn21820376"

merge.predictions.and.ground.truth <- function(preds, ground.truth) {
  preds <- merge(preds, ground.truth, all=TRUE)
  preds <- preds[!is.na(preds$measured),]
}

# Merge results across datasets into one data.frame
# Assume the files are named <method>-<dataset>-{fine,coarse}-predictions.csv
merge.predictions.across.datasets <- function(method, datasets, grain="coarse", ground.truth) {

  nms <- datasets
  names(nms) <- nms

  df <- ldply(nms,
              .fun = function(dataset) {
                       pred.file <- paste0(method, "-", dataset, "-", grain, "-predictions.csv")
                       read.table(pred.file, header=TRUE, sep=",")
                     }) 
  # Don't need to add dataset.name column, it is already in the aginome results.
  df <- df[, -1]
  df <- merge.predictions.and.ground.truth(df, ground.truth)
  df
}

score.preds <- function(preds) {
  # As described in Challenge here
  # Overall score is computed by:
  # 1. Computing correlation/RMSE independently for each dataset/cell type
  # 2. Averaging over cell type for each dataset
  # 3. Averaging over dataset
  res <- 
    ddply(preds, .variables = c("cell.type", "dataset.name"),
          .fun = function(df) {
                   cor.p <- cor(df$prediction, df$measured, method="pearson")
                   cor.s <- cor(df$prediction, df$measured, method="spearman")
                   rmse <- sqrt(mean((df$prediction - df$measured)^2))
                   data.frame(cor.p = cor.p, cor.s = cor.s, rmse = rmse)
                 })

  means.across.cell.types <-
    ddply(res, .variables = c("dataset.name"),
          .fun = function(df) {
                   data.frame(cor.p = mean(df$cor.p), cor.s = mean(df$cor.s), rmse = mean(df$rmse))
                 })

  means.across.datasets <-
    ddply(res, .variables = c("cell.type"),
          .fun = function(df) {
                   data.frame(cor.p = mean(df$cor.p), cor.s = mean(df$cor.s), rmse = mean(df$rmse))
                 })

  overall <- 
    data.frame(cor.p = mean(means.across.cell.types$cor.p),
               cor.s = mean(means.across.cell.types$cor.s),
               rmse = mean(means.across.cell.types$rmse))

  list(all = res, means.across.datasets = means.across.datasets, overall = overall)

}

coarse.gt <- read.table(synGet(coarse.gt.synId)$path, header=TRUE, sep=",")
fine.gt <- read.table(synGet(fine.gt.synId)$path, header=TRUE, sep=",")

datasets <- c("AA", "AB", "AE", "AF", "DS1", "DS2", "DS3", "DS4")
aginome.coarse.preds <- merge.predictions.across.datasets("aginome", datasets, "coarse", coarse.gt)
aginome.fine.preds <- merge.predictions.across.datasets("aginome", datasets, "fine", fine.gt)

aginome.coarse.scores <- score.preds(aginome.coarse.preds)
aginome.fine.scores <- score.preds(aginome.fine.preds)

## Read in the rerun predicitons (i.e., where the coarse- and fine-grained datasets are the same)
synId <- "syn22320329"
obj <- synGet(synId, downloadFile=TRUE)
res.all <- read.table(obj$path, sep=",", header=TRUE, as.is=TRUE, stringsAsFactors=FALSE)
res.all <- subset(res.all, submission == "1")
res.all <- res.all[, c("dataset.name", "subchallenge", "method.name", "cell.type", "sample.id", "prediction")]
coarse.orig.preds <- subset(res.all, subchallenge == "coarse")
fine.orig.preds <- subset(res.all, subchallenge == "fine")

aginome.coarse.orig.preds <- subset(coarse.orig.preds, method.name == "Aginome-XMU")
aginome.coarse.orig.preds <- merge.predictions.and.ground.truth(aginome.coarse.orig.preds, coarse.gt)
aginome.coarse.orig.scores <- score.preds(aginome.coarse.orig.preds)

aginome.fine.orig.preds <- subset(fine.orig.preds, method.name == "Aginome-XMU")
aginome.fine.orig.preds <- merge.predictions.and.ground.truth(aginome.fine.orig.preds, fine.gt)
aginome.fine.orig.scores <- score.preds(aginome.fine.orig.preds)

stop("stop")


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
  ofile <- paste0(files[i, "dataset.name"], "_symbol_tpm.csv")
  write.table(file=ofile, mat, row.names=TRUE, col.names=TRUE, quote=FALSE)
}

stop("stop")

datasets <- names(files)
names(datasets) <- datasets

results <- list()
results[["coarse"]] <- list()
results[["coarse"]][["xmu"]] <- "xmu-coarse-challenge-admixtures-prediction.csv"
results[["fine"]] <- list()
results[["fine"]][["xmu"]] <- "xmu-fine-challenge-admixtures-prediction.csv"

