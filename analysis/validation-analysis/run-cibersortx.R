suppressPackageStartupMessages(library(pacman))
suppressPackageStartupMessages(p_load(synapser))
suppressPackageStartupMessages(p_load(plyr))
suppressPackageStartupMessages(p_load(dplyr))

suppressPackageStartupMessages(p_load("foreach"))
suppressPackageStartupMessages(p_load("parallel"))

suppressPackageStartupMessages(library("optparse"))

suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(reshape2))

cat(paste0("Run on a machine with Docker installed\n"))
cat(paste0("I ran on tdaws2: 10.23.19.191\n"))

## To re-run validation data
## Rscript ./run-cibersortx.R --token=c615ba7853a9c42acd63b1aaf2231c1a --username=brian.white@sagebase.org
##   --input-file-synId=syn22267272 --input-folder-synId=syn21821096 --output-folder-synId=syn22320184 --prefix=validation

## To run fine-grained in silico admixtures
## Rscript ./run-cibersortx.R --token=c615ba7853a9c42acd63b1aaf2231c1a --username=brian.white@sagebase.org
##   --input-file-synId=syn22332679 --input-folder-synId=syn21647466 --output-folder-synId=syn22331159 --prefix=fine-in-silico-spikeins

## To run coarse-grained in silico admixtures
## Rscript ./run-cibersortx.R --token=c615ba7853a9c42acd63b1aaf2231c1a --username=brian.white@sagebase.org
##   --input-file-synId=syn22332680 --input-folder-synId=syn21647466 --output-folder-synId=syn22331159 --prefix=coarse-in-silico-spikeins

## To run purified admixtures
## Rscript ./run-cibersortx.R --token=c615ba7853a9c42acd63b1aaf2231c1a --username=brian.white@sagebase.org
##   --input-file-synId=syn22331293 --input-folder-synId=syn21782473 --output-folder-synId=syn21576641 --prefix=purified

## To run specificity analysis (which should be same data as "purified admixtures" above):
## Rscript ./run-cibersortx.R --token=c615ba7853a9c42acd63b1aaf2231c1a --username=brian.white@sagebase.org
##   --input-file-synId=syn22392156 --input-folder-synId=syn22392130 --output-folder-synId=syn22725783 --prefix=specificity-coarse

## Rscript ./run-cibersortx.R --token=c615ba7853a9c42acd63b1aaf2231c1a --username=brian.white@sagebase.org
##   --input-file-synId=syn22392155 --input-folder-synId=syn22392130 --output-folder-synId=syn22725783 --prefix=specificity-fine


option_list <- list(
    make_option(c("--prefix"), action="store",
                default=NULL,
                help="Prefix of CIBERSORTx output files"),
    make_option(c("--output-folder-synId"), action="store",
                default=NULL,
                help="Synapse ID of folder in which results should be stored"),
    make_option(c("--input-folder-synId"), action="store",
                default=NULL,
                help="Synapse ID of folder holding input data"),
    make_option(c("--input-file-synId"), action="store",
                default=NULL,
                help="Synapse ID of input.csv file describing data"),
    make_option(c("--token"), action="store",
                default=NULL,
                help="Token to pass to CIBERSORTx"),
    make_option(c("--username"), action="store",
                default="brian.white@sagebase.org",
                help="Username to pass to CIBERSORTx [%default%]")    
)

descr <- "\
Run CIBERSORTx
"

parser <- OptionParser(usage = "%prog [options]", option_list=option_list, description=descr)

arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options

token <- opt$token
if(is.null(token)) { stop("Must pass token to run cibersortX") }

username <- opt$username
if(is.null(username)) { stop("Must pass username to run cibersortX") }

output.folder.synId <- opt$`output-folder-synId`
if(is.null(output.folder.synId)) { stop("Must pass output folder Synapse ID to indicate where results should be stored\n") }

input.folder.synId <- opt$`input-folder-synId`
if(is.null(input.folder.synId)) { stop("Must pass input folder Synapse ID that holds data\n") }

input.file.synId <- opt$`input-file-synId`
if(is.null(input.file.synId)) { stop("Must pass input file Synapse ID that describes datasets\n") }

prefix <- opt$`prefix`
if(is.null(prefix)) { stop("Must specify prefix of output file\n") }


set.seed(1234)

synLogin()

num.cores <- detectCores()
if(!is.na(num.cores) && (num.cores > 1)) {
  suppressPackageStartupMessages(p_load("doMC"))
  cat(paste("Registering ", num.cores-1, " cores.\n", sep=""))
  registerDoMC(cores=(num.cores-1))
}
num.processes <- num.cores - 1

## Run CIBERSORTx on the data from the 'rerun-validation' / 'post-competitive' phase, i.e.,
## where the coarse- and fine-grained challenge use the same data.
## This is as opposed to the original competitive phase (against which the
## methods were ranked) where the coarse- and fine-grained challenges
## differed.

## Read in the "post-competitive" predictions 
obj <- synGet(input.file.synId, downloadFile=TRUE)
input.tbl <- read.table(obj$path, sep=",", header=TRUE, as.is=TRUE, stringsAsFactors=FALSE)

csx.input.dir <- paste0(getwd(), "/csx-input")
csx.output.dir <- paste0(getwd(), "/csx-output")

lm22.matrix.file <- "LM22.txt"
csx.fig.2l.matrix.file <- "cibersort-sig-matrix-2l.txt"
if(!file.exists(lm22.matrix.file)) { stop(paste0(lm22.matrix.file, " does not exist\n")) }
if(!file.exists(csx.fig.2l.matrix.file)) { stop(paste0(csx.fig.2l.matrix.file, " does not exist\n")) }

if(!dir.exists(csx.input.dir)) { dir.create(csx.input.dir) }
if(!dir.exists(csx.output.dir)) { dir.create(csx.output.dir) }

input.files <- input.tbl[, c("dataset.name", "hugo.expr.file")]

children <- synGetChildren(input.folder.synId)
l <- as.list(children)
df <- do.call(rbind.data.frame, l)

df <- merge(df, input.files, by.x = "name", by.y = "hugo.expr.file")

if(!all(input.files$hugo.expr.file %in% df$name)) { stop(paste0("Missing some input files!\n")) }

## Download each of the input files from Synapse and translate to a tsv (as required by CIBERSORTx)
indices <- 1:nrow(df)
csx.input.files <- gsub(df$name, pattern="csv", replacement="tsv")
names(csx.input.files) <- df$dataset.name

l_ply(indices,
      .fun = function(i) {
          synId <- as.character(df[i, "id"])
          name <- as.character(df[i, "name"])
          out.name <- paste0(csx.input.dir, "/", gsub(name, pattern="csv", replacement="tsv"))
	  if(!file.exists(out.name)) {
             obj <- synGet(synId, downloadFile=TRUE, downloadLocation=csx.input.dir)
             tbl <- read.table(obj$path, sep=",", header=TRUE, as.is=TRUE, stringsAsFactors=FALSE)
             write.table(file=out.name, tbl, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
	     system(paste0("rm -f ", obj$path))
	  }
	  
      })

## Run CIBERSORTx on both LM22 and on the signature matrix used in Figure 2L of the CIBERSORTx paper

## Can not run in parallel because cibersortX output files all have the
## same names and would collide with one another across datasets/runs
l_ply(csx.input.files,
      .parallel = FALSE,
      .fun = function(input) {
          ## Run CIBERSORTx on the LM22
          docker.base.cmd <-
              paste0("docker run",
                     " -v ", csx.input.dir, ":/src/data",
                     " -v ", csx.output.dir, ":/src/outdir",
                     " cibersortx/fractions",
                     " --username ", username, 
                     " --token ", token,
                     " --mixture ", input,
                     " --rmbatchBmode TRUE",
                     " --perm 1",
                     " --verbose TRUE",
                     " --QN FALSE")
          print(docker.base.cmd)

          paste0("Running CIBERSORTx with LM22 against ", input, "\n")
          cmd <- paste0(docker.base.cmd,
	                " --sigmatrix ", lm22.matrix.file)
          system(cmd)
	  
	  out.file <- gsub(input, pattern=".tsv",
	                   replacement="-lm22-CIBERSORTx_Adjusted.txt")
          system(paste0("mv -f ",
	                csx.output.dir, "/", "CIBERSORTx_Adjusted.txt ", 
                        csx.output.dir, "/", out.file))

          paste0("Running CIBERSORTx with Fig 2L matrix against ", input, "\n")
          cmd <- paste0(docker.base.cmd,
	                " --sigmatrix ", csx.fig.2l.matrix.file)
          system(cmd)

	  out.file <- gsub(input, pattern=".tsv",
	                   replacement="-fig-2l-CIBERSORTx_Adjusted.txt")
          system(paste0("mv -f ",
	                csx.output.dir, "/", "CIBERSORTx_Adjusted.txt ", 
                        csx.output.dir, "/", out.file))
      })



## HERE

fine.grained.translation.df <- tibble::tribble(
    ~cell.type, ~cibersort.cell.type,
    "memory.B.cells", "B cells memory",
    "naive.B.cells", "B cells naive",
    "memory.CD4.T.cells", "T cells CD4 memory activated",
    "memory.CD4.T.cells", "T cells CD4 memory resting",
    "naive.CD4.T.cells", "T cells CD4 naive",
    "regulatory.T.cells", "T cells regulatory (Tregs)",
    "NK.cells", "NK cells resting",
    "NK.cells", "NK cells activated",
    "neutrophils", "Neutrophils",
    "monocytes", "Monocytes",
    "myeloid.dendritic.cells", "Dendritic cells resting",
    "myeloid.dendritic.cells", "Dendritic cells activated",
    "macrophages", "Macrophages M0",
    "macrophages", "Macrophages M1",
    "macrophages", "Macrophages M2",
    "fibroblasts", "CD10",
    "endothelial.cells", "CD31"
    )

coarse.grained.translation.df <- tibble::tribble(
    ~cell.type, ~cibersort.cell.type,
    "B.cells", "B cells naive",
    "B.cells", "B cells memory",
    "CD4.T.cells", "T cells CD4 naive", 
    "CD4.T.cells", "T cells CD4 memory resting", 
    "CD4.T.cells", "T cells CD4 memory activated",
    "CD4.T.cells", "T cells regulatory (Tregs)", 
    "CD4.T.cells", "T cells follicular helper",
    "CD8.T.cells", "T cells CD8",
    "CD8.T.cells", "T cells gamma delta", 
    "NK.cells", "NK cells resting", 
    "NK.cells", "NK cells activated",
    "neutrophils", "Neutrophils",
    "monocytic.lineage", "Monocytes",
    "monocytic.lineage", "Macrophages M0",
    "monocytic.lineage", "Macrophages M1",
    "monocytic.lineage", "Macrophages M2",
    "monocytic.lineage", "Dendritic cells resting",
    "monocytic.lineage", "Dendritic cells activated",
    "fibroblasts", "CD10",
    "endothelial.cells", "CD31"
)

fine.grained.translation.df$cibersort.cell.type <-
    gsub(fine.grained.translation.df$cibersort.cell.type, pattern=" ", replacement=".")
fine.grained.translation.df$cibersort.cell.type <-
    gsub(fine.grained.translation.df$cibersort.cell.type, pattern="\\(", replacement=".")
fine.grained.translation.df$cibersort.cell.type <-
    gsub(fine.grained.translation.df$cibersort.cell.type, pattern="\\)", replacement=".")
coarse.grained.translation.df$cibersort.cell.type <-
    gsub(coarse.grained.translation.df$cibersort.cell.type, pattern=" ", replacement=".")
coarse.grained.translation.df$cibersort.cell.type <-
    gsub(coarse.grained.translation.df$cibersort.cell.type, pattern="\\(", replacement=".")
coarse.grained.translation.df$cibersort.cell.type <-
    gsub(coarse.grained.translation.df$cibersort.cell.type, pattern="\\)", replacement=".")

combine.files <- function(lm22.res, s2l.res) {

    immune.cols <- colnames(lm22.res)
    ignore.cols <- c("Mixture", "P.value", "Correlation", "RMSE")
    immune.cols <- immune.cols[!(immune.cols %in% ignore.cols)]
    m <- merge(lm22.res[, c("Mixture", immune.cols)],
               s2l.res[, c("Mixture", "CD10", "CD31", "CD45", "EPCAM")],
               by = "Mixture", all = TRUE)
    for(col in immune.cols) {
        ## Scale all of the immune columns in LM22 by the CD45 marker
        m[,col] <- m[,col] * m[,"CD45"]
    }
    colnames(m)[1] <- "sample.id"
    df <- melt(m)
    colnames(df) <- c("sample.id", "cibersort.cell.type", "prediction")
    df$sample.id <- as.character(df$sample.id)
    df$cibersort.cell.type <- as.character(df$cibersort.cell.type)
    df
}

translate.populations <- function(res, translation.df) {
    res %>%
        dplyr::inner_join(translation.df) %>%
        dplyr::select(dataset.name, sample.id, cell.type, prediction) %>% 
        dplyr::group_by(dataset.name, sample.id, cell.type) %>% 
        dplyr::summarise_all(mean)
}

aggregate.cibersortx.results <- function(datasets, lm22.res.files, s2l.res.files) {

    lm22.res.dfs <-
        llply(datasets,
              .fun = function(ds) {
                  read.table(lm22.res.files[[ds]], sep="\t", header=TRUE, stringsAsFactors = FALSE)
              })

    s2l.res.dfs <-
        llply(datasets,
              .fun = function(ds) {
                  read.table(s2l.res.files[[ds]], sep="\t", header=TRUE, stringsAsFactors = FALSE)
              })

    res <-
        ldply(datasets,
              .fun = function(ds) {
                  ret <- combine.files(lm22.res.dfs[[ds]], s2l.res.dfs[[ds]])
                  ret
              })
    colnames(res)[1] <- "dataset.name"
    
    res.fine <- as.data.frame(translate.populations(res, fine.grained.translation.df))
    res.coarse <- as.data.frame(translate.populations(res, coarse.grained.translation.df))
    res.fine$subchallenge <- "fine"
    res.coarse$subchallenge <- "coarse"
    res.fine$.id <- "CIBERSORTx"
    res.coarse$.id <- "CIBERSORTx"
    res.fine$method.name <- "CIBERSORTx"
    res.coarse$method.name <- "CIBERSORTx"
    cols <- c(".id", "dataset.name", "sample.id", "cell.type", "prediction", "subchallenge", "method.name")
    res.fine <- res.fine[, cols]
    res.coarse <- res.coarse[, cols]
    res <- rbind(res.fine, res.coarse)

    res
}

datasets <- names(csx.input.files)
names(datasets) <- datasets

lm22.all.gene.res.files <-
  llply(csx.input.files,
        .fun = function(input) {
	         out.file <- gsub(input, pattern=".tsv",
	                          replacement="-lm22-CIBERSORTx_Adjusted.txt")
                 paste0(csx.output.dir, "/", out.file)
	       })

s2l.all.gene.res.files <-
  llply(csx.input.files,
        .fun = function(input) {
	         out.file <- gsub(input, pattern=".tsv",
	                          replacement="-fig-2l-CIBERSORTx_Adjusted.txt")
                 paste0(csx.output.dir, "/", out.file)
	       })

names(lm22.all.gene.res.files) <- datasets
names(s2l.all.gene.res.files) <- datasets

print(lm22.all.gene.res.files)
print(s2l.all.gene.res.files)

csx.all.gene.res <- aggregate.cibersortx.results(datasets, lm22.all.gene.res.files, s2l.all.gene.res.files)

csx.all.gene.res <- csx.all.gene.res[, !(colnames(csx.all.gene.res) == ".id")]

print(head(csx.all.gene.res))

file <- paste0(prefix, "-csx-all-gene-predictions.tsv")
write.table(file = file, csx.all.gene.res, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

f <- File(file, parentId = output.folder.synId, synapseStore = TRUE)
synStore(f)

cat("Exiting successfully\n")
q(status=0)
