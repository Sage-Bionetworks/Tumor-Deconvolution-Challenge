## Update the following:
source("dataset-setup.R")

suppressPackageStartupMessages(p_load(tidyverse))
suppressPackageStartupMessages(p_load(synapser))
suppressPackageStartupMessages(p_load(synapserutils))
suppressPackageStartupMessages(p_load(data.table))
suppressPackageStartupMessages(p_load(magrittr))
suppressPackageStartupMessages(p_load(Biobase))
suppressPackageStartupMessages(p_load(ImmuneSpaceR))
suppressPackageStartupMessages(p_load(reshape2))
suppressPackageStartupMessages(p_load(GEOquery))

## Begin configuration

file <- "create_challenge_files.R"

## This table was manually extracted from Racle (2017) eLife
ground.truth.synId <- "syn12650217"

expr.synId <- "syn11969378"

## End configuration

script_url <- paste0(url.base, "/", dataset, "/", file)
activity_name <- "create Challenge expression and ground truth files"

source("../../scripts/utils.R")
synLogin()

if(!grepl(dataset, pattern="GSE")) {
  stop("Was expecting an GSE dataset\n")
}

gse.ids <- list(dataset)
names(gse.ids) <- gse.ids

gses <- llply(gse.ids, .parallel = TRUE, .fun = function(gse.id) getGEO(gse.id, GSEMatrix=TRUE))
gses <- unlist(gses)

stop("stop")

## Process ground truth

geo.anno <- gses %>%
  get.geo.metadata.tbl %>%
  dplyr::rename(cell_type = source_name_ch1) %>%
  dplyr::rename(time = "time:ch1") %>%
  dplyr::rename(donor = "donor:ch1") %>%
  mutate(donor = gsub(donor, pattern="Donor: ", replacement="")) %>%
  filter(time == "0 d") %>%
  select(sample, cell_type, donor) 


## End processing ground truth

## Process GEO expression

## Get the GEO expression matrices, combining into one matrix
expr.mats <- llply(gses, .fun = function(gse) exprs(gse) %>% as.data.frame)
expr.mat <- Reduce("cbind", expr.mats)
expr.mat <- drop.duplicate.columns(expr.mat)

## May need to update the following get.geo.platform.name function
platform <- get.geo.platform.name(gses)
cancer.type <- NA
data.processing <- unlist(get.geo.data.processing(gses))
print(data.processing)

stop("Set scale based on data processing\n")

## scale <- get.log.or.linear.space(data.processing)
## Update the following based on data.processing:
scale <- "Log2"
## Update the following:
native.probe.type <- "Probe"

## Extract mappings from probe to gene symbol and Ensembl ID.
probe.to.symbol.map <- get.probe.to.symbol.map(gses)
probe.to.ensg.map <- get.probe.to.ensg.map(gses)

## End processing GEO expression


## Obfuscate the sample names
sample.mapping$id <- paste0("S", 1:nrow(sample.mapping))
rownames(gt.mat.raw) <- sample.mapping$id
colnames(expr.mat) <- sample.mapping$id




set.seed(1234)
obfuscated.dataset <- paste0("DS", sum(utf8ToInt(dataset)))

gt.mat.fine <- combine.columns(gt.mat.raw, fine.grained.definitions)
gt.mat.coarse <- combine.columns(gt.mat.raw, coarse.grained.definitions)

suppressPackageStartupMessages(p_load(biomaRt))
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)
annotLookup <- getBM(
  mart = mart,
  attributes = c(
    "affy_hugene_1_0_st_v1",
    "ensembl_gene_id",
    "external_gene_name"),
  filter = "affy_hugene_1_0_st_v1",
  values = rownames(expr.mat),
  uniqueRows=TRUE)

probe.to.ensg.map <- unique(annotLookup[, c("affy_hugene_1_0_st_v1", "ensembl_gene_id")])
colnames(probe.to.ensg.map) <- c("from", "to")
probe.to.ensg.map$from <- as.character(probe.to.ensg.map$from)

probe.to.symbol.map <- unique(annotLookup[, c("affy_hugene_1_0_st_v1", "external_gene_name")])
colnames(probe.to.symbol.map) <- c("from", "to")
probe.to.symbol.map$from <- as.character(probe.to.symbol.map$from)

## Collapse probesets to gene symbol and Ensembl ID.
## Represent a gene by the corresponding probe with largest MAD.
expr.mat.symbol <- aggregate_rows(expr.mat, subset(probe.to.symbol.map, from %in% rownames(expr.mat)),
                                  fun = choose.max.mad.row, parallel = TRUE)

expr.mat.ensg <- aggregate_rows(expr.mat, subset(probe.to.ensg.map, from %in% rownames(expr.mat)),
                                fun = choose.max.mad.row, parallel = TRUE)

expr.mats <- list("native" = expr.mat, "ensg" = expr.mat.ensg, "hugo" = expr.mat.symbol)
nms <- names(expr.mats)
names(nms) <- nms
expr.mat.files <- llply(nms, .fun = function(nm) paste0(dataset, "-", nm, "-gene-expr.csv"))

expr.mat.synIds <-
  llply(nms,
        .fun = function(nm) {
                 file <- expr.mat.files[[nm]]
	         mat <- expr.mats[[nm]]
	         write.table(file = file, mat, sep = ",", col.names = TRUE, row.names = TRUE, quote = FALSE)
    	         f <- File(file, parentId = output.folder.synId, synapseStore = TRUE)
                 ss <- synStore(f, executed = script_url, forceVersion = FALSE)
                 synId <- get.synapse.id(ss)
	         synId
	       })

gt.mats <- list("fine" = gt.mat.fine, "coarse" = gt.mat.coarse)
nms <- names(gt.mats)
names(nms) <- nms
gt.mat.files <- llply(nms, .fun = function(nm) paste0(dataset, "-", nm, "-gt.tsv"))

gt.mat.synIds <-
  llply(nms,
        .fun = function(nm) {
                 file <- gt.mat.files[[nm]]
	         mat <- gt.mats[[nm]]
	         write.table(file = file, mat, sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)
    	         f <- File(file, parentId = output.folder.synId, synapseStore = TRUE)
                 ss <- synStore(f, executed = script_url, forceVersion = FALSE)
                 synId <- get.synapse.id(ss)
	         synId
	       })

mapping.mats <- list("symbol" = probe.to.symbol.map, "ensg" = probe.to.ensg.map)
nms <- names(mapping.mats)
names(nms) <- nms
mapping.mat.files <- llply(nms, .fun = function(nm) paste0(dataset, "-", nm, "-to-native-mapping.tsv"))

mapping.mat.synIds <-
  llply(nms,
        .fun = function(nm) {
                 file <- mapping.mat.files[[nm]]
	         mat <- mapping.mats[[nm]]
	         write.table(file = file, mat, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
    	         f <- File(file, parentId = output.folder.synId, synapseStore = TRUE)
                 ss <- synStore(f, executed = script_url, forceVersion = FALSE)
                 synId <- get.synapse.id(ss)
	         synId
	       })

metadata <-
  list("dataset.name" = obfuscated.dataset,
       "orig.dataset.name" = dataset,
       "native.expr.file" = expr.mat.files[["native"]],
       "native.expr.synId" = expr.mat.synIds[["native"]],
       "hugo.expr.file" = expr.mat.files[["hugo"]],
       "hugo.expr.synId" = expr.mat.synIds[["hugo"]],       
       "ensg.expr.file" = expr.mat.files[["ensg"]],
       "ensg.expr.synId" = expr.mat.synIds[["ensg"]],
       "cancer.type" = cancer.type,
       "platform" = platform,
       "scale" = scale,
       "native.probe.type" = native.probe.type,
       "data.processing" = data.processing,
       "coarse.gt.file" = gt.mat.files[["coarse"]],
       "coarse.gt.synId" = gt.mat.synIds[["coarse"]],
       "fine.gt.file" = gt.mat.files[["fine"]],
       "fine.gt.synId" = gt.mat.synIds[["fine"]],
       "symbol.to.native.mapping.file" = mapping.mat.files[["symbol"]],
       "symbol.to.native.mapping.synId" = mapping.mat.synIds[["symbol"]],       
       "ensg.to.native.mapping.file" = mapping.mat.files[["ensg"]],
       "ensg.to.native.mapping.synId" = mapping.mat.synIds[["ensg"]])

metadata.df <- data.frame("key" = names(metadata), "value" = as.character(metadata))

file <- paste0(dataset, "-metadata.tsv")
mat <- metadata.df
write.table(file = file, mat, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
f <- File(file, parentId = output.folder.synId, synapseStore = TRUE)
ss <- synStore(f, executed = script_url, forceVersion = FALSE)
synId <- get.synapse.id(ss)

