source("dataset-setup.R")

suppressPackageStartupMessages(p_load(tidyverse))
suppressPackageStartupMessages(p_load(synapser))
suppressPackageStartupMessages(p_load(synapserutils))
suppressPackageStartupMessages(p_load(data.table))
suppressPackageStartupMessages(p_load(magrittr))
suppressPackageStartupMessages(p_load(Biobase))
suppressPackageStartupMessages(p_load(ImmuneSpaceR))
suppressPackageStartupMessages(p_load(reshape2))

## Begin confirmation

dataset <- "SDY180"
file <- "create_challenge_files.R"

## End confirmation

script_url <- paste0(url.base, "/", dataset, "/", file)
activity_name <- "create Challenge expression and ground truth files"

source("../../scripts/utils.R")
synLogin()

con   <- ImmuneSpaceR::CreateConnection(dataset)

## Get the list of gene expression matrices from ImmuneSpaceR
expr_sets <- con$listGEMatrices()$name

## Specific to SDY180: drop one of the gene expression datasets, which has redundant sample names
expr_sets <- expr_sets[!(expr_sets %in% "SDY180_Other_Grp1Saline_Geo")]

expr_immunespace_df <- expr_sets %>% 
    con$getGEMatrix() %>% 
    con$mapSampleNames(EM = .) %>%
    Biobase::exprs() %>%
    matrix_to_df("Hugo") %>%
    gather(key = "sample", value = "expr",  -Hugo) %>% 
    separate(sample, sep = "_", into = c("sample", "time")) %>% 
    filter(time == "d0") %>% 
    select(-time) %>% 
    mutate(sample = str_remove_all(sample, ".180"))

## experimentData(...)

ground_truth_df <- 
    con$getDataset("fcs_analyzed_result") %>% 
    dplyr::as_tibble() %>% 
    mutate(sample = str_sub(participant_id, end = -5)) %>% 
    filter(study_time_collected == 0) %>%
    filter(study_time_collected_unit == "Days") %>%
    arrange(cohort) %>%
    distinct(participant_id, population_name_reported, .keep_all = TRUE) %>%
    as.data.frame()

sample.mapping <- dataset %>%
  get.immunespace.expression.metadata(.) %>%
  filter(study_time_collected == 0) %>%
  mutate(sample = str_sub(participant_id, end = -5)) %>%
  filter(study_time_collected == 0) %>%
  filter(study_time_collected_unit == "Days") %>%
  arrange(cohort) %>%
  distinct(participant_id, .keep_all = TRUE) %>%
  as.data.frame()
  

cols <- c("participant_id", "population_name_reported", "population_cell_number")
gt.mat.raw <- unique(ground_truth_df[, cols]) %>%
    mutate(population_cell_number = as.numeric(population_cell_number)) %>%
    acast(., formula = participant_id ~ population_name_reported)

gsms <- get.immunespace.gsms(dataset)

## gse.ids <- unique(unlist(llply(gsms, .parallel = TRUE, .fun = function(gsm.id) get.gse(gsm.id))))
gse.ids <- get.gses(gsms, sql.dir = "../")
names(gse.ids) <- gse.ids

## Get the GEO expression matrices 
gses <- llply(gse.ids, .parallel = TRUE, .fun = function(gse.id) getGEO(gse.id, GSEMatrix=TRUE))
gses <- unlist(gses)
expr.mats <- llply(gses, .fun = function(gse) exprs(gse) %>% as.data.frame)

## Combine the GEO expression matrices into one matrix
expr.mat <- Reduce("cbind", expr.mats)

expr.mat <- drop.duplicate.columns(expr.mat)

sample.mapping <- subset(sample.mapping, participant_id %in% ground_truth_df$participant_id)
sample.mapping <- subset(sample.mapping, geo_accession %in% colnames(expr.mat))
if(nrow(sample.mapping) == 0) {
  stop("No samples in common between ground truth and expression\n")
}

flag <- bidir.duplicated(sample.mapping$participant_id)
if(any(flag)) {
  print(sample.mapping[flag, ])
  stop("Redudant SDY sample ids\n")
}
if(any(duplicated(sample.mapping$geo_accession))) {
  stop("Redudant GEO sample ids\n")
}

expr.mat <- expr.mat[, sample.mapping$geo_accession]
gt.mat.raw <- gt.mat.raw[sample.mapping$participant_id, ]
rownames(gt.mat.raw) <- sample.mapping$geo_accession

## Obfuscate the sample names
sample.mapping$id <- paste0("S", 1:nrow(sample.mapping))
rownames(gt.mat.raw) <- sample.mapping$id
colnames(expr.mat) <- sample.mapping$id

fine.grained.definitions <-
  list("naive.CD4.T.cells" = c("Naive_CD4"),
       "memory.CD4.T.cells" = c("CM_CD4", "EM_CD4"),
       "naive.CD8.T.cells" = c("Naive_CD8"),
       "memory.CD8.T.cells" = c("CM_CD8", "EM_CD8"),
       "naive.B.cells" = c("naive_B"),
       "memory.B.cells" = c("IgDp_memory_B", "IgDn_memory_B"),
       "regulatory.T.cells" = c("CD25p_pCD4"),
       "myeloid.dendritic.cells" = c("CD11c_pWBC"),
       "monocytes" = c("CD14p"),
       "neutrophils" = c("Neutros"),
       "NK.cells" = c("CD56br_pLY")
      )

coarse.grained.definitions <-
  list("CD4.T.cells" = c("CD4"),
       "CD8.T.cells" = c("CD8"),
       "B.cells" = c("CD19"),
       "monocytic.lineage" = c("CD14p", "CD11c_pWBC"),
       "neutrophils" = c("Neutros"),
       "NK.cells" = c("CD56br_pLY")
      )

platform <- get.geo.platform.name(gses)
cancer.type <- NA
data.processing <- unlist(get.geo.data.processing(gses))
## scale <- get.log.or.linear.space(data.processing)
scale <- "Linear"
native.probe.type <- "Probe"

obfuscated.dataset <- paste0("DS", sum(utf8ToInt(dataset)))



gt.mat.fine <- combine.columns(gt.mat.raw, fine.grained.definitions)
gt.mat.coarse <- combine.columns(gt.mat.raw, coarse.grained.definitions)

## Extract mappings from probe to gene symbol and Ensembl ID.
probe.to.symbol.map <- get.probe.to.symbol.map(gses)
probe.to.ensg.map <- get.probe.to.ensg.map(gses)

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



# 
# expression_manifest_df <- tibble(
#     path = c("expression_log.tsv", "expression_linear.tsv"),
#     parent = upload_id,
#     executed = script_url,
#     activityName = activity_name,
#     dataset = dataset,
#     used = used,
#     file_type = "expression",
#     expression_type = "microarray",
#     microarray_type = "Illumina Human HT-12 V3 BeadChip",
#     expression_space = c("log2", "linear")
# )
# 
# # is cell number the unit?
# ground_truth_manifest_df <- tibble(
#     path = "ground_truth.tsv",
#     parent = upload_id,
#     executed = script_url,
#     activityName = activity_name,
#     dataset = dataset,
#     used = used,
#     file_type = "ground truth",
#     unit = "cell number",
#     cell_types = str_c(colnames(ground_truth_df)[-1], collapse = ";")
# )
# 
# write_tsv(log_expr_df, "expression_log.tsv")
# write_tsv(linear_expr_df, "expression_linear.tsv")
# write_tsv(ground_truth_df, "ground_truth.tsv")
# 
# write_tsv(expression_manifest_df, "expression_manifest.tsv")
# write_tsv(ground_truth_manifest_df, "ground_truth_manifest.tsv")
# 
# syncToSynapse("expression_manifest.tsv")
# syncToSynapse("ground_truth_manifest.tsv")
