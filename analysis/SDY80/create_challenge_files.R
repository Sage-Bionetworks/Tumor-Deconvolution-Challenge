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

## Begin configuration

file <- "create_challenge_files.R"

## End configuration

script_url <- paste0(url.base, "/", dataset, "/", file)
activity_name <- "create Challenge expression and ground truth files"

source("../../scripts/utils.R")
synLogin()

if(!grepl(dataset, pattern="SDY")) {
  stop("Was expecting an SDY dataset\n")
}
num.dataset.chars <- nchar(gsub(dataset, pattern="SDY", replacement=""))

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

ground_truth_df <- 
    con$getDataset("fcs_analyzed_result") %>%
    dplyr::as_tibble()

cols <- c("cell_number_unit", "study_time_collected", "study_time_collected_unit")
names(cols) <- cols
tbls <-
  llply(cols,
        .fun = function(col) {
                 tbl <- as.data.frame(table(ground_truth_df[, col]))
                 print(col)
                 print(tbl)
		 for(val in as.character(tbl[,1])) {
		   sub <- ground_truth_df[ground_truth_df[, col] == val, ]
		   pops <- sort(unique(sub$population_name_reported))
		   ## cat(paste0(col, " = ", val, ": ", paste(pops, collapse = ", "), "\n"))
		   cat(paste0(col, " = ", val, ": ", length(pops), " populations\n"))
		 }
                 tbl
               })

col <- "cell_number_unit"
target.cell.number.unit <- as.character(tbls[[col]][1, which.max(tbls[[col]][, "Freq"])[1]])
col <- "study_time_collected_unit"
target.study.time.collected.unit <- as.character(tbls[[col]][1, which.max(tbls[[col]][, "Freq"])[1]])

cat(paste0("target.cell.number.unit = ", target.cell.number.unit, "\n"))
cat(paste0("target.study.time.collected.unit = ", target.study.time.collected.unit, "\n"))

## Confirm that participated_id ends in .XXX
print(head(ground_truth_df$participant_id))

pops <- sort(unique(ground_truth_df$population_name_reported))
print(pops)
pops <- pops[!grepl(pattern="pSTAT", pops)]
if(length(pops) == 0) {
  stop("All populations are phospho\n")
}

## stop("Examine the populations reported and data units (i.e., relation of population to base/parent population)\n")

## Update the following: filtering by study_time_collected, study_time_collected_unit, cell_number_unit, and removing participant_id
ground_truth_df <- ground_truth_df %>%
    mutate(sample = str_sub(participant_id, end = -(num.dataset.chars + 2))) %>% 
    filter(study_time_collected == 0) %>%
    filter(cell_number_unit == target.cell.number.unit) %>%    
    filter(study_time_collected_unit == target.study.time.collected.unit) %>%
    arrange(cohort) %>%
    distinct(participant_id, population_name_reported, .keep_all = TRUE) %>%
    mutate(population_cell_number = as.numeric(population_cell_number)) %>%
    as.data.frame()

sub <- subset(as.data.frame(ground_truth_df), participant_id == ground_truth_df$participant_id[1])
o <- order(nchar(sub$population_definition_reported))
sub <- sub[o, ]
print(head(sub))

ground_truth_df %>%
  select(population_definition_reported, population_name_reported) %>%
  unique() %>%
  arrange(desc(population_definition_reported)) %>%
  print()

##  arrange(desc(nchar(population_definition_reported))) %>%

population.value.col <- "population_cell_number"

## Update the following:
cols <- c("participant_id", "population_name_reported", population.value.col)
gt.mat.raw <- unique(ground_truth_df[, cols]) %>%
    acast(., formula = participant_id ~ population_name_reported)

## Update the following:
population.hierarchy <-
  list(
       list("parent" = "ID80, CD19+ of viable CD45+ (Total B cells)",
            "relative" = "ID106, IgD+CD27- of CD20+ B cells* (Naive B)",
	    "absolute" = "B.cells"),
       list("parent" = "B.cells",
            "relative" = "ID94, IgD-CD27+ of CD20+ B cells* (IgD-CD27+ memory B)",
	    "absolute" = "IgD-CD27+ memory B"),
       list("parent" = "B.cells",
            "relative" = "ID101, IgD+CD27+ of CD20+ B cells* (IgD+CD27+ memory B)",
	    "absolute" = "IgD+CD27+ memory B"),
       list("parent" = "B.cells",
            "relative" = "ID113, IgD-CD27- of CD20+ B cells* (IgD-CD27- memory B)",
	    "absolute" = "IgD-CD27- memory B"),
       list("parent" = "B.cells",
            "relative" = "ID106, IgD+CD27- of CD20+ B cells* (Naive B)",
	    "absolute" = "naive B"),
       list("parent" = "ID1, CD3+ of viable CD45+ cells (Total T cells)",
            "relative" = "ID4, CD8+ of total T cells",
	    "absolute" = "CD8.T.cells"),
       list("parent" = "CD8.T.cells",
            "relative" = "ID53, CD45RA+ of CD8+ T cells",
	    "absolute" = "CD45RA+.CD8.T.cells"),
       list("parent" = "CD8.T.cells",
            "relative" = "ID56, CD45RA- of CD8+ T cells (CD45RA- memory CD8+ T)",
	    "absolute" = "CD45RA-.memory.CD8.T.cells"),
       list("parent" = "CD45RA+.CD8.T.cells",
            "relative" = "ID55, CD27- of CD45RA+CD8+ T cells (EMRA CD8+ T)",
	    "absolute" = "EMRA.CD8.T.cells"),
       list("parent" = "CD45RA+.CD8.T.cells",
            "relative" = "ID54, CD27+ of CD45RA+CD8+ T cells (Naive CD8+ T)",
	    "absolute" = "naive.CD8.T.cells"),
       list("parent" = "ID1, CD3+ of viable CD45+ cells (Total T cells)",
            "relative" = "ID2, CD4+ of total T cells",
	    "absolute" = "CD4.T.cells"),
       list("parent" = "CD4.T.cells",
            "relative" = "ID59, CD25hi FoxP3+ of CD4+ T cells (Treg)",
	    "absolute" = "Tregs"),
       list("parent" = "CD4.T.cells",
            "relative" = "ID35, CD45RA- of CD4+ T cells (Total memory CD4+ T)",
	    "absolute" = "memory.CD4.T.cells"),
       list("parent" = "CD4.T.cells",
            "relative" = "ID34, CD45RA+ of CD4+ T cells (Naive T)",
	    "absolute" = "naive.CD4.T.cells")
      )

## Update the following:
final.pops <- c("B.cells", "IgD-CD27+ memory B", "IgD+CD27+ memory B", "IgD-CD27- memory B", "naive B",
                "CD8.T.cells", "CD45RA+.CD8.T.cells", "CD45RA-.memory.CD8.T.cells", "EMRA.CD8.T.cells", "naive.CD8.T.cells",
		"CD4.T.cells", "Tregs", "memory.CD4.T.cells", "naive.CD4.T.cells",
		"ID64, CD14+ of viable CD45+ cells (Total Monocytes)")

mat <- gt.mat.raw
if(!is.null(population.hierarchy) && !is.null(final.pops)) {
  mat <- propagate.relative.population.frequencies(mat / 100, population.hierarchy, final.pops)
}

## Update the following:
map <- list("monocytes" = c("ID64, CD14+ of viable CD45+ cells (Total Monocytes)"),
            "B.cells" = c("B.cells"),
	    "memory.B.cells" = c("IgD-CD27+ memory B", "IgD+CD27+ memory B", "IgD-CD27- memory B"),
	    "naive.B.cells" = c("naive B"),
	    "CD8.T.cells" = c("CD8.T.cells"),
	    "memory.CD8.T.cells" = c("CD45RA-.memory.CD8.T.cells", "EMRA.CD8.T.cells"),
	    "naive.CD8.T.cells" = c("naive.CD8.T.cells"),
	    "CD4.T.cells" = c("CD4.T.cells"),
	    "memory.CD4.T.cells" = c("memory.CD4.T.cells"),
	    "naive.CD4.T.cells" = c("naive.CD4.T.cells"),
	    "Tregs" = "Tregs")
if(!is.null(map)) {
  mat <- combine.columns(mat, map)
}

gt.mat.raw <- mat

## Update the following:
fine.grained.definitions <-
  list("memory.B.cells" = c("memory.B.cells"),
       "naive.B.cells" = c("naive.B.cells"),
       "memory.CD4.T.cells" = c("memory.CD4.T.cells"),
       "naive.CD4.T.cells" = c("naive.CD4.T.cells"),
       "regulatory.T.cells" = c("Tregs"),
       "memory.CD8.T.cells" = c("memory.CD8.T.cells"),
       "naive.CD8.T.cells" = c("naive.CD8.T.cells"),
##       "NK.cells" = c("X"),
##       "neutrophils" = c("X"),
       "monocytes" = c("monocytes")
##       "myeloid.dendritic.cells" = c("X"),
##       "macrophages" = c("X"),
##       "fibroblasts" = c("X"),
##       "endothelial.cells" = c("X")
      )

## Update the following:
coarse.grained.definitions <-
  list("B.cells" = c("B.cells"),
       "CD4.T.cells" = c("CD4.T.cells"),
       "CD8.T.cells" = c("CD8.T.cells"),
##       "NK.cells" = c("X"),
##       "neutrophils" = c("X"),
       "monocytic.lineage" = c("monocytes")
##       "fibroblasts" = c("X"),
##       "endothelial.cells" = c("X")
      )

sample.mapping <- dataset %>%
  get.immunespace.expression.metadata(.)

cols <- c("study_time_collected", "study_time_collected_unit")
for(col in cols) {
  print(col)
  print(table(sample.mapping[, col]))
}

## Update the following: filtering by study_time_collected, study_time_collected_unit, and removing participant_id
sample.mapping <- sample.mapping %>%
  filter(study_time_collected == 0) %>%
  mutate(sample = str_sub(participant_id, end = -(num.dataset.chars + 2))) %>%
  filter(study_time_collected_unit == target.study.time.collected.unit) %>%
  arrange(cohort) %>%
  distinct(participant_id, .keep_all = TRUE) %>%
  as.data.frame()
  
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

## May need to update the following get.geo.platform.name function
platform <- get.geo.platform.name(gses)
cancer.type <- NA
data.processing <- unlist(get.geo.data.processing(gses))
print(data.processing)

## stop("Set scale based on data processing\n")

## scale <- get.log.or.linear.space(data.processing)
## Update the following based on data.processing:
scale <- "Log2"
## Update the following:
native.probe.type <- "Probe"

obfuscated.dataset <- paste0("DS", sum(utf8ToInt(dataset)))

gt.mat.fine <- combine.columns(gt.mat.raw, fine.grained.definitions)
gt.mat.coarse <- combine.columns(gt.mat.raw, coarse.grained.definitions)

## Extract mappings from probe to gene symbol and Ensembl ID.
##probe.to.symbol.map <- get.probe.to.symbol.map(gses)
##probe.to.ensg.map <- get.probe.to.ensg.map(gses)

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

