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

## stop("Examine the populations reported and data units (i.e., relation of population to base/parent population)\n")

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
print(sub)

ground_truth_df %>%
  select(population_definition_reported, population_name_reported) %>%
  unique() %>%
  arrange(desc(population_name_reported)) %>%
  print()

normalizing.population <- NA
population.value.col <- "population_cell_number"

## Update the following: normalizing.pop
if(FALSE) {
  
if(!is.na(normalizing.population)) {
  normalizing_df <- ground_truth_df %>%
    filter(population_name_reported == normalizing.population) %>%
    mutate(denominator = population_cell_number) %>%
    select(-population_name_reported, -population_definition_reported, -population_cell_number)
  ground_truth_df <- merge(ground_truth_df, normalizing_df) %>%
    mutate(population_fraction = as.numeric(population_cell_number) / as.numeric(denominator))
  population.value.col <- "population_fraction"
}
} # end FALSE

## Update the following:
cols <- c("participant_id", "population_name_reported", population.value.col)
gt.mat.raw <- unique(ground_truth_df[, cols]) %>%
    acast(., formula = participant_id ~ population_name_reported)

if(FALSE) {
## Not doing this because cell unit is cells, not fraction
population.hierarchy <-
  list(
       list("parent" = "CD45+ cells/uL",
            "relative" = "B lym CD19+,Freq. of,WBC CD45+",
            "absolute" = "CD19+.B.cells"),
       list("parent" = "CD19+.B.cells",
            "relative" = "Q2: CD19+,, CD20+,Freq. of,B lym CD19+",
            "absolute" = "CD19+CD20+.B.cells"),
       list("parent" = "CD19+CD20+.B.cells",
            "relative" = "Naive B cell,Freq. of,Q2: CD19+, CD20+",
            "absolute" = "naive.B.cells"),
       list("parent" = "CD19+CD20+.B.cells",
            "relative" = "Memory B cell,Freq. of,Q2: CD19+, CD20+",
            "absolute" = "memory.B.cells")
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

} # end FALSE

## Update the following:
## fine.grained.definitions <-
##   list("memory.B.cells" = c("X"),
##        "naive.B.cells" = c("X"),
##        "memory.CD4.T.cells" = c("X"),
##        "naive.CD4.T.cells" = c("X"),
##        "regulatory.T.cells" = c("X"),
##        "memory.CD8.T.cells" = c("X"),
##        "naive.CD8.T.cells" = c("X"),
##        "NK.cells" = c("X"),
##        "neutrophils" = c("X"),
##        "monocytes" = c("X"),
##        "myeloid.dendritic.cells" = c("X"),
##        "macrophages" = c("X"),
##        "fibroblasts" = c("X"),
##        "endothelial.cells" = c("X")
##       )

fine.grained.definitions <-
  list("naive.B.cells" = c("Naive B cell,Freq. of,Q2: CD19+, CD20+"),
       "memory.B.cells" = c("Memory B cell,Freq. of,Q2: CD19+, CD20+")
      )

## Update the following:
## coarse.grained.definitions <-
##   list("B.cells" = c("X"),
##        "CD4.T.cells" = c("X"),
##        "CD8.T.cells" = c("X"),
##        "NK.cells" = c("X"),
##        "neutrophils" = c("X"),
##        "monocytic.lineage" = c("X"),
##        "fibroblasts" = c("X"),
##        "endothelial.cells" = c("X")
##       )

coarse.grained.definitions <-
  list("B.cells" = c("B lym CD19+,Freq. of,WBC CD45+")
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
  
## Create a data frame holding the ground truth with columns sample, cell.type, and measured
gt.df <- gt.mat.raw %>%
  melt() %>%
  set_colnames(c("sample", "cell.type", "measured"))

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

expr.mat <- drop.duplicate.columns(expr.mat) %>% rownames_to_column(var = "Gene")

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

expr.mat <- expr.mat[, c("Gene", sample.mapping$geo_accession)]
gt.df <- subset(gt.df, sample %in% sample.mapping$participant_id)
gt.df <- merge(gt.df, sample.mapping[, c("participant_id", "geo_accession")], by.x = c("sample"), by.y = c("participant_id")) %>%
  dplyr::select("sample" = "geo_accession", "cell.type", "measured") 

## May need to update the following get.geo.platform.name function
platform <- get.geo.platform.name(gses)
cancer.type <- NA
data.processing <- unlist(get.geo.data.processing(gses))
print(data.processing)
## scale <- get.log.or.linear.space(data.processing)
## Update the following based on data.processing:
scale <- "Linear"
## Update the following:
native.probe.type <- "Probe"
normalization <- "average"

set.seed(1234)
obfuscated.dataset <- paste0("DS", sum(utf8ToInt(dataset)))

## Obfuscate the sample names
obfuscate.sample.names <- TRUE

ret <- subset.and.rename.samples(expr.mat, gt.df, obfuscate.sample.names)
expr.mat <- ret[["expr"]]
gt.df <- ret[["gt"]]
samples.map <- ret[["map"]]

## Spread ground truth into a matrix
gt.mat.raw <- gt.df %>%    
    spread(key = "cell.type", value = "measured") %>%
    dplyr::rename(sample.id = sample)

## Rename and combine ground truth columns into those required for
## coarse- and fine-grained challenges.
gt.df.coarse <- NA
if(length(coarse.grained.definitions) > 0) {
  gt.df.coarse <- map.and.format.populations(gt.mat.raw, coarse.grained.definitions)
}

gt.df.fine <- NA
if(length(fine.grained.definitions) > 0) {
  gt.df.fine <- map.and.format.populations(gt.mat.raw, fine.grained.definitions)
}

## Extract mappings from probe to gene symbol and Ensembl ID.
probe.to.symbol.map <- get.probe.to.symbol.map(gses)
probe.to.ensg.map <- get.probe.to.ensg.map(gses)

## Translate probes to symbols and to ensembl IDs
compression.fun <- "choose.max.mad.row"
## compression.fun <- "colMeans"
## compression.fun <- "choose.max.row"
symbol.compression.fun <- compression.fun
ensg.compression.fun <- compression.fun
expr.mat.symbol <- translate.genes(expr.mat, probe.to.symbol.map, fun = symbol.compression.fun)
expr.mat.ensg <- translate.genes(expr.mat, probe.to.ensg.map, fun = ensg.compression.fun)

expr.mats <- list("native" = expr.mat, "ensg" = expr.mat.ensg, "hugo" = expr.mat.symbol)
gt.mats <- list("fine" = gt.df.fine, "coarse" = gt.df.coarse)
mapping.mats <- list("symbol" = probe.to.symbol.map, "ensg" = probe.to.ensg.map)

metadata <-
  list("dataset.name" = obfuscated.dataset,
       "orig.dataset.name" = dataset,
       "cancer.type" = cancer.type,
       "platform" = platform,
       "scale" = scale,
       "native.probe.type" = native.probe.type,
       "data.processing" = data.processing,
       "normalization" = normalization,
       "symbol.compression.function" = symbol.compression.fun,
       "ensg.compression.function" = ensg.compression.fun)

identifier <- dataset
if(obfuscate.sample.names) {
  identifier <- obfuscated.dataset
}
## NB: the name of the metadata file uses (and _must_ use) the original dataset name
metadata.file.name <- paste0(dataset, "-metadata.tsv")
upload.data.and.metadata.to.synapse(identifier, expr.mats, gt.mats, mapping.mats, metadata,
                                    output.folder.synId, metadata.file.name,
                                    executed = script_url, used = NULL, sample.mapping = samples.map)

