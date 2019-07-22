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

## This table was manually extracted from Racle (2017) eLife.
## It has the ground truth for the PBMCs--the other populations are pure.
ground.truth.synId <- "syn12650217"

expr.synId <- "syn11969378"

## End configuration

script_url <- paste0(url.base, "/", dataset, "/", file)
activity_name <- "create Challenge expression and ground truth files"

source("../../scripts/utils.R")
synLogin()

## Store data and metadata in this folder: staging/<dataset>
output.folder.synId <- create.folder(staging.folder.synId, dataset)

metadata.file.name <- paste0(dataset, "-metadata.tsv")

if(!grepl(dataset, pattern="GSE")) {
  stop("Was expecting an GSE dataset\n")
}

gse.ids <- list(dataset)
names(gse.ids) <- gse.ids

gses <- llply(gse.ids, .parallel = TRUE, .fun = function(gse.id) getGEO(gse.id, GSEMatrix=TRUE))
gses <- unlist(gses)

## Process ground truth

## Set the proportions of all of the purified cell types to 100% (for that
## samples respective cell type)
geo.anno <- gses %>%
  get.geo.metadata.tbl %>%
  dplyr::rename(cell.type = source_name_ch1) %>%
  dplyr::rename(time = "time:ch1") %>%
  dplyr::rename(donor = "donor:ch1") %>%
  mutate(donor = gsub(donor, pattern="Donor: ", replacement="")) %>%
  mutate(ABV2 = ifelse(cell.type == "myeloid DC", "mDC", 
                         ifelse(cell.type == "Monocytes", "Mono", 
                                ifelse(cell.type == "Neutrophils", "Neut",
				       ifelse(cell.type == "B cells", "B",
				              ifelse(cell.type == "T cells", "T",
					             ifelse(cell.type == "NK cells", "NK",
					                    ifelse(cell.type == "PBMC", "PBMC", cell.type)))))))) %>%
  mutate(days = gsub(x=time, pattern="\\s*(\\d+)\\s+d", replacement="\\1")) %>%
  mutate(sample = str_c(donor, "_", ABV2, "_d", days)) %>% 
  filter(time == "0 d") %>%
  dplyr::select(sample, cell.type) %>%
  mutate(cell.type = as.character(cell.type))  

geo.anno <- subset(geo.anno, cell.type != "PBMC")
geo.anno$measured <- 100

## Get the cell type proportions for the PBMCs
pbmc.gt <- create_df_from_synapse_id(ground.truth.synId) %>%
  melt() %>%
  as.data.frame() %>%
  dplyr::rename(cell.type = variable) %>%
  dplyr::rename(measured = value) %>%
  mutate(cell.type = as.character(cell.type))

## Merge the purified and PBMC proportions
all.gt <- rbind(geo.anno, pbmc.gt)

## Fill in any missing values as zero (i.e., cell types other than the
## purified cell type in the purified samples)
all.mat <- acast(all.gt, formula = sample ~ cell.type, fill = 0)

gt.df <- all.mat %>%
  melt() %>%
  dplyr::rename(sample = Var1) %>%
  dplyr::rename(cell.type = Var2) %>%
  dplyr::rename(measured = value)
  
## End processing ground truth

## Process GEO expression

expr.mat <- expr.synId %>%
    create_df_from_synapse_id(unzip = T, skip = 3) %>% 
    dplyr::select(-c(`Gene Type`, Description, `Gene Symbol`)) %>%
    dplyr::rename(Gene = "Gene ID") %>%
    as.data.frame()

expr.mat <- drop.duplicate.columns(expr.mat) 

## May need to update the following get.geo.platform.name function
platform <- get.geo.platform.name(gses)
cancer.type <- NA
data.processing <- "TMM"
print(data.processing)

## scale <- get.log.or.linear.space(data.processing)
## Update the following based on data.processing:
scale <- "Log"
## Update the following:
native.probe.type <- "ENSG"

## Extract mappings from probe to gene symbol and Ensembl ID.
probe.to.symbol.map <- get.ensg.to.sym.map(as.character(expr.mat$Gene))
map <- data.frame(from = as.character(expr.mat$Gene), to = as.character(expr.mat$Gene))
probe.to.ensg.map <- map

## End processing GEO expression

## Define the coarse- and fine-grained population mapping

fine.grained.definitions <-
  list(
##       "memory.B.cells" = c("Memory_B_cells"),
##       "naive.B.cells" = c("Naive_B_cells"),
##       "memory.CD4.T.cells" = c("Resting_memory_CD4_T_cells", "Activated_memory_CD4_T_cells"),
##       "naive.CD4.T.cells" = c("Naive_CD4_T_cells"),
##       "regulatory.T.cells" = c("Tregs"),
##       "memory.CD8.T.cells" = c("memory.CD8.T.cells"), 
##       "naive.CD8.T.cells" = c("naive.CD8.T.cells"),
       "NK.cells" = c("NK cells"),
       "neutrophils" = c("Neutrophils"),
       "monocytes" = c("Monocytes"),
       "myeloid.dendritic.cells" = c("myeloid DC")
##       "macrophages" = c("X"),
##       "fibroblasts" = c("Fibroblasts"),
##       "endothelial.cells" = c("Endothelial")
      )

## Update the following:
coarse.grained.definitions <-
  list("B.cells" = c("B cells"),
##       "CD4.T.cells" = c("Naive_CD4_T_cells", "Resting_memory_CD4_T_cells", "Activated_memory_CD4_T_cells"),
##       "CD8.T.cells" = c("CD8_T_cells"),
       "NK.cells" = c("NK cells"),
       "neutrophils" = c("Neutrophils"),
       "monocytic.lineage" = c("Monocytes", "myeloid DC")
##       "fibroblasts" = c("Fibroblasts"),
##       "endothelial.cells" = c("Endothelial")
      )
      
## End defining the coarse- and fine-grained population mapping

## May need to update the following get.geo.platform.name function
set.seed(1234)
obfuscated.dataset <- paste0("DS", sum(utf8ToInt(dataset)))

## Obfuscate the sample names
obfuscate.sample.names <- FALSE

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

## Translate probes to symbols and to ensembl IDs
expr.mat.symbol <- translate.genes(expr.mat, probe.to.symbol.map, fun = choose.max.mad.row)

## expr.mat.ensg <- translate.genes(expr.mat, probe.to.ensg.map, fun = choose.max.mad.row)
expr.mat.ensg <- expr.mat

expr.mats <- list("native" = expr.mat, "ensg" = expr.mat.ensg, "hugo" = expr.mat.symbol)
gt.mats <- list("fine" = gt.df.fine, "coarse" = gt.df.coarse)
mapping.mats <- list("symbol" = probe.to.symbol.map, "ensg" = probe.to.ensg.map)

## Get the fastq files, the mappings are stored in the following manifest
fastq.manifest.synId <- "syn12663606"

fastq.manifest <- create_df_from_synapse_id(fastq.manifest.synId) %>%
  mutate(ABV2 = ifelse(cell_type == "myeloid DC", "mDC", 
                         ifelse(cell_type == "Monocytes", "Mono", 
                                ifelse(cell_type == "Neutrophils", "Neut",
				       ifelse(cell_type == "B cells", "B",
				              ifelse(cell_type == "T cells", "T",
					             ifelse(cell_type == "NK cells", "NK",
					                    ifelse(cell_type == "PBMC", "PBMC", cell_type)))))))) %>%
  mutate(sample = str_c(patient, "_", ABV2, "_d", day)) %>% 
  filter(day == 0) %>%
  mutate(fastq1 = str_c(SRR_id, "_1.fastq.gz")) %>%
  mutate(fastq2 = str_c(SRR_id, "_2.fastq.gz")) %>%  
  dplyr::select(sample, cell_type, fastq1, fastq2) 

if(!all(samples.map$from %in% fastq.manifest$sample)) {
  stop("Missing some fastqs\n")
}

fastq.manifest <- subset(fastq.manifest, sample %in% samples.map$from)

fastq1s <- as.character(fastq.manifest$fastq1)
names(fastq1s) <- as.character(fastq.manifest$sample)

fastq2s <- as.character(fastq.manifest$fastq2)
names(fastq2s) <- as.character(fastq.manifest$sample)

metadata <-
  list("dataset.name" = obfuscated.dataset,
       "orig.dataset.name" = dataset,
       "cancer.type" = cancer.type,
       "platform" = platform,
       "scale" = scale,
       "native.probe.type" = native.probe.type,
       "data.processing" = data.processing)

upload.data.and.metadata.to.synapse(dataset, expr.mats, gt.mats, mapping.mats, metadata, output.folder.synId, metadata.file.name,
                                    executed = script_url, used = NULL, fastq1s = fastq1s, fastq2s = fastq2s)


