

# synapse ---------------------------------------------------------------------

require(synapser)
require(data.table)
require(magrittr)
require(dplyr)
require(stringr)

create_df_from_synapse_id <- function(syn_id, location = NULL, unzip = F, ...){
    path <- download_from_synapse(syn_id, location)
    ## This doesn't work on OSX, where zcat expects the file to have .Z suffix
    ##    if(unzip) path <- stringr::str_c("zcat ", path)
    if(unzip) path <- stringr::str_c("zcat < ", path)
    path %>% 
        data.table::fread(...) %>% 
        dplyr::as_tibble() 
}

download_from_synapse <- function(syn_id, location = NULL){
    path = synapser::synGet(syn_id, downloadLocation = location)$path
    return(path)
}

upload_file_to_synapse <- function(
    path, synapse_id, 
    annotation_list = NULL, 
    activity_obj = NULL, 
    ret = "entity"){
    
    entity <- synapser::File(
        path = path, 
        parent = synapse_id, 
        annotations = annotation_list)
    entity <- synapser::synStore(entity, activity = activity_obj)
    if(ret == "entity") return(entity)
    if(ret == "syn_id") return(entity$properties$id)
}

get_file_df_from_synapse_dir_id <- function(syn_id){
    str_c('select id, name from file where parentId=="', syn_id, '"') %>% 
        synapser::synQuery() %>%
        magrittr::use_series("results") %>% 
        purrr::map(data.frame) %>% 
        dplyr::bind_rows() %>% 
        tibble::as_tibble() 
}


get.synapse.id.of.folder.content_ <- function(folder.synId, file.or.dir.name) {
  children <- synGetChildren(folder.synId)
  children <- as.list(children)
  for(kid in children) {
    if(kid$name == file.or.dir.name) { return(kid$id) }
  }
  return(NULL)
}

get.synapse.id.of.folder.content <- function(folder.synId, file.or.dir.name, exit.on.failure = TRUE) {
  synId <- get.synapse.id.of.folder.content_(folder.synId, file.or.dir.name)
  if(!is.null(synId)) { return(synId) }  
  children <- synGetChildren(folder.synId)
  children <- as.list(children)
  for(kid in children) {
    if(kid$name == file.or.dir.name) { return(kid$id) }
  }
  msg <- paste0("Could not find ", file.or.dir.name, " in ", folder.synId, "\n")
  if(exit.on.failure) {
    stop(msg)
  }
  warn(msg)
}

get.synapse.id <- function(obj) {
  obj$properties$id
}

create.folder <- function(folder.synId, folder.name) {
  synId <- get.synapse.id.of.folder.content_(folder.synId, folder.name)
  if(!is.null(synId)) { return(synId) }
  folder <- Folder(folder.name, parent=folder.synId)
  obj <- synStore(folder)
  synId <- get.synapse.id(obj)
  synId
}

get.synapse.file.location <- function(obj) {
  obj$path
}

upload.data.and.metadata.to.synapse <- function(dataset, expr.mats, gt.mats, mapping.mats, metadata, output.folder.synId, metadata.file.name, executed = NULL, used = NULL, fastq1s = NULL, fastq2s = NULL, fastq1.synIds = NULL, fastq2.synIds = NULL, sample.mapping = NULL) {

  nms <- names(expr.mats)
  names(nms) <- nms
  expr.mat.files <- llply(nms, .fun = function(nm) paste0(dataset, "-", nm, "-gene-expr.csv"))

  expr.mat.synIds <-
    llply(nms,
          .fun = function(nm) {
                   file <- expr.mat.files[[nm]]
  	           mat <- expr.mats[[nm]]
		   sep <- ","
		   if(grepl(file, pattern="tsv")) { sep <- "\t" }
  	           write.table(file = file, mat, sep = sep, col.names = TRUE, row.names = FALSE, quote = FALSE)
    	           f <- File(file, parentId = output.folder.synId, synapseStore = TRUE)
                   ss <- synStore(f, executed = executed, used = used, forceVersion = FALSE)
                   synId <- get.synapse.id(ss)
	           synId
	         })
		 
  for(nm in nms) {
    metadata[[paste0(nm, ".expr.file")]] = expr.mat.files[[nm]]
    metadata[[paste0(nm, ".expr.synId")]] = expr.mat.synIds[[nm]]
##    metadata[[paste0(nm, ".expr.file")]] = "foo.tsv"
##    metadata[[paste0(nm, ".expr.synId")]] = "syn17091415"
  }

  metadata[["n.samples"]] = ncol(expr.mats[[1]])

  nms <- names(gt.mats)
  names(nms) <- nms
  gt.mat.files <-
    llply(nms,
          .fun = function(nm) {
	           if((class(gt.mats[[nm]]) == "logical") && is.na(gt.mats[[nm]])) { return(NA) }
                   paste0(dataset, "-", nm, "-gt.csv")
		 })

  gt.mat.synIds <-
    llply(nms,
          .fun = function(nm) {
	           mat <- gt.mats[[nm]]
		   if((class(mat) == "logical") && is.na(mat)) { return(NA) }
                   file <- gt.mat.files[[nm]]
		   sep <- ","
		   if(grepl(file, pattern="tsv")) { sep <- "\t" }
  	           write.table(file = file, mat, sep = sep, col.names = TRUE, row.names = FALSE, quote = FALSE)
    	           f <- File(file, parentId = output.folder.synId, synapseStore = TRUE)
                   ss <- synStore(f, executed = script_url, forceVersion = FALSE)
                   synId <- get.synapse.id(ss)
	           synId
	         })

  for(nm in nms) {
    metadata[[paste0(nm, ".gt.file")]] = gt.mat.files[[nm]]
    metadata[[paste0(nm, ".gt.synId")]] = gt.mat.synIds[[nm]]    
  }

  for(nm in nms) {
    if((class(gt.mats[[nm]]) == "logical") && is.na(gt.mats[[nm]])) {
      metadata[[paste0("n.", nm, ".pops")]] = NA
      metadata[[paste0(nm, ".pops")]] = NA
      metadata[[paste0("n.", nm)]] = NA      
    } else {
      metadata[[paste0("n.", nm, ".pops")]] = length(unique(gt.mats[[nm]]$cell.type))
      metadata[[paste0(nm, ".pops")]] = paste(sort(as.character(unique(gt.mats[[nm]]$cell.type))), collapse=", ")
      metadata[[paste0("n.", nm)]] = nrow(gt.mats[[nm]])
    }
  }

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

  for(nm in nms) {
    metadata[[paste0(nm, ".to.native.mapping.file")]] = mapping.mat.files[[nm]]
    metadata[[paste0(nm, ".to.native.mapping.synId")]] = mapping.mat.synIds[[nm]]    
  }

  metadata[["sample.mapping.file"]] <- NA
  metadata[["sample.mapping.synId"]] <- NA
  if(!is.null(sample.mapping)) {
    sample.mapping.file <- paste0(dataset, "-sample-mapping.tsv")
    mat <- sample.mapping
    write.table(file = sample.mapping.file, mat, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
    f <- File(sample.mapping.file, parentId = output.folder.synId, synapseStore = TRUE)
    ss <- synStore(f, executed = script_url, forceVersion = FALSE)
    sample.mapping.synId <- get.synapse.id(ss)
    metadata[["sample.mapping.file"]] <- sample.mapping.file
    metadata[["sample.mapping.synId"]] <- sample.mapping.synId
  }
  
  metadata[["fastq.samples"]] <- NA
  metadata[["fastq1.files"]] <- NA
  metadata[["fastq2.files"]] <- NA
  metadata[["fastq1.synIds"]] <- NA
  metadata[["fastq2.synIds"]] <- NA
  if(!is.null(fastq1s) && !is.null(fastq2s)) {
    nms1 <- names(fastq1s)
    nms2 <- names(fastq2s)
    if(!(all(nms1 == nms2))) { stop("FASTQ names mismatch\n") }
    nms1 <- names(fastq1.synIds)
    nms2 <- names(fastq2.synIds)
    if(!(all(nms1 == nms2))) { stop("FASTQ names mismatch\n") }
    nms1 <- names(fastq1s)
    nms2 <- names(fastq1.synIds)
    if(!(all(nms1 == nms2))) { stop("FASTQ names mismatch\n") }
    metadata[["fastq.samples"]] <- paste0(nms1, collapse=",")
    metadata[["fastq1.files"]] <- paste0(fastq1s, collapse=",")
    metadata[["fastq2.files"]] <- paste0(fastq2s, collapse=",")        
    metadata[["fastq1.synIds"]] <- paste0(fastq1.synIds, collapse=",")
    metadata[["fastq2.synIds"]] <- paste0(fastq2.synIds, collapse=",")        
  }

  metadata.df <- data.frame("key" = names(metadata), "value" = as.character(metadata))

  mat <- metadata.df
  write.table(file = metadata.file.name, mat, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  f <- File(metadata.file.name, parentId = output.folder.synId, synapseStore = TRUE)
  ss <- synStore(f, executed = script_url, forceVersion = FALSE)
  synId <- get.synapse.id(ss)

  cat(paste0("Wrote metadata file ", metadata.file.name, " to ", synId, "\n"))
}

download.and.parse.dataset.metadata <- function(metadata.synId) {
  obj <- synGet(metadata.synId, downloadFile = TRUE)
  file <- get.synapse.file.location(obj)
  tbl <- read.table(file, sep = "\t", header = TRUE, as.is = TRUE)
  metadata <- as.list(tbl$value)
  names(metadata) <- tbl$key
  metadata
}

get.dataset.name <- function(metadata) {
  metadata[["orig.dataset.name"]]
}

## Download the hugo-based expression matrix and convert to linear scale
## by exponentiating (if needed)
download.hugo.expression.matrix <- function(metadata) {
  synId <- metadata[["hugo.expr.synId"]]
  scale <- metadata[["scale"]]
  print(scale)
  obj <- synGet(synId, downloadFile = TRUE)
  file <- get.synapse.file.location(obj)
  sep <- ","
  if(grepl(file, pattern="tsv")) { sep <- "\t" }
  mat <- read.table(file, sep = sep, header = TRUE, as.is = TRUE)
  if(grepl(scale, pattern="log", ignore.case=TRUE)) {
    pow <- NA
    if(scale %in% c("Log2", "LOG2")) { pow <- 2 }
    if(scale %in% c("Log10", "LOG10")) { pow <- 10 }
    if(is.na(pow)) { stop(paste0("Could not parse base from log: ", scale, "\n")) }
    col <- colnames(mat)[1]
    mat <- mat %>%
      column_to_rownames(var = col) %>%
      raise_to_power(x = pow, power = .) %>%
      rownames_to_column(var = col)
  }
  mat
}

download.coarse.grained.ground.truth <- function(metadata) {
  synId <- metadata[["coarse.gt.synId"]]
  if(is.na(synId)) { return(NA) }
  obj <- synGet(synId, downloadFile = TRUE)
  file <- get.synapse.file.location(obj)
  read.table(file, sep = ",", header = TRUE, as.is = TRUE)
}

download.fine.grained.ground.truth <- function(metadata) {
  synId <- metadata[["fine.gt.synId"]]
  if(is.na(synId)) { return(NA) }
  obj <- synGet(synId, downloadFile = TRUE)
  file <- get.synapse.file.location(obj)
  read.table(file, sep = ",", header = TRUE, as.is = TRUE)
}

# data_frame / matrix ---------------------------------------------------------

require(magrittr)
require(tibble)
require(plyr)
require(purrr)

transpose_df <- function(df, id_column, new_col){
    df %>% 
        df_to_matrix(id_column) %>%  
        t() %>% 
        matrix_to_df(new_col)
}

df_to_matrix <- function(df, id_column){
    df %>% 
        data.frame() %>% 
        tibble::column_to_rownames(id_column) %>% 
        as.matrix()
}

matrix_to_df <- function(matrix, new_col){
    matrix %>% 
        data.frame() %>% 
        tibble::rownames_to_column(new_col) %>% 
        tibble::as_tibble()
}


aggregate_matrix <- function(
    matrix, grouping_df, group_by, group_across, 
    aggregate_fun = mean,
    split_by_cols = T,
    apply_margin = 1,
    combine_fun = cbind,
    parallel = T){
    
    group_list <- create_group_list(grouping_df, group_by, group_across)
    matrix %>%
        split_matrix(group_list, split_by_cols, parallel) %>% 
        plyr::llply(apply, apply_margin, aggregate_fun, .parallel = parallel) %>% 
        purrr::invoke(combine_fun, .) 
}

library(plyr)
## Translate/compress genes from one name space (e.g., probe ids) to another (e.g., symbols)
compressGenes <- function(e, mapping, from.col = "from", to.col = "to", fun = mean)
{
  e$to    <- mapping[match(rownames(e), mapping[, from.col]), to.col]
  e           <- e[!is.na(e$to),]
  e           <- ddply(.data = e, .variables = "to", .fun = function(x){apply(x[,-ncol(x)],2,fun)},.parallel = T)
  rownames(e) <- e$to
  e           <- e[,-1]
  return(e)
}

suppressPackageStartupMessages(library(pacman))

get.gse <- function(gsm.id) {
  suppressPackageStartupMessages(p_load(GEOquery))
  gse <- getGEO(gsm.id, GSEMatrix = FALSE, AnnotGPL = FALSE, parseCharacteristics = FALSE)
  gse@header$series_id
}

get.immunespace.gsms <-	function(dataset) {
  suppressPackageStartupMessages(p_load(ImmuneSpaceR))
  con <- ImmuneSpaceR::CreateConnection(dataset)
  df <- as.data.frame(con$getDataset("gene_expression_files"))
  df <- subset(df, !is.na(geo_accession))
  df$geo_accession
}

get.immunespace.geo.sample.mapping <-	function(dataset) {
  suppressPackageStartupMessages(p_load(ImmuneSpaceR))
  con <- ImmuneSpaceR::CreateConnection(dataset)
  df <- as.data.frame(con$getDataset("gene_expression_files"))
  df <- subset(df, !is.na(geo_accession))
  df <- unique(df[, c("participant_id", "geo_accession")])
  df
}

get.immunespace.expression.metadata <-	function(dataset) {
  suppressPackageStartupMessages(p_load(ImmuneSpaceR))
  con <- ImmuneSpaceR::CreateConnection(dataset)
  df <- as.data.frame(con$getDataset("gene_expression_files"))
  df <- subset(df, !is.na(geo_accession))
  df
}

connect.geo.db <- function(sql.dir) {
  suppressPackageStartupMessages(p_load("GEOmetadb"))
  sfile <- paste0(sql.dir, "/", "GEOmetadb.sqlite")
  if(!file.exists(sfile)) getSQLiteFile(destdir = sql.dir)
  con <- dbConnect(SQLite(), sfile)
  con
}

get.gsm.gse.tbl <- function(sql.con) {
  sql <- paste0("SELECT DISTINCT gse_gsm.gse, gse_gsm.gsm FROM gse_gsm")
  gse.gsm.tbl <- dbGetQuery(sql.con, sql)
  gse.gsm.tbl
}

get.gses <- function(gsm.ids, sql.dir = "..") {
  con <- connect.geo.db(sql.dir)
  gse.gsm.tbl <- get.gsm.gse.tbl(con)
  ret <- unique(subset(gse.gsm.tbl, gsm %in% gsm.ids)$gse)
  dbDisconnect(con)
  ret
}

ensure.consistency.of.gpls <- function(gses) {
  annotations <- llply(gses, .fun = function(gse) gse@annotation)
  if(!(all(annotations == annotations[[1]]))) { stop("Inconsistent GPLs\n") }
}

get.annotation <- function(gses) {
  ensure.consistency.of.gpls(gses)
  getGEO(gses[[1]]@annotation)
}

get.annotations <- function(gses) {
    llply(gses, .fun = function(gse) getGEO(gse@annotation))
}

get.geo.metadata.tbls <- function(gses) {
  llply(gses,
        .fun = function(gse) {
	         gse %>%
		   phenoData() %>%
		   pData() %>%
		   as.data.frame()
		 })
}		

get.geo.metadata.tbl <- function(gses) {
  tbls <- get.geo.metadata.tbls(gses)
  if(length(tbls) != 1) { stop("get.geo.metadata.tbl was only expecting a single GSE\n") }
  tbls[[1]]
}

translate.platform.name <- function(platform.name) {
  if(platform.name == "[HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array") {
      platform.name <- "Affymetrix HG-U133 Plus 2.0"
  } else if(platform.name == "[HG-U133A] Affymetrix Human Genome U133A Array") {
      platform.name <- "Affymetrix HG-U133A"
  } else if(platform.name == "[HT_HG-U133_Plus_PM] Affymetrix HT HG-U133+ PM Array Plate") {
    platform.name <- "Affymetrix HG-U133 Plus PM"
  } else if(platform.name == "[PrimeView] Affymetrix Human Gene Expression Array") {
    platform.name <- "Affymetrix Human Gene PrimeView"  
  } else if(platform.name == "[HuGene-1_0-st] Affymetrix Human Gene 1.0 ST Array [transcript (gene) version]") {
    platform.name <- "Affymetrix Human Gene 1.0 ST"
  } else if(platform.name == "[HuGene-1_1-st] Affymetrix Human Gene 1.1 ST Array [transcript (gene) version]") {
    platform.name <- "Affymetrix Human Gene 1.1 ST"
  } else if(platform.name == "Illumina HumanHT-12 V3.0 expression beadchip") {
    platform.name <- "Illumina HumanHT-12 V3.0"
  } else if(platform.name == "Illumina HumanHT-12 V4.0 expression beadchip") {
    platform.name <- "Illumina HumanHT-12 V4.0"
  } else if(platform.name == "Illumina HiSeq 2000 (Homo sapiens)") {
    platform.name <- "Illumina HiSeq 2000"
  } else if(platform.name == "Illumina HiSeq 4000 (Homo sapiens)") {
    platform.name <- "Illumina HiSeq 4000"
  } else if(platform.name =- "Illumina NovaSeq 6000 (Homo sapiens)") {
    platform.name <- "Illumina NovaSeq 6000"
  } else if(platform.name == "Illumina NextSeq 500 (Homo sapiens)") {
    platform.name <- "Illumina NextSeq 500"      
  } else {
    stop(paste0("Unknown array type ", platform.name))
  }
}

get.geo.platform.name <- function(gses) {
    gpls <- get.annotations(gses)
    platform.names <-
        unique(Reduce("c",
                      llply(gpls,
                            .fun = function(gpl) {
                                translate.platform.name(Meta(gpl)$title)
                            })))
    platform.names
}

get.geo.data.processing <- function(gses) {
    tbls <- get.geo.metadata.tbls(gses)
    data.processing <-
        unique(Reduce("c", lapply(tbls,
                                  function(tbl) {
                                      flag <- grepl(colnames(tbl), pattern="data.processing")
                                      data.processing <- as.character(unlist(unique(tbl[, flag])))
                                      data.processing
                                  })))
    data.processing <- paste(data.processing, collapse = "; ")
    return(data.processing)
  
  data.processing <-
    llply(gses,
          .fun = function(gse) {
                   gse %>%
		     phenoData() %>%
                     pData() %>%
                     dplyr::pull(data_processing) %>%
                     unique() %>%
                     as.character()
	         })
  unique(data.processing)
}

get.log.or.linear.space <- function(data.processing) {
  cat(paste0(dataset, " data processing: ", paste(data.processing, collapse=", "), "\n"))
  type <- NA
  if(all(data.processing %in% c("gcRMA", "RMA"))) {
    cat(paste0("Based on ", paste(data.processing, collapse=", "),
               " processing, assuming data are in log2 space\n"))
    type <- "Linear"	      
  } else {
    stop(paste0("Data were not processed via RMA, but with ",
                paste(data.processing, collapse=", "),
		".  Can not confirm in log2 space\n\n"))
  }
  type
}

get.probe.to.symbol.map.hugene.1.0.st <- function(probes) {
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
    values = probes,
    uniqueRows=TRUE)

  probe.to.symbol.map <- unique(annotLookup[, c("affy_hugene_1_0_st_v1", "external_gene_name")])
  colnames(probe.to.symbol.map) <- c("from", "to")
  probe.to.symbol.map$from <- as.character(probe.to.symbol.map$from)

  probe.to.symbol.map
}

get.probe.to.ensg.map.hugene.1.0.st <- function(probes) {
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
    values = probes,
    uniqueRows=TRUE)

  probe.to.ensg.map <- unique(annotLookup[, c("affy_hugene_1_0_st_v1", "ensembl_gene_id")])
  colnames(probe.to.ensg.map) <- c("from", "to")
  probe.to.ensg.map$from <- as.character(probe.to.ensg.map$from)

  probe.to.ensg.map
}

get.probe.maps.biomart <- function(probes, array.name = "affy_hugene_1_0_st_v1") {
  suppressPackageStartupMessages(p_load(biomaRt))
  mart <- useMart("ENSEMBL_MART_ENSEMBL")
  mart <- useDataset("hsapiens_gene_ensembl", mart)
  annotLookup <- getBM(
    mart = mart,
    attributes = c(
      array.name,
      "ensembl_gene_id",
      "external_gene_name"),
    filter = array.name,
    values = probes,
    uniqueRows=TRUE)

  probe.to.ensg.map <- unique(annotLookup[, c(array.name, "ensembl_gene_id")])
  colnames(probe.to.ensg.map) <- c("from", "to")
  probe.to.ensg.map$from <- as.character(probe.to.ensg.map$from)

  probe.to.symbol.map <- unique(annotLookup[, c(array.name, "external_gene_name")])
  colnames(probe.to.symbol.map) <- c("from", "to")
  probe.to.symbol.map$from <- as.character(probe.to.symbol.map$from)

  list("ensg" = probe.to.ensg.map, "symbol" = probe.to.symbol.map)
}

get.probe.maps.hugene.1.0.st <- function(probes) {
  get.probe.maps.biomart(probes, array.name = "affy_hugene_1_0_st_v1")
}

get.probe.to.symbol.map.hugene.1.1.st <- function(probes) {
  suppressPackageStartupMessages(p_load(AnnotationDbi))
  suppressPackageStartupMessages(p_load(hugene11sttranscriptcluster.db))
  query_df <-
      AnnotationDbi::select(
          hugene11sttranscriptcluster.db,
          keys=keys(hugene11sttranscriptcluster.db,keytype="PROBEID"),
          columns=c("SYMBOL"),
          keytype="PROBEID") %>%
      as_tibble() %>%
      set_colnames(c("from", "to")) %>%
      drop_na() %>%
      as.data.frame() %>%
      mutate(from = as.character(from)) %>%
      mutate(to = as.character(to))

  query_df <- subset(query_df, from %in% probes)
  query_df
}

get.probe.to.ensg.map.hugene.1.1.st <- function(probes) {
  suppressPackageStartupMessages(p_load(AnnotationDbi))
  suppressPackageStartupMessages(p_load(hugene11sttranscriptcluster.db))
  query_df <-
      AnnotationDbi::select(
          hugene11sttranscriptcluster.db,
          keys=keys(hugene11sttranscriptcluster.db,keytype="PROBEID"),
          columns=c("ENSEMBL"),
          keytype="PROBEID") %>%
      as_tibble() %>%
      set_colnames(c("from", "to")) %>%
      drop_na() %>%
      as.data.frame() %>%
      mutate(from = as.character(from)) %>%
      mutate(to = as.character(to))

  query_df <- subset(query_df, from %in% probes)
  query_df
}


get.probe.to.symbol.map <- function(gses) {
  platform <- get.geo.platform.name(gses)
  if(platform == "Affymetrix Human Gene 1.0 ST") {
    stop("Call get.probe.to.symbol.map.hugene.1.0.st\n")
  }
  gpl <- get.annotation(gses)
  tbl <- Table(gpl)
  sym.col <- c("Symbol", "Gene Symbol")
  sym.col <- sym.col[sym.col %in% colnames(tbl)]
  if(length(sym.col) != 1) { stop("Could not find symbol col in: ", paste(colnames(tbl), collapse=", "), "\n") }
  mapping <- tbl[, c("ID", sym.col)]
  colnames(mapping) <- c("from", "to")
  mapping <- na.omit(mapping)
  mapping <- subset(mapping, to != "")
  ## Sometimes mappings have multiple ids per probe.
  ## e.g., 11715138_s_at CGB /// CGB1 /// CGB2 /// CGB5 /// CGB7 /// CGB8
  ## Create a row for each symbol
  re.mapping <-
    ldply(1:nrow(mapping), .parallel = TRUE,
          .fun = function(i) {
                   str <- as.character(mapping[i, "to"])
                   symbols <- unlist(strsplit(str, split="[ ]*///[ ]*"))
		   df <- data.frame(from = as.character(mapping[i, "from"]), to = symbols, stringsAsFactors = FALSE)
                 })
  re.mapping
}

get.probe.to.entrez.map <- function(gses) {
  gpl <- get.annotation(gses)
  tbl <- Table(gpl)  
  entrez.col <- c("Entrez_Gene_ID", "Entrez Gene", "ENTREZ_GENE_ID")
  entrez.col <- entrez.col[entrez.col %in% colnames(tbl)]
  if(length(entrez.col) != 1) { stop("Could not find entrez col in: ", paste(colnames(tbl), collapse=", "), "\n") }
  mapping <- tbl[, c("ID", entrez.col)]
  colnames(mapping) <- c("from", "to")
  mapping <- na.omit(mapping)
  mapping <- subset(mapping, to != "")
  ## Sometimes mappings have multiple ids per probe.
  ## e.g., 11715138_s_at CGB /// CGB1 /// CGB2 /// CGB5 /// CGB7 /// CGB8
  ## Create a row for each symbol
  re.mapping <-
    ldply(1:nrow(mapping), .parallel = TRUE,
          .fun = function(i) {
                   str <- as.character(mapping[i, "to"])
                   symbols <- unlist(strsplit(str, split="[ ]*///[ ]*"))
		   df <- data.frame(from = as.character(mapping[i, "from"]), to = symbols, stringsAsFactors = FALSE)
                 })
  re.mapping
}

get.probe.to.ensg.map <- function(gses) {
  platform <- get.geo.platform.name(gses)
  if(platform == "Affymetrix Human Gene 1.0 ST") {
    stop("Call get.probe.to.ensg.map.hugene.1.0.st\n")  
  }
  gpl <- get.annotation(gses)
  tbl <- Table(gpl)  
  ensembl.col <- c("Ensembl")
  ensembl.col <- ensembl.col[ensembl.col %in% colnames(tbl)]
  if(length(ensembl.col) == 1) {
    mapping <- tbl[, c("ID", ensembl.col)]
  } else {
    suppressPackageStartupMessages(p_load(org.Hs.eg.db))
    entrez.to.ensg <- as.data.frame(org.Hs.egENSEMBL)
    mapping <- get.probe.to.entrez.map(gses)
    mapping <- merge(mapping, entrez.to.ensg, by.x = "to", by.y = "gene_id")
    mapping <- mapping[, c("from", "ensembl_id")]
  }
  colnames(mapping) <- c("from", "to")
  ## Sometimes mappings have multiple ids per probe, not all (or any) of which are ENSG.
  ## e.g., "ENSG00000206282 /// ENSG00000224841 /// ENSG00000228736 /// ENSG00000237441 /// ENSG00000237825 /// OTTHUMG00000013077 /// OTTHUMG00000031072 /// OTTHUMG00000031342 /// OTTHUMG00000149093 /// OTTHUMG00000149608"
  ## Creat a row for each ENSG
  flag <- grepl(pattern="ENSG", mapping$to)
  mapping <- mapping[flag,]
  re.mapping <-
    ldply(1:nrow(mapping), .parallel = TRUE,
          .fun = function(i) {
                   str <- as.character(mapping[i, "to"])
		   reg.res <- gregexpr(pattern="ENSG\\d+", str)[[1]]
                   ensgs <-
		     unlist(llply(1:length(reg.res), function(j) substr(str, reg.res[j], reg.res[j] + attr(reg.res, "match.length")[j] - 1)))
		   df <- data.frame(from = as.character(mapping[i, "from"]), to = ensgs, stringsAsFactors = FALSE)
		   df
                 })
  re.mapping
}

get.ensg.to.sym.map <- function(gene.ids) {
  suppressPackageStartupMessages(p_load(mygene))
  bm <- queryMany(gene.ids, scopes="ensembl.gene", fields=c("symbol"), species="human")
  bm <- as.data.frame(bm[,c("symbol", "query")])
  names(bm) <- c("symbol", "ensg")
  bm <- bm[!is.na(bm$ensg),]
  bm <- bm[!is.na(bm$symbol),]
  bm <- bm[!(bm$symbol %in% c("")),]
  bm <- bm[!(bm$ensg %in% c("")),]
  bm <- bm[, c("ensg", "symbol")]
  colnames(bm) <- c("from", "to")
  bm$to <- as.character(bm$to)
  bm$from <- as.character(bm$from)  
  bm
}

get.symbol.to.ensg.map <- function(symbols) {
  suppressPackageStartupMessages(p_load(mygene))
  dummy <- data.frame(query = symbols, ensembl = NA)
  bm <- tryCatch({queryMany(symbols, scopes="symbol", fields=c("ensembl.gene"), species="human")}, error = function(e) { return(dummy) })
##  bm <- tryCatch({queryMany(symbols, scopes="symbol", fields=c("ensembl.gene"), species="human", return.as = "records")}, error = function(e) { return(dummy) })
##  bm <- ldply(bm.lst, .fun = function(entry) { if(!("ensembl" %in% names(entry)) || !("gene" %in% names(entry$ensembl))) { return(NULL)}; data.frame(query=entry$query, ensembl=entry$ensembl$gene) })
  if(FALSE) {
  flag <- grepl(pattern="ensembl", colnames(bm))
  if(length(which(flag)) != 1) {
    stop(paste0("Could not find ensembl col in: ", paste(colnames(bm), collapse=" "), "\n"))
  }
  ensg.col <- colnames(bm)[flag]
  }
  ensg.col <- "ensembl"
  bm <- bm[, c("query", ensg.col)]
  lst <- bm[,ensg.col]
  names(lst) <- bm$query
  bm <- ldply(lst, .fun = function(comp) data.frame(unlist(comp)))
  colnames(bm) <- c("from", "to")
  bm <- bm[!(bm$to %in% c("")),,drop=F]
  bm <- bm[!is.na(bm$to) & !is.na(bm$from),,drop=F]
  bm$to <- as.character(bm$to)
  bm$from <- as.character(bm$from)  
  bm
}


choose.max.mad.row <- function(mat) {
  if(nrow(mat) == 1) { return(mat) }
  mads <- unlist(apply(mat, 1, mad))
  mat[which.max(mads)[1], ,drop=F]
}

## Select the row that most often has the maximum value
choose.max.row <- function(mat) {
  if(nrow(mat) == 1) { return(mat) }
  max.indices <- unlist(apply(mat, 2, which.max))
  tmp <- as.data.frame(table(max.indices))
  indx <- as.numeric(tmp[which.max(tmp[,2])[1], 1])
  mat[indx, ,drop=F]
}

bidir.duplicated <- function(x) duplicated(x, fromLast = TRUE) | duplicated(x, fromLast = FALSE)

drop.duplicate.columns <- function(mat) {
  if(!any(duplicated(colnames(mat)))) { return(mat) }
  dups <- unique(colnames(mat)[duplicated(colnames(mat))])
  eps <- 10^-3
  for(dup in dups) {
    cols <- colnames(mat) == dup
    tmp <- mat[, cols, drop = FALSE]
    tmp <- tmp - tmp[,1]
    ## if(all(abs(tmp) < eps)) { cat("All duplicated\n") }
    if(!(all(abs(tmp) < eps))) { stop("Not all duplicated\n") }
  }
  mat[, !duplicated(colnames(mat))]
}

aggregate_rows <- function(mat, map, from.col = "from", to.col = "to", fun = "colMeans", parallel = T) {

  map <- map[map[, from.col] %in% rownames(mat), ]
  lst <- dlply(map, .variables = to.col, .fun = function(df) df[, from.col])
  rows <- llply(lst, .parallel = parallel,
                .fun = function(probes) {
                         do.call(fun, list(mat[probes, , drop=FALSE]))
                       })
  mat <- ldply(rows)
  rownames(mat) <- mat[,1]
  mat <- mat[,-1]
  mat
}



split_matrix <- function(matrix, lst, by_cols = T, parallel = T){
    fun <- ifelse(by_cols, 
                  function(names) matrix[,names, drop = F], 
                  function(names) matrix[names, , drop = F])
    plyr::llply(lst, fun, .parallel = parallel)
}


create_group_list <- function(df, group_by, group_across){
    df %>%
        magrittr::extract2(group_by) %>% 
        split(df, .) %>%
        purrr::map(extract2, group_across)
}

## Rename and add columns based on a map,
## where new column new.col is defined as the sum of the
## existing columns in map[[new.col]]
combine.columns <- function(mat, map) {
  for(new.col in names(map)) {
    old.cols <- map[[new.col]]
    missing.cols <- old.cols[!(old.cols %in% colnames(mat))]
    if(length(missing.cols) != 0) {
      stop(paste0("Missing cols ",
           paste(missing.cols, collapse = ", "),
	   " in matrix with cols ",
	   paste(colnames(mat), collapse = ", "), "\n"))
    }
    ## Exclude any columns that are entirely NA
    flag <- unlist(apply(mat[, old.cols, drop = FALSE], 1,
                         function(row) all(is.na(row))))
    if(all(flag)) { next }			 
    cat(paste0("Combining ",
        paste(map[[new.col]], collapse = ", "),
	" into ", new.col, "\n"))
    vec <- as.numeric(unlist(apply(mat[, old.cols, drop = FALSE], 1,
                                   function(row) {
				     if(all(is.na(row))) { return(NA) }
				     return(sum(row, na.rm = FALSE))
				   })))
    mat <- cbind(mat, vec)
    colnames(mat)[ncol(mat)] <- new.col
  }
  new.cols <- names(map)
  new.cols <- new.cols[new.cols %in% colnames(mat)]
  mat <- mat[, new.cols, drop = FALSE]
  mat
}

# plotting --------------------------------------------------------------------


lm_corr_eqn <- function(df, method = "pearson", display.r2 = FALSE, display.pval = FALSE){
    m <- lm(y ~ x, df);
    ct <- cor.test(df$x, df$y, method = method)
    estimate <- as.numeric(ct$estimate)
    if(display.r2 == TRUE) { estimate <- estimate*estimate }
    pval <- ct$p.value
##    cat(paste0("method = ", method, " estimate = ", estimate, " pval = ", pval, "\n"))
    eq <- NULL
    if((method == "pearson") && (display.r2 == TRUE)) { 
      if(display.pval) { 
        eq <- substitute(italic(r)^2~"="~est*","~~italic(p)~"="~pval, 
                         list(est = format(estimate, digits=3, scientific=0),
                              pval = format(pval, digits=3, scientific=0)))
      } else {
        eq <- substitute(italic(r)^2~"="~est, 
                         list(est = format(estimate, digits=3, scientific=0)))

      }
    } else if((method == "pearson") && (display.r2 == FALSE)) {
      if(display.pval) { 
        eq <- substitute(italic(r)~"="~est*","~~italic(p)~"="~pval, 
                         list(est = format(estimate, digits=3, scientific=0),
                              pval = format(pval, digits=3, scientific=0)))
      } else {
        eq <- substitute(italic(r)~"="~est, 
                         list(est = format(estimate, digits=3, scientific=0)))

      }
    } else if((method == "spearman") && (display.r2 == FALSE)) {
      if(display.pval) { 
        eq <- substitute(rho~"="~est*","~~italic(p)~"="~pval, 
                         list(est = format(estimate, digits=3, scientific=0),
                              pval = format(pval, digits=3, scientific=0)))
      } else {
        eq <- substitute(rho~"="~est, 
                         list(est = format(estimate, digits=3, scientific=0)))

      }
    } else {
      stop(paste("lm_corr_eqn does not know how to handle method = ", method,  " display.r2 = ", display.r2, "\n"))
    }
    as.character(as.expression(eq));                 
}

plot.correlation <- function(x, y, labels = NULL, colors = NULL, display.r2 = FALSE, method = "pearson", display.pval = FALSE, xoffset = 0.5, ...) {
  df <- data.frame(x = x, y = y)
  if(!is.null(labels)) {
    df$labels <- labels
  }
  g <- NULL
  if(is.null(labels)) {
    g <- ggplot(df, aes(x = x, y = y))
  } else {
    g <- ggplot(df, aes(x = x, y = y, label = labels))
  }
  if(!is.null(colors)) {
    g <- g + geom_point(aes(colour = colors))
  } else {
    g <- g + geom_point()
  }
  if(!is.null(labels)) {
    g <- g + geom_text(vjust = "inward", hjust = "inward")
##    suppressPackageStartupMessages(p_load(ggrepel))
##    g <- g + geom_text_repel(point.padding = NA, box.padding = 1)
  }
##  g <- g + theme(legend.position="none")
  g <- g + geom_smooth(data = df, aes(x = x, y = y), method='lm')
  x.min <- min(df$x, na.rm=TRUE)
  x.max <- max(df$x, na.rm=TRUE)
  y.min <- min(df$y, na.rm=TRUE)
  y.max <- max(df$y, na.rm=TRUE)

  ylimits <- NULL
if(FALSE) {
  use.ggplot.2.2.1.limit.code <- TRUE
  if(use.ggplot.2.2.1.limit.code) {
    ylimits <- ggplot_build(g)$layout$panel_ranges[[1]]$y.range
    xlimits <- ggplot_build(g)$layout$panel_ranges[[1]]$x.range
  } else {
    ylimits <- ggplot_build(g)$layout$panel_params[[1]]$y.range
    xlimits <- ggplot_build(g)$layout$panel_params[[1]]$x.range
  }
}
  xlimits <- ggplot_build(g)$layout$panel_params[[1]]$x.range
  ylimits <- ggplot_build(g)$layout$panel_params[[1]]$y.range

## to see why geom_text(size = sz) sz is different than in theme see: ratio of 14/5
## https://stackoverflow.com/questions/25061822/ggplot-geom-text-font-size-control/25062509 
##  g <- g + geom_text(x = x.min + 0.5 * (x.max - x.min), y = y.min + 1 * (y.max - y.min), label = lm_corr_eqn(df, method = method, display.r2 = display.r2, display.pval = display.pval), parse=TRUE, ...)
##  g <- g + geom_text(x = x.min + 0.5 * (x.max - x.min), y = y.min + 0.8 * (y.max - y.min), label = lm_corr_eqn(df, method = method, display.r2 = display.r2, display.pval = display.pval), parse=TRUE, ...)
##  g <- g + geom_text(x = x.min + 0.5 * (x.max - x.min), y = 0.8 * ylimits[2], label = lm_corr_eqn(df, method = method, display.r2 = display.r2, display.pval = display.pval), parse=TRUE, ...)
    g <- g + geom_text(x = xlimits[1] + xoffset * (xlimits[2] - xlimits[1]), y = ylimits[1] + 0.8 * (ylimits[2] - ylimits[1]), label = lm_corr_eqn(df, method = method, display.r2 = display.r2, display.pval = display.pval), parse=TRUE, ...)
  sz <- 25
  g <- g +theme(axis.text.x = element_text(size=sz),
             axis.text.y = element_text(size=sz),
             axis.title.x = element_text(size=sz),
             axis.title.y = element_text(size=sz),
##             title = element_text(size=sz),
             plot.title = element_text(hjust = 0.5))
  g
}

plot.all.cell.type.correlations <- function(data, title, x.col = "actual", y.col = "prediction", method = "spearman") {
  p_load(ggpubr)
  g <- ggplot(data, aes_string(x = x.col, y = y.col))
##  cors <- ddply(data, c("cell.type"), .fun = function(df) { cor = round(cor(df[,x.col], df[,y.col]), 2) })
  g <- g + ggtitle(title)
  g <- g + geom_point()
  g <- g + geom_smooth(method = "lm")
  g <- g + facet_wrap(~ cell.type, scales = "free")
##  g <- g + geom_text(data=cors, aes(label=paste("r=", cor, sep="")), x=30, y=4)
  g <- g + stat_cor(method = method)
  g
}

plot.individual.cell.type.correlations <- function(data, title, x.col = "actual", y.col = "prediction", method = "spearman") {
  gs <-
    dlply(data, .variables = c("cell.type"), .parallel = FALSE,
          .fun = function(df) {
	           cell.type <- df$cell.type[1]
		   x <- df[, x.col]
		   y <- df[, y.col]
		   methods <- c("spearman", "pearson")
		   for(m in method) {
                     ct <- cor.test(x, y, method = method)
                     estimate <- as.numeric(ct$estimate)
                     pval <- ct$p.value
                     cat(paste0("\n", title, " cell type = ", cell.type, " method = ", m, " estimate = ", estimate, " pval = ", pval, "\n"))
		   }
	           g <- plot.correlation(x, y, method = method, display.pval = TRUE)
		   g <- g + xlab("Measured") + ylab("Predicted")
                   g <- g + ggtitle(paste0(title, cell.type))
		   return(g)
                   g <- ggplot(df, aes_string(x = x.col, y = y.col))
                   g <- g + geom_point()
                   g <- g + geom_smooth(method = "lm")
                   g
		 })
  gs		   
}



# misc ------------------------------------------------------------------------

require(magrittr)
require(preprocessCore)

## expr.mat has columns Gene and the sample names
## gt.df has columns sample, cell.type, and possibly others 
subset.and.rename.samples <- function(expr.mat, gt.df, obfuscate.sample.names = TRUE) {

  samples_in_common <- intersect(colnames(expr.mat), gt.df$sample)

  map <- data.frame(from = samples_in_common, to = samples_in_common)

  if(obfuscate.sample.names) {
    new_samples_in_common <- paste0("S", 1:length(samples_in_common))
    map <- data.frame(sample = samples_in_common, new.sample = new_samples_in_common)
    gt.df <- merge(gt.df, map) %>%
      dplyr::select(-sample) %>%
      dplyr::rename(sample = new.sample)
    expr.mat <- expr.mat[, c("Gene", samples_in_common)]
    colnames(expr.mat) <- c("Gene", new_samples_in_common)
    map <- data.frame(from = samples_in_common, to = new_samples_in_common)
    samples_in_common <- new_samples_in_common
  }

  expr.mat <- expr.mat[, c("Gene", samples_in_common)]

  gt.df <- gt.df %>%
      filter(sample %in% samples_in_common)

  list("expr" = expr.mat, "gt" = gt.df, "map" = map)
}

map.populations <- function(mat, map) {
  ret <- combine.columns(mat, map) 
  rownames(ret) <- mat$sample.id

  gt.df.fine <- melt(as.matrix(ret)) %>%
    dplyr::rename(sample.id = Var1) %>%
    dplyr::rename(cell.type = Var2) %>%
    dplyr::rename(measured = value) %>%
    dplyr::select(sample.id, cell.type, measured)

  ret
}

translate.genes <- function(mat, map, fun = "choose.max.mad.row") {
  tmp.map <- subset(map, from %in% mat$Gene)
  flag <- as.character(mat$Gene) %in% as.character(map$from)
  mat <- mat[flag, ]
  rownames(mat) <- NULL
##  tmp.map <- tmp.map[!duplicated(tmp.map$to, fromLast = TRUE) & !duplicated(tmp.map$to, fromLast = FALSE), ]
  ret <- mat %>%
    column_to_rownames(var = "Gene") %>%
    aggregate_rows(., tmp.map, fun = fun, parallel = TRUE) %>%
    rownames_to_column(var = "Gene")
  ret
}

map.and.format.populations <- function(mat, map) {
  mat %>%
    map.populations(., map) %>%
    as.matrix() %>%
    melt() %>%
    dplyr::rename(sample.id = Var1) %>%
    dplyr::rename(cell.type = Var2) %>%
    dplyr::rename(measured = value) %>%        
    add_column(dataset.name = obfuscated.dataset) %>%
    dplyr::select(dataset.name, sample.id, cell.type, measured)
}

get_summary_by_matrix_cols <- function(columns, matrix, fn){
    m <- matrix[,columns]
    if(length(columns) == 1) return(m)
    apply(m, 1, fn)
}

calculate_cpm <- function(counts){
    1000000 * (counts/sum(counts))
}

zscore_matrix <- function(matrix){
    matrix %>% 
        apply(1, scale) %>% 
        t %>% 
        magrittr::set_colnames(colnames(matrix))
}

quantile_normalize_matrix <- function(matrix){
    matrix %>% 
        preprocessCore::normalize.quantiles() %>% 
        magrittr::set_colnames(colnames(matrix)) %>% 
        magrittr::set_rownames(rownames(matrix))
}

if(FALSE) {
  geo.sdys <- c("SDY113", "SDY180", "SDY305", "SDY315", "SDY522", "SDY80", "SDY984")
  names(geo.sdys) <- geo.sdys
  
  get.ground.truth <- function(dataset) {
  
      con   <- ImmuneSpaceR::CreateConnection(dataset)
      df <- con$getDataset("fcs_analyzed_result") %>% 
               dplyr::as_tibble() %>% 
               as.data.frame()
  }
  
  gtds <- llply(geo.sdys, .fun = function(ds) get.ground.truth(ds))
  
  for(gd in names(gtds)) {
    print(gd)
    print(unique(gtds[[gd]]$population_name_reported))
  }

}

propagate.relative.population.frequencies <- function(mat, hierarchical.model, pops) {
    l <- hierarchical.model
    while(length(l) > 0) {
        for(i in 1:length(l)) {
	    parent.col <- l[[i]]$parent	    
            if(!parent.col %in% colnames(mat)) { next }
	    relative.col <- l[[i]]$relative
	    if(!relative.col %in% colnames(mat)) {
              stop(paste0("Was expecting ", relative.col, " to be in matrix ",
	                  "since parent = ", parent.col, " was.\n",
			  "Cols = ", paste(colnames(mat), collapse=", "), "\n"))
            }
	    absolute.col <- l[[i]]$absolute
	    vec <- mat[, parent.col] * mat[, relative.col]
	    mat <- cbind(mat, vec)
	    colnames(mat)[ncol(mat)] <- absolute.col
            l[[i]] <- NULL
            break
        }
    }
    mat[, colnames(mat) %in% pops]
}
