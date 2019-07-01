

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

get.synapse.id <- function(obj) {
  obj$properties$id
}

get.synapse.file.location <- function(obj) {
  obj$path
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

download.hugo.expression.matrix <- function(metadata) {
  synId <- metadata[["hugo.expr.synId"]]
  obj <- synGet(synId, downloadFile = TRUE)
  file <- get.synapse.file.location(obj)
  read.table(file, sep = ",", header = TRUE, as.is = TRUE)
}

download.coarse.grained.ground.truth <- function(metadata) {
  synId <- metadata[["coarse.gt.synId"]]
  obj <- synGet(synId, downloadFile = TRUE)
  file <- get.synapse.file.location(obj)
  read.table(file, sep = "\t", header = TRUE, as.is = TRUE)
}

download.fine.grained.ground.truth <- function(metadata) {
  synId <- metadata[["fine.gt.synId"]]
  obj <- synGet(synId, downloadFile = TRUE)
  file <- get.synapse.file.location(obj)
  read.table(file, sep = "\t", header = TRUE, as.is = TRUE)
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

get.geo.platform.name <- function(gses) {
  gpl <- get.annotation(gses)
  platform.name <- Meta(gpl)$title
  if(platform.name == "[HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array") {
    platform.name <- "Affymetrix HG-U133 Plus 2.0"
  } else if(platform.name == "[PrimeView] Affymetrix Human Gene Expression Array") {
    platform.name <- "Affymetrix Human Gene PrimeView"  
  } else if(platform.name == "[HuGene-1_0-st] Affymetrix Human Gene 1.0 ST Array [transcript (gene) version]") {
    platform.name <- "Affymetrix Human Gene 1.0 ST"
  } else if(platform.name == "Illumina HumanHT-12 V3.0 expression beadchip") {
    platform.name <- "Illumina HumanHT-12 V3.0"
  } else {
    stop(paste0("Unknown array type ", platform.name))
  }
  platform.name
}

get.geo.data.processing <- function(gses) {
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

get.probe.to.symbol.map <- function(gses) {
  gpl <- get.annotation(gses)
  tbl <- Table(gpl)
  sym.col <- c("Symbol", "Gene Symbol")
  sym.col <- sym.col[sym.col %in% colnames(tbl)]
  if(length(sym.col) != 1) { stop("Could not find symbol col in: ", paste(colnames(tbl), collapse=", "), "\n") }
  mapping <- tbl[, c("ID", sym.col)]
  colnames(mapping) <- c("from", "to")
  mapping <- na.omit(mapping)
  mapping <- subset(mapping, to != "")
  mapping
}

get.probe.to.entrez.map <- function(gses) {
  gpl <- get.annotation(gses)
  tbl <- Table(gpl)  
  entrez.col <- c("Entrez_Gene_ID", "Entrez Gene")
  entrez.col <- entrez.col[entrez.col %in% colnames(tbl)]
  if(length(entrez.col) != 1) { stop("Could not find entrez col in: ", paste(colnames(tbl), collapse=", "), "\n") }
  mapping <- tbl[, c("ID", entrez.col)]
  colnames(mapping) <- c("from", "to")
  mapping <- na.omit(mapping)
  mapping <- subset(mapping, to != "")
  mapping
}

get.probe.to.ensg.map <- function(gses) {
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
  mapping
}

choose.max.mad.row <- function(mat) {
  if(nrow(mat) == 1) { return(mat) }
  mads <- unlist(apply(mat, 1, mad))
  mat[which.max(mads)[1], ,drop=F]
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

aggregate_rows <- function(mat, map, from.col = "from", to.col = "to", fun = colMeans, parallel = T) {

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
  mat <- mat[, new.cols]
  mat
}

# plotting --------------------------------------------------------------------


lm_corr_eqn <- function(df, method = "pearson", display.r2 = FALSE, display.pval = FALSE){
    m <- lm(y ~ x, df);
    ct <- cor.test(df$x, df$y, method = method)
    estimate <- as.numeric(ct$estimate)
    if(display.r2 == TRUE) { estimate <- estimate*estimate }
    pval <- ct$p.value
    cat(paste0("method = ", method, " estimate = ", estimate, " pval = ", pval, "\n"))
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

plot.all.cell.type.correlations <- function(data, title) {
  g <- ggplot(data, aes(x = actual, y = prediction))
  g <- g + ggtitle(title)
  g <- g + geom_point()
  g <- g + geom_smooth(method = "lm")
  g <- g + facet_wrap(~ cell.type, scales = "free")
  g
}

plot.individual.cell.type.correlations <- function(data, title) {
  gs <-
    dlply(data, .variables = c("cell.type"),
          .fun = function(df) {
	           g <- plot.correlation(df$actual, df$prediction, method = "spearman", display.pval = TRUE)
		   g <- g + xlab("Actual") + ylab("Prediction")
                   g <- g + ggtitle(paste0(title, df$cell.type[1]))
		   return(g)
                   g <- ggplot(df, aes(x = actual, y = prediction))
                   g <- g + geom_point()
                   g <- g + geom_smooth(method = "lm")
                   g
		 })
  gs		   
}



# misc ------------------------------------------------------------------------

require(magrittr)
require(preprocessCore)

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