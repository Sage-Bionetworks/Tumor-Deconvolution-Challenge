translate.matrix <- function(mat, translation.tbl, from.col, to.col) {
  # rownames(mat) <- mat$Gene
  # translation.tbl <- translation.tbl[translation.tbl[,from.col] %in% rownames(mat),]
  # mat <- mat[, !(colnames(mat) == "Gene")]
  # tmp <- t(mat)
  ttbl <- translation.tbl[translation.tbl[, from.col] %in% rownames(mat),]
  tmp <- (t(mat))[, unique(ttbl[, from.col])]
  ret <- aggregate_matrix(tmp, ttbl, to.col, from.col, aggregate_fun = sum)
  #ret <- cbind(Gene = colnames(ret), as.data.frame(t(ret)))
  return(t(ret))
}

# Based on https://github.com/hbc/knowledgebase/blob/master/scrnaseq/pseudobulkDE_edgeR.md
# Metadata should have cell ids as rownames and columns patient.col and cluster.col.
create.pseudobulk <- function(expr.mat, metadata, patient.col = "Patient", cluster.col = "celltype_coarse") {
  common.cells <- intersect(colnames(expr.mat), rownames(metadata))
  expr.mat <- expr.mat[, common.cells]
  metadata <- metadata[common.cells,]
  groups <- metadata[, c(patient.col, cluster.col)]
  cell.type.freq <- melt(table(groups))
  colnames(cell.type.freq)[3] <- "cnt"
  cell.type.freq <- ddply(cell.type.freq, .variables = c("Patient"), .fun = function(df) { tot = sum(df$cnt); df$frac <- df$cnt / tot; df})
  pb <- aggregate.Matrix(t(expr.mat), 
                         groupings = groups, fun = "sum") 
  lst <- list("freq" = cell.type.freq, "pb" = pb)
  return(lst)
}

## # expr.mat is assumed to be a pseudobulk matrix, with first column Gene and with other columns named as in patient.cluster.col 
## # in metadata (e.g., <patient id>_<cluster id>)
## # metadata should also have columns patient.col and cluster.col
create.admixtures.from.pseudobulk <- function(expr.mat, metadata, patient.cluster.col = "patient_cluster",
                                              patient.col = "Patient", cluster.col = "celltype_coarse") {
  common.cell.types <- intersect(colnames(expr.mat), metadata[, patient.cluster.col])
  expr.mat <- expr.mat[, common.cell.types]
  metadata <- metadata[metadata[, patient.cluster.col] %in% common.cell.types, ]
  groups <- metadata[, c(patient.col, cluster.col)]
  cell.type.freq <- melt(table(groups))
  colnames(cell.type.freq)[3] <- "cnt"
  cell.type.freq <- ddply(cell.type.freq, .variables = c("Patient"), .fun = function(df) { tot = sum(df$cnt); df$frac <- df$cnt / tot; df})
  pb <- aggregate.Matrix(t(expr.mat), 
                         groupings = groups, fun = "sum") 
  lst <- list("freq" = cell.type.freq, "pb" = pb)
  return(lst)
}

fine_cell_types <- c(
  "myeloid.dendritic.cells",
  "endothelial.cells",
  "fibroblasts",
  "macrophages",
  "memory.CD4.T.cells",
  "memory.CD8.T.cells",
  "monocytes",
  "naive.B.cells",
  "naive.CD4.T.cells",
  "naive.CD8.T.cells",
  "neutrophils",
  "NK.cells",
  "regulatory.T.cells",
  "memory.B.cells"
)

coarse_cell_types <- c(
  "B.cells",
  "CD4.T.cells",
  "CD8.T.cells",
  "NK.cells",
  "neutrophils",
  "monocytic.lineage",
  "fibroblasts",
  "endothelial.cells"
)

fine.to.coarse.trans <- list(
  "myeloid.dendritic.cells" = "monocytic.lineage",
  "endothelial.cells" = "endothelial.cells",
  "fibroblasts" = "fibroblasts",
  "macrophages" = "monocytic.lineage",
  "memory.CD4.T.cells" = "CD4.T.cells",
  "memory.CD8.T.cells" = "CD8.T.cells",
  "monocytes" = "monocytic.lineage",
  "naive.B.cells" = "B.cells",
  "naive.CD4.T.cells" = "CD4.T.cells",
  "naive.CD8.T.cells" = "CD8.T.cells",
  "fibroblasts" = "fibroblasts",
  "NK.cells" = "NK.cells",
  "regulatory.T.cells" = "CD4.T.cells",
  "memory.B.cells" = "B.cells"
)

fine.to.coarse.trans.tbl <- 
  data.frame(celltype_fine = names(fine.to.coarse.trans), celltype_coarse = unlist(unname(fine.to.coarse.trans)))

# Find cell types that cluster "incorrectly" -- i.e., strongly with other cell types
# in umap space 
# df should have samples as rownames and columns umap1, umap2, and cell.type
find.inconsistent.umap.cell.types <- function(df) {
  dst <- as.matrix(dist(as.matrix(df[,c("umap1", "umap2")])))
  flag <- unlist(lapply(1:nrow(dst), function(i) {
    dists <- dst[i,]
    dists <- dists[-i]
    lbl <- df$cell.type[i]
    lbls <- df$cell.type[-i]
    full.lbls <- rownames(df)[-i]
    o <- order(dists, decreasing=FALSE)
    lbls <- lbls[o]
    # Don't distinguish between CD4 and CD8 t cells -- they are likely to be intermixed
    lbls <- gsub(lbls, pattern="CD..", replacement="")
    lbl <- gsub(lbl, pattern="CD..", replacement="")
    dists <- dists[o]
    full.lbls <- full.lbls[o]
    n.top <- 10
    tp <- lbls[1:n.top]
    tp <- as.data.frame(table(tp))
    max.ct <- as.character(tp[which.max(tp$Freq), "tp"])
    if(!(all(max.ct == lbl))) {
      print(rownames(dst)[i])
      print(tp)
      print(full.lbls[1:n.top])
      cat("\n")
      return(TRUE)
    }
    return(FALSE)
  }))
  return(df[flag,])
}

plot.challenge.umap <- function(df, plot.labels = TRUE) {
  cbbPalette <- Blue2OrangeRed14Steps
  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  pal <- cbbPalette
  if(length(unique(df$cell.type)) > length(pal)) {
    pal <- Blue2OrangeRed14Steps
  }
  
  # g <- ggplot(data = df, aes(x = umap1, y = umap2, colour = cell.type)) + geom_point() + scale_colour_manual(values = pal)
  df$cell.type <- unlist(lapply(df$cell.type, function(str) substr(str, start=0, stop=4)))
  if(plot.labels) {
    g <- ggplot(data = df, aes(x = umap1, y = umap2, colour = cell.type)) + geom_point() + scale_colour_manual(values = pal) + geom_text(aes(label = cell.type))
  } else {
    g <- ggplot(data = df, aes(x = umap1, y = umap2, colour = cell.type)) + geom_point() + scale_colour_manual(values = pal)   
  }
  g
}

# Create umap using deconv.genes of challenge samples
challenge.umap <- function(mat) {
  #source("/projects/compsci/jgeorge/whitebr/Tumor-Deconvolution-Challenge/analysis/utils.R") # for get.deconvolution.genes
  #install.packages("remotes")
  #remotes::install_github("omnideconv/immunedeconv")
  #deconv.genes <- get.deconvolution.genes()
  deconv.gene.file <- "/projects/compsci/jgeorge/whitebr/Tumor-Deconvolution-Challenge/external/deconv-genes.tsv"
  deconv.genes <- read.table(deconv.gene.file, header=TRUE)$gene
  
  # Assume colnames are of the form <patient>_<challenge cell type>
  # Limit to challenge cell types
  splitf <- sapply(stringr::str_split(colnames(mat), 
                                      pattern = "_",n = 2), `[`, 2)
  flag <- splitf %in% c(coarse_cell_types, fine_cell_types)
  mat <- mat[, flag]
  
  # Create umap using deconvolution genes
  umap <- umap(t(mat[rownames(mat) %in% deconv.genes,]))
  df <- data.frame(umap1 = umap$layout[,1], umap2 = umap$layout[,2], sample = rownames(umap$layout))
  splitf <- sapply(stringr::str_split(df$sample,pattern = "_",n = 2), `[`, 2)
  df$cell.type <- splitf
  splitf <- sapply(stringr::str_split(df$sample,pattern = "_",n = 2), `[`, 1)
  df$patient <- splitf
  df
}

run.da505_ <- function(mat, grain="fine") {
  path <- "/projects/compsci/jgeorge/USERS/whitebr/Dream_Deconv_Challenge_Team_DA505/"
  cur.wd <- getwd()
  setwd(path)
  source("cps_functions.R")
  suppressPackageStartupMessages(p_load(glmnet))
  suppressPackageStartupMessages(p_load(e1071))
  rownames(mat) <- mat$Gene
  mat <- mat[, -1]
  if(grain=="fine") {
    res <- do_CPS_fine_glmnet(mat)
  } else {
    res <- do_CPS_coarse_glmnet(mat)
  }
  setwd(cur.wd)
  
  res <- melt(as.matrix(res))
  colnames(res) <- c("cell.type", "sample", "value")
  
  return(res)
  mat.file <- tempfile()
  out.file <- tempfile()
  write.table(mat, file=mat.file, row.names=FALSE, col.names=TRUE, sep=",", quote=FALSE)
  cmd <- paste0("Rscript ", path, "./cps_run.R --input=", mat.file, " --output=", out.file, " --", grain)
  system(cmd)
  ret <- read.table(out.file)
  file.remove(mat.file)
  file.remove(out.file)
  setwd(cur.wd)
  ret
}

run.da505.fine <- function(mat) {
  run.da505_(mat, grain="fine")
}

run.da505.coarse <- function(mat) {
  run.da505_(mat, grain="coarse")
}

# tpm.expr is expected to have a Gene column
# tpm.expr is expected to have a Gene column
run.mcpcounter.coarse <- function(tpm.expr) {
  rownames(tpm.expr) <- tpm.expr$Gene
  tpm.expr <- tpm.expr[, !(colnames(tpm.expr) == "Gene")]
  res <- MCPcounter.estimate(tpm.expr, featuresType="HUGO_symbols")
  res <- melt(res)
  colnames(res) <- c("mcpcounter.cell.type", "sample", "value")
  
  translation_df <- tibble::tribble(
    ~cell.type, ~mcpcounter.cell.type,
    "B.cells", "B lineage",
    "CD8.T.cells", "CD8 T cells",
    "NK.cells", "NK cells",
    "neutrophils", "Neutrophils",
    "monocytic.lineage", "Monocytic lineage",
    "fibroblasts", "Fibroblasts",
    "endothelial.cells", "Endothelial cells")
  
  res <- merge(res, translation_df)
  # res <- res[, c("cell.type", "sample", "value")]
  res
}

run.xcell_ <- function(tpm.expr, translation_df) {
  rownames(tpm.expr) <- tpm.expr$Gene
  tpm.expr <- tpm.expr[, !(colnames(tpm.expr) == "Gene")]
  
  res <- xCell::xCellAnalysis(tpm.expr,
                              rnaseq = TRUE,
                              cell.types.use = translation_df$xcell.cell.type)
  res <- melt(res)
  colnames(res) <- c("xcell.cell.type", "sample", "value")
  
  res <- merge(res, translation_df)
  # res <- res[, c("cell.type", "sample", "value")]
  res
}


run.xcell.fine <- function(tpm.expr) {
  
  # Translate xCell output to challenge subtypes
  translation_df <- tibble::tribble(
    ~cell.type, ~xcell.cell.type,
    "memory.B.cells", "Memory B-cells",
    "naive.B.cells", "naive B-cells",
    "memory.CD4.T.cells", "CD4+ memory T-cells",
    "naive.CD4.T.cells", "CD4+ naive T-cells",
    "regulatory.T.cells", "Tregs",
    "memory.CD8.T.cells", "CD8+ Tem",
    "naive.CD8.T.cells", "CD8+ naive T-cells",
    "NK.cells", "NK cells",
    "neutrophils", "Neutrophils",
    "monocytes", "Monocytes",
    "myeloid.dendritic.cells", "DC",
    "macrophages", "Macrophages",
    "fibroblasts", "Fibroblasts",
    "endothelial.cells", "Endothelial cells"
  )
  
  run.xcell_(tpm.expr, translation_df)
  
}

run.xcell.coarse <- function(tpm.expr) {
  
  translation_df <- tibble::tribble(
    ~cell.type, ~xcell.cell.type,
    "B.cells", "B-cells",
    "CD4.T.cells", "CD4+ T-cells",
    "CD8.T.cells", "CD8+ T-cells",
    "NK.cells", "NK cells",
    "neutrophils", "Neutrophils",
    "monocytic.lineage", "Monocytes",
    "fibroblasts", "Fibroblasts",
    "endothelial.cells","Endothelial cells"
  )
  
  run.xcell_(tpm.expr, translation_df)
}

plot.deconv.results <- function(res, gt.col = "frac") {
  # g <- ggplot(data = res, aes(x = frac, y = value)) + geom_point() + facet_wrap(~ cell.type, scales="free") + xlab("Ground Truth Fraction") + ylab("Prediction")
  g <- ggscatter(res, x = gt.col, y = "value",
                 add = "reg.line",  # Add regressin line
                 add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                 conf.int = TRUE # Add confidence interval
  )
  # Add correlation coefficient
  g <- facet(g, facet.by = "cell.type", scales="free") 
  g <- g + stat_cor(method = "pearson")
  g <- g + xlab("Ground Truth Fraction") + ylab("Prediction") + theme(plot.title = element_text(hjust = 0.5))
  g
}

format.cibersort.results <- function(res, translation_df) {
  res <- melt(res)
  colnames(res) <- c("sample", "cibersort.cell.type", "value")
  
  res <- merge(res, translation_df)
  res <- ddply(res, .variables = c("sample", "cell.type"),
               .fun = function(df) {
                 data.frame(cell.type = df[1, "cell.type"], sample = df[1, "sample"], value = sum(df$value))
               })
  # res <- res[, c("cell.type", "sample", "value")]
  res
}

cibersort_fine_translation_df <- tibble::tribble(
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
  "endothelial.cells", "CD31",
  "fibroblasts", "CD10",
  "cancer", "EPCAM"
)


cibersort_coarse_translation_df <- tibble::tribble(
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
  "endothelial.cells", "CD31",
  "fibroblasts", "CD10",
  "cancer", "EPCAM"
)

run.cibersort_ <- function(tpm.expr, sig.file, file, translation_df = NULL) {
  write.table(file = file, tpm.expr, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
  res <- CIBERSORT(sig.file, file, QN = FALSE, absmean = TRUE, abs_method = 'sig.score')
  
  if(!is.null(translation_df)) {
    res <- melt(res)
    colnames(res) <- c("sample", "cibersort.cell.type", "value")
    
    res <- merge(res, translation_df)
    res <- ddply(res, .variables = c("sample", "cell.type"),
                 .fun = function(df) {
                   data.frame(cell.type = df[1, "cell.type"], sample = df[1, "sample"], value = sum(df$value))
                 })
    # res <- res[, c("cell.type", "sample", "value")]
  }
  file.remove(file)
  res
}

run.cibersort.fine <- function(tpm.expr, sig.file, file) {
  
  # Translate cibersort output to challenge subtypes
  
  run.cibersort_(tpm.expr, sig.file, file, cibersort_fine_translation_df)
  
}

run.cibersort.coarse <- function(tpm.expr, sig.file, file) {
  
  # Translate cibersort output to challenge subtypes
  
  run.cibersort_(tpm.expr, sig.file, file, cibersort_coarse_translation_df)
  
}

run.all.deconv <- function(ensg.mat, sym.mat, cell.type.freq, grain="fine") {
  cs.res <- NULL
  
  
  cs.lm22.out.file <- tempfile(pattern = "cs-lm22")
  cs.2l.out.file <- tempfile(pattern = "csx-2l")
  
  cs.lm22.res <- run.cibersort_(sym.mat, cibersort.lm22.sig.file, cs.lm22.out.file)
  cs.2l.res <- run.cibersort_(sym.mat, cibersort.2l.sig.file, cs.2l.out.file)
  cs.lm22.res <- cbind("Mixture" = rownames(cs.lm22.res), as.data.frame(cs.lm22.res))
  cs.2l.res <- cbind("Mixture" = rownames(cs.2l.res), as.data.frame(cs.2l.res))
  
  csx.res <- as.data.frame(merge(cs.2l.res, cs.lm22.res, by=c("Mixture")))
  cell.types <- colnames(cs.lm22.res)
  cell.types <- cell.types[!(cell.types %in% c("Mixture", "P-value", "Correlation", "RMSE"))]
  # Scale by immune proportion
  for(ct in cell.types) { csx.res[, ct] <- csx.res[,ct] * csx.res[,"CD45"]}
  csx.res <- csx.res[,c("Mixture", cell.types, "CD10", "CD31", "EPCAM")]
  rownames(csx.res) <- csx.res[,1]
  csx.res <- csx.res[,-1]
  
  cs.res <- cs.lm22.res
  cs.res <- cs.res[, c("Mixture", cell.types)]
  rownames(cs.res) <- cs.res[,1]
  cs.res <- cs.res[,-1]
  if(grain == "fine") {
    csx.res <- format.cibersort.results(as.matrix(csx.res), cibersort_fine_translation_df)
    csx.res <- merge(csx.res, cell.type.freq, by.x = c("cell.type", "sample"), by.y = c("cell.type", "sample.id"))
    
    cs.res <- format.cibersort.results(as.matrix(cs.res), cibersort_fine_translation_df)
    cs.res <- merge(cs.res, cell.type.freq, by.x = c("cell.type", "sample"), by.y = c("cell.type", "sample.id"))
  } else {
    csx.res <- format.cibersort.results(as.matrix(csx.res), cibersort_coarse_translation_df)
    csx.res <- merge(csx.res, cell.type.freq, by.x = c("cell.type", "sample"), by.y = c("cell.type", "sample.id"))
    
    cs.res <- format.cibersort.results(as.matrix(cs.res), cibersort_coarse_translation_df)
    cs.res <- merge(cs.res, cell.type.freq, by.x = c("cell.type", "sample"), by.y = c("cell.type", "sample.id"))
  }
  
  da.res <- run.da505_(ensg.mat, grain=grain)
  da.res <- merge(da.res, cell.type.freq, by.x = c("cell.type", "sample"), by.y = c("cell.type", "sample.id"))
  
  mcp.res <- NULL
  if(grain == "coarse") {
    mcp.res <- run.mcpcounter.coarse(sym.mat)
    mcp.res <- merge(mcp.res, cell.type.freq, by.x = c("cell.type", "sample"), by.y = c("cell.type", "sample.id"))
  }
  
  xcell.res <- NULL
  if(grain == "coarse") {
    xcell.res <- run.xcell.coarse(sym.mat)
  } else {
    xcell.res <- run.xcell.fine(sym.mat)
  }
  xcell.res <- merge(xcell.res, cell.type.freq, by.x = c("cell.type", "sample"), by.y = c("cell.type", "sample.id"))
  
  res <- list("da505" = da.res, "mcp" = mcp.res, "xcell" = xcell.res, "cs" = cs.res, "csx" = csx.res)
  
}

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

