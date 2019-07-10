run.deconvolution.methods <- function(metadata.synId, output.folder.synId) {

  metadata <- download.and.parse.dataset.metadata(metadata.synId)
  dataset <- get.dataset.name(metadata)

  hugo.expr.mat <- download.hugo.expression.matrix(metadata)

  ## Create an expression matrix whose first column is the HUGO gene name,
  ## as required by CIBERSORT.
##  hugo.expr.gene.col.mat <- cbind(HUGO = rownames(hugo.expr.mat), hugo.expr.mat)
  hugo.expr.gene.col.mat <- hugo.expr.mat %>%
    dplyr::rename(HUGO = Gene)

  cibersort.expr.input.file <- paste0(dataset, "-hugo-expr-cibersort-input.tsv")
  write.table(file = cibersort.expr.input.file, hugo.expr.gene.col.mat, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

  gt.mat.coarse <- download.coarse.grained.ground.truth(metadata)
  gt.mat.fine <- download.fine.grained.ground.truth(metadata)  

  gt.mat.coarse <- gt.mat.coarse %>%
    select(-dataset.name)

  gt.mat.fine <- gt.mat.fine %>%
    select(-dataset.name)

  measured.col <- "measured"
  predicted.col <- "prediction"
##  gt.mat.coarse <- melt(as.matrix(gt.mat.coarse))
##  colnames(gt.mat.coarse) <- c("sample.id", "cell.type", measured.col)

##  gt.mat.fine <- melt(as.matrix(gt.mat.fine))
##  colnames(gt.mat.fine) <- c("sample.id", "cell.type", measured.col)

  cs.preds <-
    run.cibersort.models(list(cibersort.expr.input.file), list(dataset))
  coarse.cs.pred <- cs.preds[["coarse"]]
  fine.cs.pred <- cs.preds[["fine"]]  

  hugo.expr.mat <- hugo.expr.mat %>%
    column_to_rownames(var = "Gene")
    
  coarse.mcp.pred <-
    run.coarse.grained.mcpcounter.model(list(hugo.expr.mat),
                                        list(dataset)) %>% as.data.frame()

  coarse.mcp <- merge(coarse.mcp.pred, gt.mat.coarse)
  coarse.cs <- merge(coarse.cs.pred, gt.mat.coarse)
  fine.cs <- merge(fine.cs.pred, gt.mat.fine)    

  ## Plot MCPcounter-based coarse-grained predictions
  g.coarse.mcp.all <-
    plot.all.cell.type.correlations(coarse.mcp, "MCP Counter", x.col = measured.col, y.col = predicted.col)

  g.coarse.mcp.individual <-
    plot.individual.cell.type.correlations(coarse.mcp, "MCP Counter: ", x.col = measured.col, y.col = predicted.col)

  l <- c(list(g.coarse.mcp.all), g.coarse.mcp.individual)
  file <- paste0(dataset, "-mcp-counter-coarse-correlations.pdf")
  pdf(file, onefile = TRUE)
  l_ply(l, .fun = function(g) print(g))
  d <- dev.off()
  
  f <- File(file, parentId = output.folder.synId, synapseStore = TRUE)
  ss <- synStore(f, executed = script_url, forceVersion = FALSE)
  synId <- get.synapse.id(ss)

  ## Plot CIBERSORT-based coarse-grained predictions
  g.coarse.cs.all <-
    plot.all.cell.type.correlations(coarse.cs, "CIBERSORT", x.col = measured.col, y.col = predicted.col)

  g.coarse.cs.individual <-
    plot.individual.cell.type.correlations(coarse.cs, "CIBERSORT: ", x.col = measured.col, y.col = predicted.col)

  l <- c(list(g.coarse.cs.all), g.coarse.cs.individual)
  file <- paste0(dataset, "-cibersort-coarse-correlations.pdf")
  pdf(file, onefile = TRUE)
  l_ply(l, .fun = function(g) print(g))
  d <- dev.off()
  
  f <- File(file, parentId = output.folder.synId, synapseStore = TRUE)
  ss <- synStore(f, executed = script_url, forceVersion = FALSE)
  synId <- get.synapse.id(ss)

  ## Plot CIBERSORT-based fine-grained predictions
  g.fine.cs.all <-
    plot.all.cell.type.correlations(fine.cs, "CIBERSORT", x.col = measured.col, y.col = predicted.col)

  g.fine.cs.individual <-
    plot.individual.cell.type.correlations(fine.cs, "CIBERSORT: ", x.col = measured.col, y.col = predicted.col)

  l <- c(list(g.fine.cs.all), g.fine.cs.individual)
  file <- paste0(dataset, "-cibersort-fine-correlations.pdf")
  pdf(file, onefile = TRUE)
  l_ply(l, .fun = function(g) print(g))
  d <- dev.off()
  
  f <- File(file, parentId = output.folder.synId, synapseStore = TRUE)
  ss <- synStore(f, executed = script_url, forceVersion = FALSE)
  synId <- get.synapse.id(ss)


}
