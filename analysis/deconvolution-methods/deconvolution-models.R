run.deconvolution.methods <- function(metadata.synId, output.folder.synId) {

  metadata <- download.and.parse.dataset.metadata(metadata.synId)
  dataset <- get.dataset.name(metadata)

  hugo.expr.mat <- download.hugo.expression.matrix(metadata)

  gt.mat.coarse <- download.coarse.grained.ground.truth(metadata)
  gt.mat.fine <- download.fine.grained.ground.truth(metadata)  

  samples <- unique(c(gt.mat.coarse$sample, gt.mat.fine$sample))
  flag <- colnames(hugo.expr.mat) %in% c(samples, "HUGO", "Hugo", "Gene")
  hugo.expr.mat <- hugo.expr.mat[, flag]
  
  ## Create an expression matrix whose first column is the HUGO gene name,
  ## as required by CIBERSORT.
  use.matrix.with.hugo.col <- FALSE
##  hugo.expr.gene.col.mat <- cbind(HUGO = rownames(hugo.expr.mat), hugo.expr.mat)
  hugo.expr.gene.col.mat <- NULL
  if(use.matrix.with.hugo.col) {
    hugo.expr.gene.col.mat <- hugo.expr.mat %>% dplyr::rename(HUGO = Hugo)
  } else { 
    hugo.expr.gene.col.mat <- hugo.expr.mat %>% dplyr::rename(HUGO = Gene)
  }

  cibersort.expr.input.file <- paste0(dataset, "-hugo-expr-cibersort-input.tsv")
  write.table(file = cibersort.expr.input.file, hugo.expr.gene.col.mat, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

  measured.col <- "measured"
  predicted.col <- "prediction"

  cs.preds <-
    run.cibersort.models(list(cibersort.expr.input.file), list(dataset))
  coarse.cs.pred <- cs.preds[["coarse"]]
  fine.cs.pred <- cs.preds[["fine"]]  

  if(use.matrix.with.hugo.col) {
    hugo.expr.mat <- hugo.expr.mat %>% column_to_rownames(var = "Hugo") 
  } else {
    hugo.expr.mat <- hugo.expr.mat %>% column_to_rownames(var = "Gene")
  }
  coarse.mcp.pred <-
    run.coarse.grained.mcpcounter.model(list(hugo.expr.mat),
                                        list(dataset)) %>% as.data.frame()

  methods <- c("spearman", "pearson")
  
  if(!((class(gt.mat.coarse) == "logical") && is.na(gt.mat.coarse))) {
    gt.mat.coarse <- gt.mat.coarse %>%
      dplyr::select(-dataset.name)
    coarse.mcp <- merge(coarse.mcp.pred, gt.mat.coarse)
    coarse.cs <- merge(coarse.cs.pred, gt.mat.coarse)

    ## Plot MCPcounter-based coarse-grained predictions
    for(method in methods) {
      if(nrow(coarse.mcp) > 0) {
        g.coarse.mcp.all <-
          plot.all.cell.type.correlations(coarse.mcp, paste0("MCP Counter ", method), x.col = measured.col, y.col = predicted.col, method = method)
    
        g.coarse.mcp.individual <-
          plot.individual.cell.type.correlations(coarse.mcp, paste0("MCP Counter ", method, ": "), x.col = measured.col, y.col = predicted.col, method = method)
    
        l <- c(list(g.coarse.mcp.all), g.coarse.mcp.individual)
        file <- paste0(dataset, "-mcp-counter-coarse-", method, "-correlations.pdf")
        pdf(file, onefile = TRUE)
        l_ply(l, .fun = function(g) print(g))
        d <- dev.off()
      
        f <- File(file, parentId = output.folder.synId, synapseStore = TRUE)
        ss <- synStore(f, executed = script_url, forceVersion = FALSE)
        synId <- get.synapse.id(ss)
      }
  
      if(nrow(coarse.cs) > 0) {
        ## Plot CIBERSORT-based coarse-grained predictions
        g.coarse.cs.all <-
          plot.all.cell.type.correlations(coarse.cs, paste0("CIBERSORT ", method), x.col = measured.col, y.col = predicted.col, method = method)
  
        g.coarse.cs.individual <-
          plot.individual.cell.type.correlations(coarse.cs, paste0("CIBERSORT ", method, ": "), x.col = measured.col, y.col = predicted.col, method = method)
  
        l <- c(list(g.coarse.cs.all), g.coarse.cs.individual)
        file <- paste0(dataset, "-cibersort-coarse-", method, "-correlations.pdf")
        pdf(file, onefile = TRUE)
        l_ply(l, .fun = function(g) print(g))
        d <- dev.off()
      
        f <- File(file, parentId = output.folder.synId, synapseStore = TRUE)
        ss <- synStore(f, executed = script_url, forceVersion = FALSE)
        synId <- get.synapse.id(ss)
      }
    } ## for(method in methods)
  }

  if(!((class(gt.mat.fine) == "logical") && is.na(gt.mat.fine))) {
    gt.mat.fine <- gt.mat.fine %>%
      dplyr::select(-dataset.name)

    fine.cs <- merge(fine.cs.pred, gt.mat.fine)    

    for(method in methods) {

      if(nrow(fine.cs) > 0) {
        ## Plot CIBERSORT-based fine-grained predictions
        g.fine.cs.all <-
          plot.all.cell.type.correlations(fine.cs, paste0("CIBERSORT ", method), x.col = measured.col, y.col = predicted.col, method = method)
    
        g.fine.cs.individual <-
          plot.individual.cell.type.correlations(fine.cs, paste0("CIBERSORT ", method, ": "), x.col = measured.col, y.col = predicted.col, method = method)
    
        l <- c(list(g.fine.cs.all), g.fine.cs.individual)
        file <- paste0(dataset, "-cibersort-fine-", method, "-correlations.pdf")
        pdf(file, onefile = TRUE)
        l_ply(l, .fun = function(g) print(g))
        d <- dev.off()
      
        f <- File(file, parentId = output.folder.synId, synapseStore = TRUE)
        ss <- synStore(f, executed = script_url, forceVersion = FALSE)
        synId <- get.synapse.id(ss)
      }
    } ## for(method in methods)
  }

}
