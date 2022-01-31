coarse.cell.types <-
  c("B.cells", "CD4.T.cells", "CD8.T.cells", "NK.cells", "neutrophils", "monocytic.lineage",
    "fibroblasts", "endothelial.cells")

## See https://stackoverflow.com/questions/27803710/ggplot2-divide-legend-into-two-columns-each-with-its-own-title
plot.anno.heatmap.with.multiple.legends <-
    function(df, id.col, anno.columns, anno.pals) {

        suppressPackageStartupMessages(p_load("RColorBrewer"))
        df <- df[, c(id.col, anno.columns)]

        ## Assume annotations are characters
        ## NB: id.col is a factor
        for(col in c(anno.columns)) {
            df[, col] <- as.character(df[, col])
        }

        columns <- 1:length(anno.columns)
        names(columns) <- anno.columns

        color.vecs <-
            llply(columns,
                  .fun = function(idx) {
                      anno.col <- anno.columns[idx]
                      vec <- unique(df[, anno.col])
                      len <- length(vec)
                      colors <- brewer.pal(len, anno.pals[idx])
                      names(colors) <- vec
                      colors
                  })

        all.colors <- Reduce("c", color.vecs)
        names(all.colors) <- Reduce("c", unlist(lapply(color.vecs, names)))

        names(anno.columns) <- anno.columns
        anno.df <- ldply(anno.columns,
                     .fun = function(anno.col) {
                         data.frame(val = df[, anno.col], id = df[, id.col])
                     })
        colnames(anno.df)[1] <- "type"

        full.plot <-
            ggplot(anno.df, aes(y = id, x = type, fill = val)) + geom_tile() +
            scale_fill_manual(values = all.colors) +
            theme(legend.position="none")

        full.plot <- full.plot + theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
                                       axis.ticks.y = element_blank(), text = element_text(size = 18),
                                       axis.text.x = element_text(angle = 45, hjust = 1),
                                       axis.title.x = element_blank())
        
        
        legends <-
            llply(anno.columns, .parallel = FALSE,
                  .fun = function(anno.col) {
                      flag <- anno.df$type == anno.col
                      leg.colors <- all.colors[names(all.colors) %in% anno.df[flag,"val"]]
                      g <- ggplot(anno.df[flag, ], aes_string(x = "id", y = "type", fill = "val"))
                      g <- g + geom_tile()
                      g <- g + scale_fill_manual(values = leg.colors, name = anno.col)
                  })
        return(list("full.plot" = full.plot, "legends" = legends))
    }

get.method.annotations <- function() {
  synId <- "syn23395242"
  obj <- synGet(synId, downloadFile = TRUE)
  method.anno <- read.xlsx(obj$path, sheetIndex = 1)
  
  method.anno$method.type <- as.character(method.anno$method.type)
  method.anno$output.type <- as.character(method.anno$output.type)
  method.anno[, round.col] <- as.character(method.anno[, round.col])
  flag <- method.anno[, round.col] == "NA"
  method.anno[flag, round.col] <- NA
  method.anno[, subchallenge.col] <- as.character(method.anno[, subchallenge.col])
  flag <- method.anno[, subchallenge.col] == "NA"
  method.anno[flag, subchallenge.col] <- NA

  method.rename.list <-
      list("NNLS" = "NNLS",
           "summary" = "SUM",
           "other" = "OTH",
           "other regression" = "REG",
           "unknown" = "UNK",
           "SVR" = "SVR",
           "DNN" = "DNN",
           "ensemble" = "ENS",
           "NMF" = "NMF",
           "probabilistic inference" = "PI")
  method.rename.df <- data.frame(method.type = names(method.rename.list), Method = as.character(method.rename.list))

  output.rename.list <-
      list("fraction" = "Frac",
           "proportion" = "Prop",
           "normalized.score" = "Norm",
           "score" = "Score")
  output.rename.df <- data.frame(output.type = names(output.rename.list), Output = as.character(output.rename.list))

  method.anno <- merge(method.anno, method.rename.df)
  method.anno <- merge(method.anno, output.rename.df)
  method.anno$Output <- as.character(method.anno$Output)
  method.anno$Method <- as.character(method.anno$Method)
  method.anno
}


get.round.specific.annotations <- function(method.anno, round) {
  method.anno.round <-
    ddply(method.anno,
          .variables = c(method.name.col, subchallenge.col),
	  .fun = function(df) {
                   flag <- is.na(df[, round.col]) | (df[, round.col] == round)
		   if(any(flag)) {
		     ret <- df[flag,,drop=F]
		     if(nrow(ret) != 1) { print(ret); print(df); stop(paste0("Wrong number of rows for ",
		                                                  df[1,method.name.col], " in ", df[1, subchallenge.col], "!\n")); }
		     return(ret)
		   }
                   flag <- df[, round.col] == "latest"
		   ret <- df[flag,,drop=F]
		   if(nrow(ret) != 1) { print(ret); print(df); stop(paste0("Wrong number of rows2 for ",
		                                                df[1,method.name.col], " in ", df[1, subchallenge.col], "!\n")); }
		   return(ret)
                 })
  method.anno.round
}

get.cell.type.translation <- function() {
    list("naive.B.cells" = "naive B",
         "memory.B.cells" = "memory B",
         "fibroblasts" = "fibroblasts",
         "neutrophils" = "neutrophils",
         "endothelial.cells" = "endothelial",
         "monocytes" = "monocytes",
         "NK.cells" = "NK",
         "macrophages" = "macrophages",
         "memory.CD8.T.cells" = "memory CD8 T",
         "regulatory.T.cells" = "Tregs",
         "B.cells" = "B",
         "naive.CD8.T.cells" = "naive CD8 T",
         "monocytic.lineage" = "monocytic lineage",
         "naive.CD4.T.cells" = "naive CD4 T",
         "myeloid.dendritic.cells" = "myeloid DCs",
         "CD8.T.cells" = "CD8 T",
         "memory.CD4.T.cells" = "memory CD4 T",
         "CD4.T.cells" = "CD4 T")         
}

rename.cell.types <- function(df, from.col = "cell.type", to.col = "cell.type") {
    cell.type.trans <- get.cell.type.translation()
    for(nm in names(cell.type.trans)) {
        flag <- (grepl(df[, from.col], pattern = nm))
        df[flag, to.col] <- cell.type.trans[[nm]]
    }
    df
}

safe.merge <- function(x, y, ...) {

    orig.nrow <- nrow(x)
    x <- merge(x, y, ...)
    new.nrow <- nrow(x)
    if(orig.nrow != new.nrow) {
        stop("Changed rows\n")
    }
    x
}

## Assign team names and rounds to results
assign.result.team.names.and.rounds <- function(res, error.fun = stop) {
    ## Read in the final predictions from the competitive phase (i.e., when coarse- and fine-grained
    ## datasets differed)
    synId <- "syn22149603"
    obj <- synGet(synId, downloadFile=TRUE)
    res.comp <- read.table(obj$path, sep=",", header=TRUE, as.is=TRUE, stringsAsFactors=FALSE)

    ## Use the competitive phase results to translate objectIds to teams / submitters
    ## The final name in the "path" of the repo_name for the post-competitive phase is the
    ## objectId in the competitive phase (except for the baselines)
    cat("Assigning objectId\n")
    map <- data.frame(repo_name = as.character(unique(res$repo_name)), stringsAsFactors = FALSE)

    map$objectId <-
        unlist(lapply(map$repo_name,
                      function(str) {
                          strs <- unlist(strsplit(str, "/"))
                          strs[length(strs)]
                      }))

    cat("Assigning comparator\n")
    map$comparator <-
        unlist(lapply(as.character(map$objectId),
                      function(str) {
                          comps <- get.comparators()
                          for(comp in comps) {
                              if(grepl(str, pattern=comp, ignore.case=TRUE)) { return(TRUE) }
                          }
                          return(FALSE)
                      }))
    
    res <- res[, !(colnames(res) %in% c("objectId", "comparator"))]
    res <- safe.merge(res, map, by = c("repo_name"))

    cat("Ensuring\n")
    
    ## Ensure that all objectIds match between competitive and post-competitive
    flag <- grepl(res.comp$repo_name, pattern="baseline") | (res.comp$objectId %in% res$objectId)
    if(!all(flag)) {
        print(table(flag))
        print(unique(res.comp[!flag, c("repo_name", "objectId", "submitterId")]))
        l_ply(unique(res.comp[!flag, "submitterId"]),
              .fun = function(id) print(translate.submitterId(id)))
        error.fun("Some competitive objectIds are not in post-competitive results\n")
    }
    
    flag <- (res$comparator == TRUE) | (res$objectId %in% res.comp$objectId)
    if(!all(flag)) {
        error.fun("Some post-competitive objectIds are not in competitive results\n")
    }

    cat("Merging with res.comp\n")
    res <- safe.merge(res, unique(res.comp[, c("objectId", "submitterId")]), all.x=TRUE)

    ## Assign the team name

    cat("Defining team name tbl\n")
    team.name.tbl <- unique(res[, c("objectId", "subchallenge", "submitterId", "comparator")])
    flag <- !(team.name.tbl$comparator == TRUE)
    
    team.name.tbl$method.name <- NA
    team.name.tbl[flag, "method.name"] <-
        unlist(lapply(as.character(team.name.tbl[flag, "submitterId"]),
                      function(str) translate.submitterId(str)))
    
    team.name.tbl <-
        assign.baseline.names(team.name.tbl, from.col = "objectId", to.col = "method.name")
    
    
    team.name.tbl <- simplify.submitter.names(team.name.tbl, col = "method.name")

    ## Assign the round
    team.name.tbl <-
        assign.submission.rounds(team.name.tbl, object.id.col = "objectId",
                                 context.cols = c("subchallenge", "submitterId"),
                                 method.name.col = "method.name")

    cat("Merging team name tbl\n")
    print(colnames(team.name.tbl))
    res <- merge(res, team.name.tbl, all.x = TRUE)
    cat("Done merging team name tbl\n")
    
    res
}

get.comparators <- function() {
    c("cibersort", "mcp", "quantiseq", "xcell", "epic", "timer", "cibersortx")
}

get.comparators.cap <- function() {
    c("CIBERSORT", "MCP-counter", "quanTIseq", "xCell", "EPIC", "TIMER", "CIBERSORTx")
}

get.baseline.translation <- function() {
    baseline.method.trans <-
        list("baseline_method1" = "CIBERSORT",
             "baseline_method2" = "MCP-counter",
             "baseline_method3" = "quanTIseq",
             "baseline_method4" = "xCell",
             "baseline_method5" = "EPIC",
             "baseline_method6" = "TIMER",
             "baseline_method7" = "CIBERSORTx",
             "cibersortx" = "CIBERSORTx",
             "CIBERSORTx" = "CIBERSORTx",             
             "cibersort" = "CIBERSORT",
             "mcp" = "MCP-counter",
             "quantiseq" = "quanTIseq",
             "xcell" = "xCell",
             "epic" = "EPIC",
             "timer" = "TIMER")
    baseline.method.trans
}

assign.baseline.names <- function(df, from.col = "repo_name", to.col = "repo_name") {
    baseline.method.trans <- get.baseline.translation()
    for(nm in names(baseline.method.trans)) {
        flag <- (grepl(df[, from.col], pattern = nm))
        df[flag, to.col] <- baseline.method.trans[[nm]]
    }
    df
}

define.baseline.method.flag <- function(res, method.name.col) {
    baseline.method.trans <- get.baseline.translation()    
    baseline.methods <- unname(unlist(baseline.method.trans))
    baseline.method.flag <- grepl(res[, method.name.col], pattern="baseline")
    for(baseline.method in baseline.methods) {
        baseline.method.flag <- baseline.method.flag |
            grepl(res[, method.name.col], pattern=baseline.method, ignore.case=TRUE)
    }
    baseline.method.flag
}

## Assign rounds to assignments based on objectId < object Id Y then X was
## submitted before Y

## context.cols: the combination of columns that select all versions/submissions of a method
## i.e., for the raw prediction results, context.cols = c(subchallenge.col, submitter.id.col)
## For the leaderboard results (already separated by leaderboard), it is
## context.cols = c(submitter.id.col)
## method.name.col: probably repo_name
assign.submission.rounds <- function(res, object.id.col, context.cols, method.name.col,
                                     assign.latest = TRUE) {
    baseline.method.flag <-
        define.baseline.method.flag(res, method.name.col)
    tmp <- unique(res[!baseline.method.flag, c(context.cols, object.id.col)])
    ret <-
        ddply(tmp,
              .variables = context.cols,
              .fun = function(df) {
                  if(length(unique(df[, object.id.col])) != nrow(df)) {
                      stop("Was only expecting unique object IDs\n")
                  }
                  o <- order(df[, object.id.col], decreasing = FALSE)
                  df <- df[o, ]
                  ret.df <- data.frame(df, submission = 1:nrow(df))
                  if(assign.latest) {
                      ret.df <- rbind(ret.df, ret.df[nrow(ret.df),,drop=F])
                      ret.df[nrow(ret.df), "submission"] <- "latest"
                  }
                  ret.df
              })
    ret <- merge(res[!baseline.method.flag, ], ret, all.x = TRUE)

    ## Handle baselines, where all submissions of the same method
    ## should have the same method name
    tmp <- unique(res[baseline.method.flag, c(context.cols, object.id.col, method.name.col)])
    baseline.ret <-
        ddply(tmp,
              .variables = c(context.cols, method.name.col),
              .fun = function(df) {
                  if(length(unique(df[, object.id.col])) != nrow(df)) {
                      stop("Was only expecting unique object IDs\n")
                  }
                  o <- order(df[, object.id.col], decreasing = FALSE)
                  df <- df[o, ]
                  ret.df <- data.frame(df, submission = 1:nrow(df))
                  if(assign.latest) {
                      ret.df <- rbind(ret.df, ret.df[nrow(ret.df),,drop=F])
                      ret.df[nrow(ret.df), "submission"] <- "latest"
                  }
                  ret.df
              })
    baseline.ret <- merge(res[baseline.method.flag, ], baseline.ret, all.x = TRUE)    

    rbind(ret, baseline.ret)
}

simplify.submitter.names <- function(df, col = "submitter") {
    df[, col] <- as.character(df[, col])
    trans <-
        list("Northwestern Polytechnical University" = "NPU",
             "TJU and the renegade mouse" = "TJU")
    for(nm in names(trans)) {
        flag <- grepl(df[, col], pattern = nm)
        df[flag, col] <- trans[[nm]]
    }
    df
    
}

translate.submitterId <- function(submitterId) {
    tryCatch(synGetTeam(submitterId)$name,
             error = function(e) synGetUserProfile(submitterId)$userName)
}

ggplot_smooth_scatter <- function(data, mapping, pwr = 0.25, n = 200){
      p <- ggplot(data = data, mapping = mapping) + 
        stat_density2d(aes(fill=..density..^pwr), geom="tile", contour = FALSE, n = n) +
        scale_fill_continuous("", low = "white", high = "dodgerblue4")
##        scale_fill_gradientn(colours=rainbow(100))
      p
}

## Calculate the mean / variance loess fit for the genes (rows) in mat
## Expression values should be in log space
calc.mean.variance <- function(mat) {
  means <- unlist(apply(mat, 1, mean))
  sds <- unlist(apply(mat, 1, sd))
  x <- means
  y <- sqrt(sds)
  loess.fit <- loess(y~x,span=0.3,degree=1,family="symmetric",iterations=4,surface="direct")
  se <- FALSE
  loess.pred <- predict(loess.fit, data.frame(x=as.numeric(x)), se=se)

  loess.df <- data.frame(mean = x, sqrt.std = y, 
                         residual = as.numeric(residuals(loess.fit)))
  if(se) {
    loess.df$sqrt.std.pred = as.numeric(loess.pred$fit)
  } else {
    loess.df$sqrt.std.pred = as.numeric(loess.pred)
  }

  ## See from stats:::predLoess that se.fit is s(x) = s * \sqrt(\sum_{i=1}^n l_i^2(x)) on page 4 of Cleveland and Grosse (1991) 
  ## Statistics and Computing.  There the authors state that [ g^(x) - g(x) ] / s(x) has a t-distribution.
  ## g(x) is the original function and g^(x) is the loess prediction, so that g^(x) - g(x) is the residual.
  ## We considered selecting highly-variable genes based on the standardized residuals, but this penalizes areas of the
  ## curve in which the prediction is poor because of few points.  So, just defined highly variable points based
  ## on residual.
  ##
  ## Define standardized residuals
  if(se) { 
    loess.df$se.fit <- as.numeric(less.pred$se.fit)
    loess.df$residual.std <- loess.df$residual / loess.df$se.fit
  }
  loess.df <- loess.df[order(loess.df$mean), ]

  loess.df
}


## Assumes that expressed genes are represented by the maximal-density
## peak to the right of zero.
z_fpkm<-function(log.expr){
  if(all(log.expr>0)) stop('Input not log2 transformed.')
  pos.log.expr <- log.expr[log.expr > 0]
  density.peak <- max(density(pos.log.expr,na.rm=T)$y)
  my<-density(pos.log.expr,na.rm=T)$x[which.max(density(pos.log.expr,na.rm=T)$y)]
  U<-mean(log.expr[log.expr>my],na.rm=T)
  sigma<-(U-my)*(.5*pi)^.5
  z<-(log.expr-my)/sigma
  z[z< -3]<-NA
  return(list("mean" = my, "sd" = sigma, "z" = z, "density.peak" = density.peak))
}

z_fpkm_rightmost_peak <-function(log.expr, min.boundary = 0, max.boundary = 6){
  if(all(log.expr>0)) stop('Input not log2 transformed.')
  use.derivative.method <- TRUE
  ## I use this first derivative approach, because otherwise I need to closely
  ## bracket the peak (between min/max boundaries).  e.g., the max density
  ## in TCGA sample TCGA−AB−2929−03A−01T−0735−13 occurs at 0, whereas the
  ## desired peak occurs at 5. The latter is more easily found using the
  ## derivative approach. The max density approach would require setting the
  ## min boundary in the trough between peaks or up the left hand side of the
  ## max density peak.
  my <- NA
  density.peak <- NA
  den <- density(log.expr, na.rm=TRUE)
  df <- data.frame(x = den$x, y = den$y)
  if(use.derivative.method) {
    sp <- splinefun(den$x, den$y)
    df$deriv = sp(den$x, deriv=1)
    df <- subset(df, (x > min.boundary) & (x < max.boundary))
    for(i in (nrow(df)-1):1) {
      if(sign(df$deriv[i+1]) != sign(df$deriv[i])) {
        my <- df$x[i]
        density.peak <- df$y[i]
        break
      }
    }
  } else {
    df <- subset(df, (x > min.boundary) & (x < max.boundary))
    flag <- which.max(df$y)
    density.peak <- df$y[flag]
    my <- df$x[flag]
  }
  U<-mean(log.expr[log.expr>my],na.rm=T)
  sigma<-(U-my)*(.5*pi)^.5
  z<-(log.expr-my)/sigma
  z[z< -3]<-NA
  return(list("mean" = my, "sd" = sigma, "z" = z, "density.peak" = density.peak))
}

calc.density.of.expressed.genes <- function(log.expr, num.sds = NULL, z = NULL) {
  if(is.null(z)) {
    z <- z_fpkm(log.expr)
  }
  den <- density(log.expr)
  data.density <- data.frame(x = den$x, y = den$y)
  dn <- dnorm(den$x, mean = z$mean, sd = z$sd)
  density.peak <- den$y[which.min((den$x - z$mean)^2)]
  fit <- data.frame(x = den$x, y = density.peak * dn / max(dn)) 
  lst <- list("density" = data.density, "fit" = fit, "z" = z)
  lst
}

plot.density.of.expressed.genes <- function(log.expr, num.sds = NULL, z = NULL) {
  suppressPackageStartupMessages(p_load(latex2exp))
  lst <- calc.density.of.expressed.genes(log.expr, num.sds, z)
  data.density <- lst[["density"]]
  fit <- lst[["fit"]]
  z <- lst[["z"]]
  g <- ggplot()
  g <- g + geom_line(data = data.density, aes(x = x, y = y))
  g <- g + geom_line(data = fit, aes(x = x, y = y), linetype = "dashed", col = "blue")
  if(!is.null(num.sds)) {
    g <- g + geom_vline(xintercept = z$mean - num.sds * z$sd, linetype = "dashed")
  }
##  g <- g + xlab("Log CPM")
  g <- g + xlab(TeX("$Log_2(CPM)$"))
  g <- g + ylab("Density")
  g
}


get.deconvolution.genes <- function() {
    suppressPackageStartupMessages(p_load("immunedeconv"))
    suppressPackageStartupMessages(p_load("MCPcounter"))

    xcell.deconv.genes <-
    sort(unique(unlist(lapply(xCell.data$signatures, function(x) x@geneIds))))

    lm22 <- NULL
    files <- c("/home/bwhite/LM22.txt", "/Users/Brian/Downloads/new-cibersort-code/LM22.txt")
    for(file in files) {
        if(file.exists(file)) {
            lm22 <- read.table(file, sep="\t", header=TRUE)
        }
    }
    cibersort.deconv.genes <- as.character(lm22$Gene.symbol)
    mcp.deconv.genes <-
        sort(unique(
            read.table(curl:::curl("http://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/genes.txt"), 
                       sep = "\t", stringsAsFactors = FALSE, header = TRUE, 
                       colClasses = "character", check.names = FALSE)$`HUGO symbols`))
    
    epic.deconv.genes <- sort(unique(c(TRef$sigGenes, BRef$sigGenes)))
    
    signame <- "TIL10"
    sig.mat.file <- system.file("extdata", "quantiseq",
                                paste0(signame, 
                                       "_signature.txt"), package = "immunedeconv", mustWork = TRUE)
    df <- read.table(sig.mat.file, sep="\t", header=TRUE)
    quantiseq.deconv.genes <- as.character(unique(df$ID))

    deconv.genes <-
        sort(unique(c(cibersort.deconv.genes, mcp.deconv.genes, epic.deconv.genes, quantiseq.deconv.genes)))

    deconv.genes
}

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

fix.string <- function(str) {
    str <- gsub(str, pattern="\\.", replacement=" ")
    str <- unlist(strsplit(str, split=" "))
    str <- lapply(str, function(x) firstup(x))
    str <- paste0(str, collapse=" ")
    str
}

get.in.silico.metadata <- function() {

    synId <- "syn21773017"
    metadata.file <- synGet(synId, downloadFile = TRUE)$path
    df <- read.table(metadata.file, sep = ",", header = TRUE, stringsAsFactors = FALSE)
    return(df)
    
    l <-
        list(
            "A" = list("mixture.type" = "Random", "tumor.type" = NA, "subchallenge" = "fine"),
            "B" = list("mixture.type" = "Biological", "tumor.type" = NA, "subchallenge" = "fine"),            
            "G" = list("mixture.type" = "Random", "tumor.type" = NA, "subchallenge" = "coarse"),
            "H" = list("mixture.type" = "Biological", "tumor.type" = NA, "subchallenge" = "coarse"),
            "C" = list("mixture.type" = "Random", "tumor.type" = "CRC", "subchallenge" = "fine"),
            "E" = list("mixture.type" = "Biological", "tumor.type" = "CRC", "subchallenge" = "fine"),            
            "D" = list("mixture.type" = "Random", "tumor.type" = "BRCA", "subchallenge" = "fine"),
            "F" = list("mixture.type" = "Biological", "tumor.type" = "BRCA", "subchallenge" = "fine"),            
            "I" = list("mixture.type" = "Random", "tumor.type" = "CRC", "subchallenge" = "coarse"),
            "K" = list("mixture.type" = "Biological", "tumor.type" = "CRC", "subchallenge" = "coarse"),            
            "J" = list("mixture.type" = "Random", "tumor.type" = "BRCA", "subchallenge" = "coarse"),
            "L" = list("mixture.type" = "Biological", "tumor.type" = "BRCA", "subchallenge" = "coarse")
        )
    df <- ldply(l, .fun = function(entry) as.data.frame(entry, stringsAsFactors = FALSE))
    colnames(df)[1] <- "dataset.name"
    df
}

get.validation.metadata <- function() {
    ## Our original admixture specification includes the vendor for each sample
    synId <- "syn21577258"
    obj <- synGet(synId, downloadFile = TRUE)
    bm1 <- read.xlsx(obj$path, sheet = 1)
    bm1 <- bm1[-c(1:3),]
    bm2 <- read.xlsx(obj$path, sheet = 2)
    bm2 <- bm2[-c(1:3),]
    rm1 <- read.xlsx(obj$path, sheet = 3)
    rm1 <- rm1[-c(1:3),]
    rm2 <- read.xlsx(obj$path, sheet = 4)
    rm2 <- rm2[-c(1:3),]
    
    lst <- list(bm1, bm2, rm1, rm2)

    df <-
        do.call(rbind, lapply(lst,
                              function(entry) {
                                  label = ifelse(entry$CRC > 0, "CRC",
                                          ifelse(entry$breast > 0, "BRCA", ""))                              
                                  data.frame(id = entry[,1], tumor.type = label)
                              }))
    df$batch <- "empty"
    flag <- df$id %in% bm1[,1]
    df$batch[flag] <- "BM1"
    flag <- df$id %in% bm2[,1]
    df$batch[flag] <- "BM2"
    flag <- df$id %in% rm1[,1]
    df$batch[flag] <- "RM1"
    flag <- df$id %in% rm2[,1]
    df$batch[flag] <- "RM2"
    
    df$mixture.type <- "empty"
    flag <- ( df$id %in% bm1[,1] ) | ( df$id %in% bm2[,1] )
    ##    df$mixture.type[flag] <- "BM"
    df$mixture.type[flag] <- "Biological"
    flag <- ( df$id %in% rm1[,1] ) | ( df$id %in% rm2[,1] )
    ## df$mixture.type[flag] <- "RM"
    df$mixture.type[flag] <- "Random"
    
    df$dataset <- "empty"
    flag <- ( df$tumor.type == "BRCA" ) & ( df$mixture.type == "Biological" )
    df$dataset[flag] <- "DS1"
    flag <- ( df$tumor.type == "CRC" ) & ( df$mixture.type == "Biological" )
    df$dataset[flag] <- "DS2"
    flag <- ( df$tumor.type == "BRCA" ) & ( df$mixture.type == "Random" )
    df$dataset[flag] <- "DS3"
    flag <- ( df$tumor.type == "CRC" ) & ( df$mixture.type == "Random" )
    df$dataset[flag] <- "DS4"

    df
}

## This is the ground truth we _specified_ to the Stanford core, not that we actually had created
get.specified.validation.ground.truth <- function() {
    ## Our original admixture specification includes the vendor for each sample
    synId <- "syn21577258"
    obj <- synGet(synId, downloadFile = TRUE)
    bm1 <- read.xlsx(obj$path, sheet = 1)
    bm1 <- bm1[-c(1:3),]
    bm2 <- read.xlsx(obj$path, sheet = 2)
    bm2 <- bm2[-c(1:3),]
    rm1 <- read.xlsx(obj$path, sheet = 3)
    rm1 <- rm1[-c(1:3),]
    rm2 <- read.xlsx(obj$path, sheet = 4)
    rm2 <- rm2[-c(1:3),]
    
    lst <- list(bm1, bm2, rm1, rm2)

    df <-
        do.call(rbind, lapply(lst,
                              function(entry) {
                                  rownames(entry) <- entry[,1]
                                  entry <- entry[,-1]
                                  ret <- melt(as.matrix(entry))
                                  colnames(ret) <- c("sample.id", "cell.type", "measured")
                                  ret$sample.id <- as.character(ret$sample.id)
                                  ret$cell.type <- as.character(ret$cell.type)
                                  ret$measured <- as.numeric(as.character(ret$measured))
                                  ret
                              }))

    md <- get.validation.metadata()
    colnames(md) <- c("sample.id", "tumor.type", "batch", "mixture.type", "dataset.name")
    df <- merge(df, md)
    
    df
}

get.purified.sample.translation.table <- function() {
    fine.grained.translation <- list(
        "Naive_B_cells_1" = "naive.B.cells",
        "Macrophages_2" = "macrophages",
        "Dendritic_cells_1" = "myeloid.dendritic.cells",
        "Macrophages_1" = "macrophages",
        "Fibroblasts" = "fibroblasts",
        "Endothelial_cells" = "endothelial.cells",
        "Memory_CD8_T_cells_1" = "memory.CD8.T.cells",
        "Monocytes_1" = "monocytes",
        "Neutrophils_2" = "neutrophils",
        "Dendritic_cells_2" = "myeloid.dendritic.cells",
        "NK_cells_2" = "NK.cells",
        "Monocytes_2" = "monocytes",
        "NK_cells_1" = "NK.cells",
        "Memory_CD8_T_cells_2" = "memory.CD8.T.cells",
        "Memory_CD4_T_cells_2" = "memory.CD4.T.cells",
        "Naive_CD4_T_cells_1" = "naive.CD4.T.cells",
        "Tregs" = "regulatory.T.cells",
        "Memory_CD4_T_cells_1" = "memory.CD4.T.cells",
        "Naive_CD8_T_cells_2" = "naive.CD8.T.cells",
        "Naive_CD4_T_cells_2" = "naive.CD4.T.cells")

    coarse.grained.translation <- list(
        "Naive_B_cells_1" = "B.cells",
        "Macrophages_2" = "monocytic.lineage",
        "Dendritic_cells_1" = "monocytic.lineage",
        "Macrophages_1" = "monocytic.lineage",
        "Fibroblasts" = "fibroblasts",
        "Endothelial_cells" = "endothelial.cells",
        "Memory_CD8_T_cells_1" = "CD8.T.cells",
        "Monocytes_1" = "monocytic.lineage",
        "Neutrophils_2" = "neutrophils",
        "Dendritic_cells_2" = "monocytic.lineage",
        "NK_cells_2" = "NK.cells",
        "Monocytes_2" = "monocytic.lineage",
        "NK_cells_1" = "NK.cells",
        "Memory_CD8_T_cells_2" = "CD8.T.cells",
        "Memory_CD4_T_cells_2" = "CD4.T.cells",
        "Naive_CD4_T_cells_1" = "CD4.T.cells",
        "Tregs" = "CD4.T.cells",
        "Memory_CD4_T_cells_1" = "CD4.T.cells",
        "Naive_CD8_T_cells_2" = "CD8.T.cells",
        "Naive_CD4_T_cells_2" = "CD4.T.cells")

    fine.grained.tbl <- data.frame(sample = names(fine.grained.translation),
                                   fine.grained.cell.type = unname(unlist(fine.grained.translation)))

    coarse.grained.tbl <- data.frame(sample = names(coarse.grained.translation),
                                   coarse.grained.cell.type = unname(unlist(coarse.grained.translation)))
    tbl <- merge(fine.grained.tbl, coarse.grained.tbl, by = "sample")
    tbl
}




    
## This is the ground truth that the Stanford core actually created
## NB: this function returns ground truth specified at the level of cell types
## (e.g., naive.B.cells) not sample names (e.g., Naive_B_cells_2)
get.actual.validation.ground.truth <- function() {

    ## Get the final admixture ratios
    ## Note these rows to not sum to 1, so need to renormalize
    ## This is the old/first "ratio" file that John Coller sent us "Final_Mixture_Ratios.xlsx"
    ## Instead use the subsequent file he sent "Final_Mixture_Fractions.xlsx" in which
    ## he appears to have just decided the rows by the row sum (which is what I do for the ratio file below).
    if(FALSE) {
        synId <- "syn21577248"
        obj <- synGet(synId, downloadFile = TRUE)
        stop("Need to read _both_ sheets -- this code is not doing that. See below for fractions file")
        ratios <- read.xlsx(obj$path, sheet = 1, startRow = 2)
        ratios <- ratios[-1,]
        nms <- ratios$X1
        ratios <- ratios[,-1]
        ratios <-
            ldply(1:nrow(ratios), .fun = function(i) ratios[i,] / sum(ratios[i,]))
        rownames(ratios) <- nms
    }
    
    synId <- "syn21598638"
    obj <- synGet(synId, downloadFile = TRUE)

    eps <- 10^-4
    
    old.col.names <- c("breast", "CRC", "Fibroblasts", "Endothelial_cells", "Dendritic_cells",
                       "Monocytes", "Macrophages", "NK_cells", "Tregs", "Naive_CD4_T_cells",
                       "Memory_CD4_T_cells", "Memory_CD8_T_cells", "Naive_B_cells")
    new.col.names <- c("BRCA", "CRC", "fibroblasts", "endothelial.cells", "DC", "monocytes",
                       "macrophages", "NK.cells", "regulatory.T.cells", "naive.CD4.T.cells",
                       "memory.CD4.T.cells", "memory.CD8.T.cells", "naive.B.cells")


    
    ratios <-
        ldply(1:2,
              .fun = function(sheetIndex) {
                  df <- read.xlsx(obj$path, sheet = sheetIndex, startRow = 2)
                  df <- df[-1,]
                  nms <- df$X1
                  df <- df[,-1]
                  df <- df[, old.col.names]
                  colnames(df) <- new.col.names
                  df$sample <- nms
                  df
              })
    
    rownames(ratios) <- ratios$sample
    ratios <- ratios[, !(colnames(ratios) == "sample")]
    for(col in 1:ncol(ratios)) { ratios[,col] <- as.numeric(ratios[,col]) }
    
    l_ply(1:nrow(ratios),
          .fun = function(i) if(abs(sum(as.numeric(ratios[i,])) - 1) > eps) { stop("Rows don't sum to one") })
    
    melted.ratios <- melt(as.matrix(ratios))
    colnames(melted.ratios) <- c("sample.id", "cell.type", "measured")
    melted.ratios$cell.type <- as.character(melted.ratios$cell.type)

    md <- get.validation.metadata()
    colnames(md) <- c("sample.id", "tumor.type", "batch", "mixture.type", "dataset.name")
    melted.ratios <- merge(melted.ratios, md)
    
    return(melted.ratios)
}

## This is the ground truth that the Stanford core actually created
## NB: this function returns ground truth specified at the level of
## sample names (e.g., Naive_B_cells_2) not cell types (e.g., naive.B.cells) 
get.actual.validation.ground.truth.sample.names <- function() {
    ## Get the final admixture ratios
    ## Note these rows to not sum to 1, so need to renormalize
    ## This is the old/first "ratio" file that John Coller sent us "Final_Mixture_Ratios.xlsx"
    ## Instead use the subsequent file he sent "Final_Mixture_Fractions.xlsx" in which
    ## he appears to have just decided the rows by the row sum (which is what I do for the ratio file below).
    if(FALSE) {
        synId <- "syn21577248"
        obj <- synGet(synId, downloadFile = TRUE)
        stop("Need to read _both_ sheets -- this code is not doing that. See below for fractions file")
        ratios <- read.xlsx(obj$path, sheet = 1, startRow = 2)
        ratios <- ratios[-1,]
        nms <- ratios$X1
        ratios <- ratios[,-1]
        ratios <-
            ldply(1:nrow(ratios), .fun = function(i) ratios[i,] / sum(ratios[i,]))
        rownames(ratios) <- nms
    }

    synId <- "syn21598638"
    obj <- synGet(synId, downloadFile = TRUE)
    old.col.names <- c("breast", "CRC", "Fibroblasts", "Endothelial_cells", "Dendritic_cells",
                       "Monocytes", "Macrophages", "NK_cells", "Tregs", "Naive_CD4_T_cells",
                       "Memory_CD4_T_cells", "Memory_CD8_T_cells", "Naive_B_cells")
    
    exclude.cols <- c("Tregs", "Endothelial_cells", "breast", "CRC", "Fibroblasts")
    
    ## Leave off the second sheet (note 1:1 below) since Naive_B_cells_2
    ## appear to have been used in admixtures, but not sequenced.
    sample.ratios <-
        ldply(1:1,
              .fun = function(sheetIndex) {
                  df <- read.xlsx(obj$path, sheet = sheetIndex, startRow = 2)
                  df <- df[-1,]
                  nms <- df$X1
                  df <- df[,-1]
                  df <- df[, old.col.names]
                  colnames(df) <-
                      unlist(lapply(colnames(df),
                                    function(str) ifelse(str %in% exclude.cols,
                                                         str,
                                                         paste0(str, "_", sheetIndex))))
                  rownames(df) <- nms
                  for(col in 1:ncol(df)) { df[,col] <- as.numeric(df[,col]) }		 
                  m <- melt(as.matrix(df))
                  colnames(m) <- c("sample", "cell.type", "actual")
                  m
                  
              })
    sample.ratios$cell.type <- as.character(sample.ratios$cell.type)
    sample.ratios$sample <- as.character(sample.ratios$sample)
    
    sample.ratios[sample.ratios$cell.type == "breast", "cell.type"] <- "Breast"
    eps <- 10^-4
    chk <- ddply(sample.ratios, .variables = c("sample"), .fun = function(df) sum(df$actual))
    if(any(abs(chk$V1 - 1) > eps)) { stop("Something doesn't sum to one") }
    sample.ratios
}

## Our original admixture specification includes the vendor for each sample
## Return a data frame with columns sample and vendor, given the vendor associated with each sample
get.vendor.sample.assignments <- function() {
    synId <- "syn21577258"
    obj <- synGet(synId, downloadFile = TRUE)
    vendors1 <- read.xlsx(obj$path, sheet = 1)[3,,drop=F]
    vendors2 <- read.xlsx(obj$path, sheet = 2)[3,,drop=F]
    vendors1 <- vendors1[,-1,drop=F]
    vendors2 <- vendors2[,-1,drop=F]
    exclude.cols <- c("Tregs", "Endothelial_cells", "breast", "CRC", "Fibroblasts")
    colnames(vendors1) <-
        unlist(lapply(colnames(vendors1),
                      function(str) ifelse(str %in% exclude.cols,
                                           str,
                                           paste0(str, "_1"))))
    colnames(vendors2) <-
        unlist(lapply(colnames(vendors2),
                      function(str) ifelse(str %in% exclude.cols,
                                           str,
                                           paste0(str, "_2"))))
    vendors <- cbind(vendors1, vendors2[, !(colnames(vendors2) %in% exclude.cols), drop=F])
    vendors <- t(rename.samples(vendors))
    rownames(vendors)[grepl(rownames(vendors), pattern="breast")] <- "Breast"
    colnames(vendors) <- c("vendor")
    vendors <- data.frame(sample = rownames(vendors), vendors)
    vendors
}

plot.admixtures <- function(mat) {

    hc <- hclust(dist(mat))
    labels <- hc$labels[hc$order]
    if("CRC" %in% labels) { labels <- c("CRC", labels[labels != "CRC"]) }
    if("breast" %in% labels) { labels <- c("breast", labels[labels != "breast"]) }
    if("Breast" %in% labels) { labels <- c("Breast", labels[labels != "Breast"]) }
    if("BRCA" %in% labels) { labels <- c("BRCA", labels[labels != "BRCA"]) }        
        
    mat <- mat[labels,]
    cols <- c("CRC", "breast", "Breast", "BRCA")
    cols <- cols[cols %in% labels]
    if(length(cols) == 2) {
        o <- order(mat[cols[1],,drop=T], mat[cols[2],,drop=T])
        mat <- mat[, o]
    } else {
        hc <- hclust(dist(t(mat)))
        mat <- mat[,hc$labels[hc$order]]
    }
    
    row.levels <- rownames(mat)
    col.levels <- colnames(mat)

    df <- melt(mat)
    df$Var1 <- factor(df$Var1, levels = row.levels)
    df$Var2 <- factor(df$Var2, levels = col.levels)
    
    g <- ggplot(data = df, aes(x = Var2, y = Var1, fill = value))
    g <- g + geom_tile()
    g <- g + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                   text = element_text(size=15))
    g <- g + scale_fill_gradient2("Proportion", limits = c(0,1), low = "red", high = "blue", mid = "white", na.value = "black")
    g <- g + xlab("Admixture") + ylab("Cell Type")
    g <- g + theme(axis.text.x = element_blank())
    g
}

# whether or not to order methods in the heatmap using both a summary statistic (e.g., mean)
# that excludes nulls and one that does not or only to use that excludes or includes ...
use.include.and.exclude.null.in.summary <- FALSE
# ... and if only one, which to use
na.rm.summary <- TRUE


calculate.method.levels <- function(df, id.var = "modelId", cell.type.var = "cell.type", cor.var = "cor.p",
                                    row.summary.fun = "mean", col.summary.fun = "max",
                                    order.decreasing = FALSE) {

    orig.df <- df
    df <- df[, c(id.var, cell.type.var, cor.var)]
    df[, id.var] <- as.character(df[, id.var])
    df[, cell.type.var] <- as.character(df[, cell.type.var])

    # Calculate the method summary over cell types (e.g., mean) with and without excluding NAs (if use.include.and.exclude.null.in.summary == TRUE)
    # e.g., "mean" = mean without excluding NAs
    #       "nnmean" = non-null mean (mean excluding NAs)
    # i.e., always prepend "nn" to summary
    # Order NAs from highest to lowest based on nnmean and then non-NAs from highest to lowest based on mean.
    # Admittedly, this is a strange rule, but I like the end result in the heatmaps -- it makes a block of
    # the comparators at the top that have NAs, with the next score the top performing method (comparator or
    # or otherwise). This makes that top performer stand out more than it would if the blanks/NAs were 
    # scattered throughout.
    nnrow.summary.fun <- paste0("nn", row.summary.fun)
    method.summaries <- ddply(df, .variables = c(id.var),
                          .fun = function(tmp) {
                              if(use.include.and.exclude.null.in.summary) {
                                nret <- data.frame(cell.type = row.summary.fun, cor = do.call(row.summary.fun, list(tmp[, cor.var], na.rm=FALSE)))
                                nnret <- data.frame(cell.type = nnrow.summary.fun, cor = do.call(row.summary.fun, list(tmp[, cor.var], na.rm=TRUE)))
			        ret <- rbind(nret, nnret)
                              } else {
                                ret <- data.frame(cell.type = row.summary.fun, cor = do.call(row.summary.fun, list(tmp[, cor.var], na.rm=na.rm.summary)))
                              }
                              colnames(ret) <- c(cell.type.var, cor.var)
                              ret
                          })

    if(use.include.and.exclude.null.in.summary) {
      na.methods <- method.summaries[is.na(method.summaries[, cor.var]),id.var]
      na.summaries <- method.summaries[(method.summaries[, cell.type.var] == nnrow.summary.fun) & (method.summaries[, id.var] %in% na.methods),]
      na.summaries <- na.summaries[order(na.summaries[, cor.var], decreasing = order.decreasing),]
      flag1 <- !is.na(method.summaries[, cor.var])
      flag2 <- !(method.summaries[, id.var] %in% na.summaries[, id.var])
      flag3 <- (method.summaries[, cell.type.var] == row.summary.fun)
      non.na.summaries <- method.summaries[flag1 & flag2 & flag3, ]
      non.na.summaries <- non.na.summaries[order(non.na.summaries[, cor.var], decreasing = order.decreasing),]
      method.summaries <- rbind(non.na.summaries, na.summaries)
    } else {
      method.summaries <- method.summaries[order(method.summaries[, cor.var], decreasing = order.decreasing),]
    }
    method.levels <- method.summaries[method.summaries[, id.var] != col.summary.fun, id.var]

    method.levels
}

calculate.cell.type.levels <- function(df, id.var = "modelId", cell.type.var = "cell.type", cor.var = "cor.p",
                                       row.summary.fun = "mean", col.summary.fun = "max",
                                       order.decreasing = FALSE) {

    orig.df <- df
    df <- df[, c(id.var, cell.type.var, cor.var)]
    df[, id.var] <- as.character(df[, id.var])
    df[, cell.type.var] <- as.character(df[, cell.type.var])

    na.rm <- FALSE
    
    cell.type.summaries <- ddply(df, .variables = cell.type.var,
                             .fun = function(tmp) {
                                 ret <- data.frame(id = col.summary.fun, cor = do.call(col.summary.fun, list(tmp[, cor.var], na.rm=TRUE)))
                                 colnames(ret)[1] <- id.var
                                 colnames(ret)[2] <- cor.var                                 
                                 ret
                             })
    cell.type.summaries <- cell.type.summaries[order(cell.type.summaries[, cor.var], decreasing = order.decreasing),]
    
    cell.type.levels <- cell.type.summaries[cell.type.summaries[, cell.type.var] != row.summary.fun, cell.type.var]

    cell.type.levels
}


plot.cell.type.correlation.heatmap <- function(df, show.corr.text = FALSE, id.var = "modelId", cell.type.var = "cell.type", cor.var = "cor.p",
                                               cor.type.label = "Pearson\nCorrelation", limits = c(-1, 1),
                                               pval.var = NULL, row.summary.fun = "mean", col.summary.fun = "max",
                                               order.decreasing = FALSE,
                                               method.levels = NULL,
                                               cell.type.levels = NULL,
                                               formatter = function(x) formatC(x, format="f", digits=2)) {
    orig.df <- df
    df <- df[, c(id.var, cell.type.var, cor.var)]
    df[, id.var] <- as.character(df[, id.var])
    df[, cell.type.var] <- as.character(df[, cell.type.var])

    ## Add NAs for any missing entries
    df <- acast(df, as.formula(paste0(id.var, " ~ ", cell.type.var)), value.var = cor.var, fill = NA)
    df <- reshape2::melt(as.matrix(df))
    colnames(df) <- c(id.var, cell.type.var, cor.var)
    df[, id.var] <- as.character(df[, id.var])
    df[, cell.type.var] <- as.character(df[, cell.type.var])
    
    # na.rm <- FALSE
    cell.type.summaries <- ddply(df, .variables = cell.type.var,
                             .fun = function(tmp) {
                                 ret <- data.frame(id = col.summary.fun, cor = do.call(col.summary.fun, list(tmp[, cor.var], na.rm=TRUE)))
                                 colnames(ret)[1] <- id.var
                                 colnames(ret)[2] <- cor.var                                 
                                 ret
                             })
    cell.type.summaries <- cell.type.summaries[order(cell.type.summaries[, cor.var], decreasing = order.decreasing),]
    
    # Calculate the method summary over cell types (e.g., mean) with and without excluding NAs (if use.include.and.exclude.null.in.summary == TRUE)
    # e.g., "mean" = mean without excluding NAs
    #       "nnmean" = non-null mean (mean excluding NAs)
    # i.e., always prepend "nn" to summary
    nnrow.summary.fun <- paste0("nn", row.summary.fun)
    method.summaries <- ddply(df, .variables = c(id.var),
                          .fun = function(tmp) {
                              if(use.include.and.exclude.null.in.summary) {
                                ret <- data.frame(cell.type = row.summary.fun, cor = do.call(row.summary.fun, list(tmp[, cor.var], na.rm=FALSE)))
                                nnret <- data.frame(cell.type = nnrow.summary.fun, cor = do.call(row.summary.fun, list(tmp[, cor.var], na.rm=TRUE)))
                                ret <- rbind(ret, nnret)
                              } else {
                                ret <- data.frame(cell.type = row.summary.fun, cor = do.call(row.summary.fun, list(tmp[, cor.var], na.rm=na.rm.summary)))
                              }
                              colnames(ret) <- c(cell.type.var, cor.var)
                              ret
                          })
    
    method.summaries <- method.summaries[order(method.summaries[, cor.var], decreasing = order.decreasing),]


    if(is.null(cell.type.levels)) {
        cell.type.levels <- c(grep(x = cell.type.summaries[, cell.type.var], pattern = row.summary.fun, values = TRUE), nnrow.summary.fun, row.summary.fun)
    } else {
        cell.type.levels <- c(cell.type.levels, nnrow.summary.fun, row.summary.fun)
    }
    
    ## method.levels <- c(col.summary.fun, method.summaries[method.summaries[, id.var] != col.summary.fun, id.var])
    if(is.null(method.levels)) {
        method.levels <- c(col.summary.fun, method.summaries[method.summaries[, id.var] != col.summary.fun, id.var])
    } else {
        method.levels <- c(col.summary.fun, method.levels)
    }
    
    
    df <- rbind(df, cell.type.summaries[, c(id.var, cell.type.var, cor.var)])
    df <- rbind(df, method.summaries[, c(id.var, cell.type.var, cor.var)])    

    ## df$cor.label <- formatC(df[, cor.var], format="f", digits=digits)
    ## formatter <- function(x) formatC(x, format="f", digits=2)
    df$cor.label <- formatter(df[, cor.var])    
    if(!is.null(pval.var)) {
        p_load(gtools)
        df <- merge(df, orig.df[, c(id.var, cell.type.var, pval.var)], all.x = TRUE)
        flag <- is.na(df[, pval.var])
        df[flag, pval.var] <- 1
        df$cor.label <- paste0(stars.pval(df[, pval.var]), "\n", df$cor.label)
    }
    df[, cell.type.var] <- factor(df[, cell.type.var], levels = cell.type.levels)
    df[, id.var] <- factor(df[, id.var], levels = method.levels)
    g <- ggplot(data = df, aes_string(y = id.var, x = cell.type.var, fill = cor.var))
    g <- g + geom_tile()
    if(show.corr.text) {
        g <- g + geom_text(aes(label = cor.label))
    }
    g <- g + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                   text = element_text(size=15))
    ## g <- g + ylab("Method") + xlab("")
    g <- g + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
    ## g <- g + scale_fill_continuous("Pearson\ncorrelation", limits = c(-1,1))
    ## g <- g + scale_fill_gradient2("Pearson\ncorrelation", limits = c(-1,1),
    ##                               low = "red", high = "blue", mid = "white", na.value = "black")
##    g <- g + scale_fill_gradient2(paste0(cor.type.label, "\nCorrelation"),
##                                  limits = c(-1,1), low = "red", high = "blue", mid = "white", na.value = "black")
    g <- g + scale_fill_gradient2(cor.type.label, limits = limits,
                                  low = "red", high = "blue", mid = "white", na.value = "black")
    ## g <- g + theme(text = element_text(size=20))
    g
}

proportion.labels <- function(x) ifelse(x == 0, "0", ifelse(x == 1, "1", x))

plot.strip.plots <- function(df, id.var = "modelId", cell.type.var = "cell.type", var = "cor.p",
                             label = "Pearson\nCorrelation", digits = 2, limits = c(-1, 1),
                             pval.var = NULL, col.summary.fun = "max", row.summary.fun = "mean", order.decreasing = FALSE,
                             method.levels = NULL,
                             cell.type.levels = NULL) {
    orig.df <- df
    df <- df[, c(id.var, cell.type.var, var)]
    df[, id.var] <- as.character(df[, id.var])
    df[, cell.type.var] <- as.character(df[, cell.type.var])

    na.rm <- TRUE
    if(is.null(cell.type.levels)) {
        cell.type.summaries <- ddply(df, .variables = cell.type.var,
                                     .fun = function(tmp) {
                                         ret <- data.frame(id = col.summary.fun, cor = do.call(col.summary.fun, list(tmp[, var], na.rm=TRUE)))
                                         colnames(ret)[1] <- id.var
                                         colnames(ret)[2] <- var                                 
                                         ret
                                     })
        cell.type.summaries <- cell.type.summaries[order(cell.type.summaries[, var], decreasing = order.decreasing),]
    
        ##    cell.type.levels <- c(cell.type.summaries[cell.type.summaries[, cell.type.var] != row.summary.fun, cell.type.var], row.summary.fun)
        cell.type.levels <- cell.type.summaries[cell.type.summaries[, cell.type.var] != row.summary.fun, cell.type.var]
    } 
        
    if(is.null(method.levels)) {
        method.summaries <- ddply(df, .variables = c(id.var),
                                  .fun = function(tmp) {
                                      ret <- data.frame(cell.type = row.summary.fun, cor = do.call(row.summary.fun, list(tmp[, var], na.rm=na.rm)))
                                      colnames(ret) <- c(cell.type.var, var)
                                      ret
                                  })
        
        method.summaries <- method.summaries[order(method.summaries[, var], decreasing = order.decreasing),]

        method.levels <- method.summaries[method.summaries[, id.var] != col.summary.fun, id.var]
    }


    df$cor.label <- formatC(df[, var], format="f", digits=digits)
    if(!is.null(pval.var)) {
        p_load(gtools)
        df <- merge(df, orig.df[, c(id.var, cell.type.var, pval.var)], all.x = TRUE)
        flag <- is.na(df[, pval.var])
        df[flag, pval.var] <- 1
        df$cor.label <- paste0(stars.pval(df[, pval.var]), "\n", df$cor.label)
    }
    cat("ggplot'ing\n")
    ## This re-ordering o works with as.table = FALSE to respect the
    ## the ordering we want from cell.type.levels
    ntot <- length(cell.type.levels)
    nrow <- 3
    ncol <- ceiling(ntot / nrow)
    o <- unlist(llply(1:(nrow-1), .fun = function(i) (i*ncol):(1+((i-1)*ncol))))
    o <- c(o, ntot:(max(o)+1))
    df[, cell.type.var] <- factor(df[, cell.type.var], levels = cell.type.levels[o])
    df[, id.var] <- factor(df[, id.var], levels = method.levels)
    g <- ggplot(data = df, aes_string(y = id.var, x = var))
    g <- g + geom_boxplot(outlier.shape = NA)
    g <- g + facet_wrap(cell.type.var, as.table = FALSE, nrow = nrow)
    ##    g <- g + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6), text = element_text(size=15))
    strip.text.sz <- 15
    if(ntot > 10) { strip.text.sz <- 8 }
    g <- g + theme(axis.text.y = element_text(size = 6), text = element_text(size=15), strip.text = element_text(size = strip.text.sz))
    
    g <- g + scale_x_continuous(labels = proportion.labels)
    g <- g + xlab(label) + ylab("")
    g
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
  sz <- 25
  g <- g + geom_text(x = xlimits[1] + xoffset * (xlimits[2] - xlimits[1]), y = ylimits[1] + 0.8 * (ylimits[2] - ylimits[1]), label = lm_corr_eqn(df, method = method, display.r2 = display.r2, display.pval = display.pval), parse=TRUE, ...)
  g <- g +theme(text = element_text(size = sz),
             axis.text.x = element_text(size=sz),
             axis.text.y = element_text(size=sz),
             axis.title.x = element_text(size=sz),
             axis.title.y = element_text(size=sz),
             title = element_text(size=sz),
             plot.title = element_text(hjust = 0.5, size=sz))
  g
}

limit.matrix.to.protein.coding <- function(mat, use.symbols = TRUE) {

    use.biomaRt <- TRUE
    exclude.mt <- FALSE
    # Need to do follow this
    # https://support.bioconductor.org/p/p132709/#p133562
    # to resolve an error in _filter
    #Install previous version of dplyr (0.8) 
    # devtools::install_github("hadley/dplyr@v0.8.0")

    #Install previous version of dbplyr (1.3)
    #devtools::install_url(https://cran.r-project.org/src/contrib/Archive/dbplyr/dbplyr_1.3.0.tar.gz)

    #Install AnnotationHub for Bioconductor 3.10 
    #BiocManager::install("AnnotationHub", version = "3.10", dependencies = TRUE)
    if(use.biomaRt) {
        suppressPackageStartupMessages(p_load(biomaRt))
        ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
        pc.tbl <- getBM(attributes = c("ensembl_gene_id","external_gene_name","description","chromosome_name"),
                        filters = 'biotype', values = c('protein_coding'), mart = ensembl)
        if(exclude.mt) {
            pc.tbl <- subset(pc.tbl, chrosome_name != "MT")
        }
        keys <- as.character(pc.tbl$external_gene_name)
        if(!use.symbols) { keys <- as.character(pc.tbl$ensembl_gene_id) }
        mat <- mat[rownames(mat) %in% keys,]
        mat <- sweep(mat, 2, colSums(mat), `/`) * 10^6
        return(mat)
    }
    if(exclude.mt) { stop("Have not implemented exclude.mt for annotationHub -- but is easy\n") }
    suppressPackageStartupMessages(p_load(AnnotationHub))
    suppressPackageStartupMessages(p_load(ensembldb))
    ah <- AnnotationHub()
    flag <- (ah$species == "Homo sapiens") & (ah$genome == "GRCh38") & (ah$dataprovider == "Ensembl") & (ah$rdataclass == "EnsDb")
    ah2 <- ah[flag, ]
    ## as.data.frame(mcols(ah2))[1:10,c("title"),drop=FALSE]
    edb <- ah2[["AH73881"]]

    ## keytypes(edb)
    ## columns(edb)
    keys <- keys(edb, "GENENAME")
    columns <- c("GENEID", "ENTREZID", "GENEBIOTYPE")
    tbl <- ensembldb::select(edb, keys, columns, keytype = "GENENAME")
    pc.tbl <- subset(tbl, GENEBIOTYPE == "protein_coding")

    keys <- as.character(pc.tbl$GENENAME)
    if(!use.symbols) { keys <- as.character(pc.tbl$GENEID) }
    mat <- mat[rownames(mat) %in% keys,]
    mat <- sweep(mat, 2, colSums(mat), `/`) * 10^6
    mat
}


## Create in silico admixtures from the purified profiles
## sample.ratios is a data frame with columns sample.col (indicating the name of a sample
## admixture), cell.type.col (indicating a cell type within the sample), and fraction.col
## (indicating the fraction of the corresponding cell type within the sample). Note that
## the cell types indicated in cell.type.col should exist in mat, which provides the
## profiles that are mixed together.
create.in.silico.admixtures <- function(mat, sample.ratios,
                                        sample.col = "sample", cell.type.col = "cell.type",
                                        fraction.col = "actual") {
    insilico.admixtures <-
        dlply(sample.ratios,
              .variables = sample.col,
              .fun = function(df) {
                  cols <- as.character(df[, cell.type.col])
                  fracs <- as.numeric(df[, fraction.col])
                  tmp.mat <- mat[, cols]
                  ret <- as.matrix(tmp.mat) %*% fracs
                  colnames(ret) <- df[1, sample.col]
                  ret
              })
    insilico.admixtures <- do.call(cbind, insilico.admixtures)
    insilico.admixtures
}

