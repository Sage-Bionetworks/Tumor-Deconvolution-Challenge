fit.proportions <- function(mat, top, prefix, plot.pdf = TRUE, cnts = NULL) {
    glist <- list()

    flag <- unlist(apply(mat, 1, function(row) any(is.na(row))))
    mat <- mat[!flag, ]

    flag <- unlist(apply(mat, 1, function(row) !(all(row == 0))))
    mat <- mat[flag, ]

    flag <- unlist(apply(mat, 2, function(col) !(all(col == col[1]))))
    mat <- mat[, flag]

    ## Calculate PCA of covariance matrix _before_ normalizing
    if(!(all(rowSums(mat) == 1))) {
        pca <- prcomp(cov(mat, use = "pairwise.complete.obs"), center = FALSE, scale = FALSE)
        cat(paste("Eigenvalues of covariance matrix before normalization:\n"))
        print(pca$sdev)
    }

    pca <- prcomp(cov(log(eps+mat), use = "pairwise.complete.obs"), center = FALSE, scale = FALSE)
    cat(paste("Eigenvalues of covariance of log matrix before normalization:\n"))
    print(pca$sdev)

    
    ## Adjust proportions to sum to 100
    sums <- unlist(apply(mat, 1, function(row) sum(row)))
    for(i in 1:nrow(mat)) {
        mat[i, ] <- mat[i, ] / sums[i]
    }

    ## Calculate PCA of covariance matrix _after_ normalizing
    pca <- prcomp(cov(mat, use = "pairwise.complete.obs"), center = FALSE, scale = FALSE)
    cat(paste("Eigenvalues of covariance matrix after normalization:\n"))
    print(pca$sdev)
    
    means <- colMeans(mat * 100)
    if(!is.null(cnts)) {
        ## Compute weighted mean
        print(rownames(mat))
        print(cnts)
        res <- cov.wt(x = mat * 100, wt = cnts[rownames(mat)]/sum(cnts[rownames(mat)]), cor = TRUE, center = TRUE,
                      method = "unbiased")
        means <- res$center
    }
    o <- order(means, decreasing = FALSE)
    means <- means[o]
    mat <- mat[, o]
    cors <- cor(mat * 100)
    covs <- cov(mat * 100)
    if(!is.null(cnts)) {
        ## Compute weighted covariance and correlation
        res <- cov.wt(x = mat * 100, wt = cnts[rownames(mat)]/sum(cnts[rownames(mat)]), cor = TRUE, center = TRUE,
                      method = "unbiased")
        cors <- res$cor
        covs <- res$cov
    }

    eps <- 10^-5
    ## Adjust proportions to sum to 100
    eps.mat <- eps + mat
##    eps.mat <- eps.mat[, !(colnames(eps.mat) == "MastCell")]
    exclude.zero.rows <- FALSE
    flag <- unlist(apply(mat, 1, function(row) any(row == 0)))
    print(table(flag))
    if(exclude.zero.rows) {
        eps.mat <- mat[!flag, ]
    }
    eps.sums <- unlist(apply(eps.mat, 1, function(row) sum(row)))
    for(i in 1:nrow(eps.mat)) {
        eps.mat[i, ] <- eps.mat[i, ] / eps.sums[i]
    }
    log.means <- colMeans(log(eps.mat * 100))
    log.covs <- cov(log(eps.mat * 100))
    if(!is.null(cnts)) {
        res <- cov.wt(x = log(eps.mat * 100), wt = cnts[rownames(mat)]/sum(cnts[rownames(mat)]), cor = TRUE, center = TRUE,
                      method = "unbiased")
        log.means <- res$center
        log.covs <- res$cov
    }
    ##    print(log.covs)
    ##    print(exp(log.means))
    ##    print(exp(log.covs))
    
    mvln.samples <- exp(mvrnorm(n = 1000, mu = log.means, Sigma = log.covs))
    colnames(mvln.samples) <- colnames(mat)
    mvln.cors <- cor(mvln.samples)
    mvln.covs <- cor(mvln.samples)
##    pca <- prcomp(mvln.covs)
##    print(pca$sdev)

    ## Proportions can be negative, but they sum to 100.
    ## This must be because cors is degenerate.
    ## Generate n random samples from the multivariate normal
    mvn.samples <- mvrnorm(n = 1000, mu = means, Sigma = covs)
    colnames(mvn.samples) <- colnames(mat)
    mvn.cors <- cor(mvn.samples)
##    pca <- prcomp(mvn.cors)
##    print(pca$sdev)
    
    ## install.packages("dirichlet", repos=c("http://R-Forge.R-project.org"))
    
    ## Fit Dirichlet
    d.fit <- fit.dirichlet(mat, type="mm")
    
    ## Generate n random samples from the Dirichlet distribution
    d.samples <- rdirichlet(n = 1000, d.fit$p) * 100
    colnames(d.samples) <- colnames(mat)
    d.cors <- cor(d.samples)
    
    ## Fit Generalized Dirichlet
    gd.fit <- fit.genDirichlet(mat, type="mm")
    
    ## Generate n random samples from the Generalized Dirichlet distribution
    ## For some reason, p must be in decreasing order
    o <- order(gd.fit$p, decreasing=TRUE)
    gd.samples <- rGenDirichlet(n = 1000, p=gd.fit$p[o], k=gd.fit$k[o]) * 100
    colnames(gd.samples) <- colnames(mat)[o]
    gd.samples <- gd.samples[, colnames(mat)] 
    gd.cors <- cor(gd.samples)
    
    cell.types <- colnames(mat)
    names(cell.types) <- cell.types
    
    bounds <- apply(mat * 100, 2, function(vec) calc.CI(vec, 0.95))
    ## gd.plt <- plot.means.and.bounds(colMeans(gd.samples), diag(cov(gd.samples)), main = "Dirichlet")
    emp.ret <- plot.means.and.bounds(colMeans(mat * 100),
                                    lower.bounds = bounds[1, ],
                                    upper.bounds = bounds[2, ],
                                    main = "Empirical")
    emp.plt <- emp.ret$g
    emp.plt <- emp.plt + xlab("")
    emp.plt <- emp.plt + theme(axis.text.x=element_blank())
    
    bounds <- apply(gd.samples, 2, function(vec) calc.CI(vec, 0.95))
    ## gd.plt <- plot.means.and.bounds(colMeans(gd.samples), diag(cov(gd.samples)), main = "Dirichlet")
    gd.ret <- plot.means.and.bounds(colMeans(gd.samples),
                                    lower.bounds = bounds[1, ],
                                    upper.bounds = bounds[2, ],
                                    main = "Generalized Dirichlet")
    gd.plt <- gd.ret$g
    gd.plt <- gd.plt + xlab("")
    gd.plt <- gd.plt + theme(axis.text.x=element_blank())
    
    bounds <- apply(d.samples, 2, function(vec) calc.CI(vec, 0.95))
    d.ret <- plot.means.and.bounds(colMeans(d.samples),
                                   lower.bounds = bounds[1, ],
                                   upper.bounds = bounds[2, ],
                                   main = "Dirichlet")
    ## d.plt <- plot.means.and.bounds(colMeans(d.samples), diag(cov(d.samples)), main = "Dirichlet")
    d.plt <- d.ret$g
    d.plt <- d.plt + xlab("")
    d.plt <- d.plt + theme(axis.text.x=element_blank())
    
    bounds <- apply(mvn.samples, 2, function(vec) calc.CI(vec, 0.95))
    mvn.ret <- plot.means.and.bounds(colMeans(mvn.samples),
                                     lower.bounds = bounds[1, ],
                                     upper.bounds = bounds[2, ],
                                     main = "Multivariate Normal")
    ## mvn.plt <- plot.means.and.bounds(means, diag(covs), main = "Multivariate Normal")
    mvn.plt <- mvn.ret$g
    
    bounds <- apply(mvln.samples, 2, function(vec) calc.CI(vec, 0.95))
    mvln.ret <- plot.means.and.bounds(colMeans(mvln.samples),
                                     lower.bounds = bounds[1, ],
                                     upper.bounds = bounds[2, ],
                                     main = "Multivariate Log Normal")

    emp.ret$df$distribution <- "Empirical"
    d.ret$df$distribution <- "Dirichlet"
    gd.ret$df$distribution <- "Generalized Dirichlet"
    mvn.ret$df$distribution <- "Normal"
    mvln.ret$df$distribution <- "Log Normal"
    all.df <- rbind(emp.ret$df, d.ret$df, gd.ret$df, mvn.ret$df)
    lvls <- c("Empirical", "Dirichlet", "Generalized Dirichlet", "Normal")
    use.mvln <- TRUE
    if(use.mvln) {
        all.df <- rbind(all.df, mvln.ret$df)
        lvls <- c(lvls, "Log Normal")
    }
    all.df$distribution <- factor(all.df$distribution, levels = lvls)
    
    if(plot.pdf) {
        pdf(paste0(prefix, "-distribution-fits.pdf"), onefile=FALSE)
    }
##    glist[["fit"]] <- do.call("grid.arrange", c(list(emp.plt, d.plt, gd.plt, mvn.plt), "top" = top))
    ##    glist[["fit"]] <- arrangeGrob(grobs = list(emp.plt, d.plt, gd.plt, mvn.plt), "top" = top, nrow = 4, ncol = 1)

    g <- ggplot(data = all.df, aes(x = cell.type, y = mean))
    g <- g + geom_bar(stat="identity", color="black")
    g <- g + facet_wrap( ~ distribution)
    g <- g + geom_errorbar(aes(ymin = lower, ymax = upper), width = .2)
    g <- g + geom_text(aes(x = cell.type, y = 100, label = round(mean, digits=2)))
    g <- g + xlab("Cell Type")
    g <- g + ylab("Proportion")
    g <- g + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    lower.bounds <- all.df$lower
    upper.bounds <- all.df$upper
    g <- g + ylim(c(min(lower.bounds-0.5, 0), max(upper.bounds+0.5, 110)))
    g <- g + ggtitle(top)

    glist[["fit"]] <- g
    ## plot(g)
    if(plot.pdf) {
        plot(glist[["fit"]])
        d <- dev.off()
    }

    if(nrow(cors) > 2) {
        cor.vec <- cors[upper.tri(cors, diag=FALSE)]
        lst <- list()
        ## lst <- list("dirichlet" = d.cors, "generalized dirichlet" = gd.cors, "multivariate normal" = mvn.cors)
        lst[["multivariate normal"]] <- mvn.cors
        lst[["multivariate log-normal"]] <- mvln.cors
        df <- ldply(lst,
                    .fun = function(l) {
                        vec <- l[upper.tri(l, diag=FALSE)]
                        cr <- format(cor(cor.vec, vec, method = "pearson"), digits = 2)
                        data.frame(correlation = cr)
                    })
        colnames(df) <- c("distribution", "correlation of cell type correlations\n(sample vs distribution)")
        df <- df[order(as.numeric(as.character(df[,2]))), ]
        if(plot.pdf) {
            pdf(paste0(prefix, "-cor-of-cors.pdf"))
            grid.table(df, rows = NULL)
            d <- dev.off()
        }
        
        plts <- llply(names(lst),
                      .fun = function(nm) {
                          l <- lst[[nm]]
                          vec <- l[upper.tri(l, diag=FALSE)]
                          g <- plot.correlation(cor.vec, vec, display.r2 = TRUE)
                          g <- g + xlab("Empirical cell/cell correlations")
                          g <- g + ylab(paste0(nm, "\ncell/cell correlations"))
                          g
                      })

        if(plot.pdf) {
            pdf(paste0(prefix, "-cor-of-cors-plots.pdf"), onefile=FALSE)
        }
        ## glist[["cor"]] <- do.call("grid.arrange", c(plts, "top" = top))
        glist[["cor"]] <- arrangeGrob(grobs = plts, "top" = top)
        ## arrangeGrob(grobs = plts, "top" = top, nrow = 2, ncol = 2)
        if(plot.pdf) {
            plot(glist[["cor"]])
            d <- dev.off()
        }
    }
    glist
}

