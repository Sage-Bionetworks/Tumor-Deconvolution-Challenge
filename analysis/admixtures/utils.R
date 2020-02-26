unif.broken.stick.proportion <- function(n, min.prop = 0, max.prop = 1, step = 0.01) {
    sq <- seq(from = min.prop, to = max.prop - step, by = step)
    sam <- sample(sq, size = n-1, replace = FALSE)
    sam <- c(0,sort(sam), max.prop)
    ret <- unlist(lapply(1:(length(sam)-1), function(i) sam[i+1]-sam[i]))
    eps <- 10^-4
    if(abs(sum(ret) - max.prop) > eps) { stop(paste0("unif.broken.stick.proportion sum is not 1, but = ", sum(ret), ": ", paste(ret, collapse=","), "\n")) }
    ret
}

generate.random.uniform.admixtures <- function(populations, n, min.prop = 0, max.prop = 1, lbs = NULL, ubs = NULL, step = 0.01) {

    ## Generate the random admixtures
    mat <- ldply(1:n, .fun = function(i) unif.broken.stick.proportion(length(populations), min.prop = min.prop, max.prop = max.prop, step = step))
    if(is.null(lbs)) { lbs <- rep(0, length(populations)) }
    if(is.null(ubs)) { ubs <- rep(1, length(populations)) }    
    
    flag <- unlist(apply(mat, 1, function(row) any(row < lbs) | any(row > ubs)))
    while(any(flag)) {
        num <- length(which(flag))
        mat[flag,] <- ldply(1:num, .fun = function(i) unif.broken.stick.proportion(length(populations), min.prop = min.prop, max.prop = max.prop, step = step))
        flag <- unlist(apply(mat, 1, function(row) any(row < lbs) | any(row > ubs)))
    }
    ## Generate the admixture with only tumor content
##    cancer.only.admixture <- rep(0, length(populations))
##    cancer.only.admixture[populations == tumor.type] <- 1
##    mat <- rbind(mat, cancer.only.admixture)
    mat <- t(mat)
    rownames(mat) <- populations
    colnames(mat) <- NULL
    mat
}

generate.random.dirichlet.admixtures <- function(proportions, n, tumor.type, min.prop = 0) {

    ## Generate the random admixtures
    mat <- rdirichlet(n=n-1, alpha=proportions*5)
    flag <- unlist(apply(mat, 1, function(row) any(row < min.prop)))
    while(any(flag)) {
        mat[flag,] <- rdirichlet(n=length(which(flag)), alpha=proportions*5)
        flag <- unlist(apply(mat, 1, function(row) any(row < min.prop)))
    }
    
    ## Generate the admixture with only tumor content
    cancer.only.admixture <- rep(0, length(proportions))
    cancer.only.admixture[names(proportions) == tumor.type] <- 1

    mat <- rbind(mat, cancer.only.admixture)
    mat <- t(mat)
    rownames(mat) <- names(proportions)
    colnames(mat) <- NULL
    mat
}

plot.admixture.correlations <- function(admixtures, ...) {

    fc <- as.matrix(cor(admixtures, method = "spearman"))
    
    mar <- c(1, 1, 1, 1)
    ##    corrplot(fc, method = "ellipse", type = "upper", order = "hclust", tl.cex = 0.6, mar = mar, ...)
    corrplot(fc, method = "ellipse", type = "upper", order = "original", tl.cex = 0.6, mar = mar, ...)
}

plot.frequencies <- function(admixtures) {
    m <- melt(admixtures)
    colnames(m) <- c("Population", "Admixture", "Proportion")
    g <- ggplot()
    m$Population <- factor(m$Population, levels = rownames(admixtures))
    ##    g <- g + geom_point(data = m, aes(x = Population, y = Proportion))
    ## g <- g + geom_beeswarm(data = m, aes(x = Population, y = log(Proportion)))
    ## g <- g + geom_violin(data = m, aes(x = Population, y = Proportion))
##    g <- g + geom_boxplot(data = m, aes(x = Population, y = log(Proportion)))
##    g <- g + geom_beeswarm(data = m, aes(x = Population, y = Proportion))    
    g <- g + geom_boxplot(data = m, aes(x = Population, y = Proportion))    
    g <- g + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    g
}

sample.dist <- function(n, dist.name, params, step = 0.01) {
    ret <- NA
    if(dist.name == "dir") {
        ret <- rdirichlet(n, params)
    } else if(dist.name == "unif") {
        cols <- colnames(params)
        names(cols) <- cols
        ret <- ldply(cols,
                     .fun = function(col) {
                         runif(n,
                               min = as.numeric(params["min", col]),
                               max = as.numeric(params["max", col]))
                     })
        pops <- as.character(ret$.id)
        ret <- t(ret[, colnames(ret) != ".id"])
        colnames(ret) <- pops
    } else if(dist.name == "unif-step") {
        cols <- colnames(params)
        names(cols) <- cols
        ret <- ldply(cols,
                     .fun = function(col) {
                         sample(seq(from = as.numeric(params["min", col]),
                                    to = as.numeric(params["max", col]), by = step), size=n, replace=TRUE)
                     })
        pops <- as.character(ret$.id)
        ret <- t(ret[, colnames(ret) != ".id"])
        colnames(ret) <- pops
    } else {
        stop(paste0("Unknown distribution: ", dist.name, "\n"))
    }
    return(ret)
}

sample.hierarchical.model_ <- function(n, hierarchical.model) {
    df <- data.frame(root = rep(1, n))
    l <- hierarchical.model
    while(length(l) > 0) {
        for(i in 1:length(l)) {
            if(!l[[i]]$parent %in% colnames(df)) { next }
            mat <- sample.dist(n, l[[i]]$dist, l[[i]]$params)
            mat <- mat * df[, l[[i]]$parent]
            colnames(mat) <- names(l[[i]]$params)
            if(is.null(names(l[[i]]$params))) {
                colnames(mat) <- colnames(l[[i]]$params)
            }
            df <- cbind(df, mat)
            l[[i]] <- NULL
            break
        }
    }
    df
}

sample.flat.unif.model <- function(n, flat.model, min.prop = 0.01, constraint.col = "tumor.fraction") {
    factor <- 10000
    df <- sample.dist(n*factor, "unif-step", flat.model)
    df[, constraint.col] <- 1 - rowSums(df[, !(colnames(df) == constraint.col), drop=F])
    flag <- unlist(apply(df, 1, function(row) all(row >= min.prop)))
    flag <- flag & df[, constraint.col] >= flat.model["min", constraint.col] &
        df[, constraint.col] < flat.model["max", constraint.col]
    while(length(which(flag)) < n) {
        print(length(which(flag)))
        tmp <- sample.dist(length(which(!flag)*factor), "unif-step", flat.model)
        tmp[, constraint.col] <- 1 - rowSums(tmp[, !(colnames(tmp) == constraint.col), drop=F])
        tmp <- tmp[order(tmp[, constraint.col], decreasing=TRUE),]
        flag <- unlist(apply(df, 1, function(row) all(row >= min.prop)))
        flag <- flag & df[, constraint.col] >= flat.model["min", constraint.col] &
            df[, constraint.col] < flat.model["max", constraint.col]
        df[!flag,] <- tmp[1:length(which(!flag)),,drop=F]
        flag <- unlist(apply(df, 1, function(row) all(row >= min.prop)))
        flag <- flag & df[, constraint.col] >= flat.model["min", constraint.col] &
            df[, constraint.col] < flat.model["max", constraint.col]
    }
    df <- df[flag,]
    df <- df[1:n,]
}
    
sample.hierarchical.model <- function(n, hierarchical.model, pops, min.prop = 0.01) {
    df <- sample.hierarchical.model_(n, hierarchical.model)
    df <- df[, colnames(df) %in% pops]
    df <- df / rowSums(df)

    flag <- unlist(apply(df, 1, function(row) any(row < min.prop)))
    while(any(flag)) {
        print(length(which(flag)))
        tmp <- sample.hierarchical.model_(length(which(flag)), hierarchical.model)
        tmp <- tmp[, colnames(tmp) %in% pops]
        tmp <- tmp / rowSums(tmp)
        df[flag,] <- tmp
        flag <- unlist(apply(df, 1, function(row) any(row < min.prop)))
    }
    t(df)
}

flatten.hierarchical.unif.model <- function(hierarchical.model, pops) {
    df <- data.frame(root = c(min=1, max=1))
    l <- hierarchical.model
    while(length(l) > 0) {
        for(i in 1:length(l)) {
            if(!l[[i]]$parent %in% colnames(df)) { next }
            mat <- l[[i]]$params
            mat["min",] <- mat["min",] * df["min", l[[i]]$parent]
            mat["max",] <- mat["max",] * df["max", l[[i]]$parent]
            colnames(mat) <- colnames(l[[i]]$params)
            df <- cbind(df, mat)
            l[[i]] <- NULL
            break
        }
    }
    df[, colnames(df) %in% pops]
}
