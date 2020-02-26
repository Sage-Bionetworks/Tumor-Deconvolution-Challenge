zero.out.populations.old <- function(tbl) {

  pops.to.zero.out <- list(
    "Dendritic_cells" = "Monocytes",
    "Monocytes" = "Macrophages",
    "Macrophages" = "Dendritic_cells",
    "Naive_CD4_T_cells" = "Memory_CD4_T_cells",
    "Memory_CD4_T_cells" = "Naive_CD4_T_cells",
    "Naive_CD8_T_cells" = "Memory_CD8_T_cells",
    "Memory_CD8_T_cells" = "Naive_CD8_T_cells",
    "Naive_B_cells" = "Memory_B_cells",
    "Memory_B_cells" = "Naive_B_cells",
    "Tregs" = "Naive_CD4_T_cells",
    "Naive_CD4_T_cells" = "Tregs"
  )

  next.pop <- 1
  for(i in 1:nrow(tbl)) {
    val <- tbl[i, names(pops.to.zero.out)[next.pop]]
    ## Don't zero out very small values
    if(val < 0.02) { next }
    tbl[i, names(pops.to.zero.out)[next.pop]] <- 0
    tbl[i, pops.to.zero.out[[next.pop]]] <- tbl[i, pops.to.zero.out[[next.pop]]] + val
    next.pop <- next.pop + 1
    if(next.pop > length(pops.to.zero.out)) { next.pop <- 1 }
  }

  tbl
}

zero.out.populations <- function(tbl, pops.to.zero.out) {

  next.pop <- 1
  remaining.indices <- 1:nrow(tbl)
  while(length(remaining.indices) > 0) {
    
  }
  for(i in 1:nrow(tbl)) {
    src.val <- tbl[i, names(pops.to.zero.out)[next.pop]]
    target.val <- tbl[i, names(pops.to.zero.out)[next.pop]]    
    ## Don't zero out very small values
    if(val < 0.02) { next }
    tbl[i, names(pops.to.zero.out)[next.pop]] <- 0
    tbl[i, pops.to.zero.out[[next.pop]]] <- tbl[i, pops.to.zero.out[[next.pop]]] + val
    next.pop <- next.pop + 1
    if(next.pop > length(pops.to.zero.out)) { next.pop <- 1 }
  }

  tbl
}

## Use greedy approach to find m extremal pts / columns of mat,
## using the approach described in Appendix C of White and Shalloway (2009) Phys Rev E
find.extremal.points.shalloway <- function(mat, m) {
  if(ncol(mat) <= 2) { return(mat) }

  ## Find the two mutually most distant points
##  d <- as.matrix(dist(t(mat)))
##  max.indices <- which(d == max(d), arr.ind=TRUE)
##  r.indices <- as.numeric(max.indices[1,])
  r.indices <- c()
  max.dst <- 0
  nr <- nrow(mat)
  nc <- ncol(mat)  
  for(i in 1:(nc-1)) {
    vi <- mat[,i]
    for(j in 2:i) {
      vj <- mat[,j]
      dst <- sum((vi-vj)^2)
      if(dst > max.dst) { max.dst <- dst; r.indices <- c(i,j); }
    }
  }
  q <- 2

  ## Define projection matrix
  Pq <- diag(nr)
  r1 <- r.indices[1]
  qp <- r.indices[q]
  diff.vec <- as.numeric(mat[, qp]) - as.numeric(mat[, r1])
##  denom <- sqrt(sum(diff.vec^2))^2
  denom <- sum(diff.vec^2)
  for(i in 1:nr) {
    for(j in 1:nr) {
      Pq[i,j] <- Pq[i,j] - ( diff.vec[i] * diff.vec[j] ) / denom
    }
  }

  while(q < m) {

    ## Select the q+1 item as that which maximizes the
    ## projected distance
    dsts <- unlist(apply(mat, 2, function(col) { diff.vec <- as.numeric(col) - as.numeric(mat[,r1]); vec <- Pq %*% diff.vec; sum(vec^2) }))
    print(r.indices)		   
##    dsts[r.indices] <- 0
    q <- q + 1
    r.indices[q] <- which.max(dsts)
    ## Update the projection matrix
    qp <- r.indices[q]
    diff.vec <- as.numeric(mat[, qp]) - as.numeric(mat[, r1])
##    denom <- sqrt(sum(diff.vec^2))^2
    denom <- sum(diff.vec^2)
    for(i in 1:nr) {
      for(j in 1:nr) {
        Pq[i,j] <- Pq[i,j] - ( diff.vec[i] * diff.vec[j] ) / denom
      }
    }
  }
  mat[, r.indices]
}

distance.squared.from.pt.to.hyperplane <- function(hyperplane, pt) {
  w <- hyperplane[1:(length(hyperplane)-1)]
  b <- hyperplane[length(hyperplane)]
  num <- sum(w*pt) + b
  denom <- sum(w*w)
  num / denom
}

library(OjaNP) # for hyperplane
find.extremal.points.hyp <- function(mat, m) {
  if(ncol(mat) <= 2) { return(mat) }

  ## Find the two mutually most distant points
##  d <- as.matrix(dist(t(mat)))
##  max.indices <- which(d == max(d), arr.ind=TRUE)
##  r.indices <- as.numeric(max.indices[1,])
  r.indices <- c()
  max.dst <- 0
  nr <- nrow(mat)
  nc <- ncol(mat)  
  for(i in 1:(nc-1)) {
    vi <- mat[,i]
    for(j in 2:i) {
      vj <- mat[,j]
      dst <- sum((vi-vj)^2)
      if(dst > max.dst) { max.dst <- dst; r.indices <- c(i,j); }
    }
  }
  q <- 2

  while(q < m) {

    ## Select the q+1 item that maximizes the distance to the hyperplane
    ## defined by the q items
    hyp <- hyperplane(t(mat[,r.indices]))

    dsts <- unlist(apply(mat, 2,
                   function(col) { distance.squared.from.pt.to.hyperplane(hyp, col) }))
    q <- q + 1
    r.indices[q] <- which.max(dsts)
  }
  mat[, r.indices]
}

find.extremal.points <- function(mat, m) {
  if(ncol(mat) <= 2) { return(mat) }

  ## Find the two mutually most distant points
##  d <- as.matrix(dist(t(mat)))
##  max.indices <- which(d == max(d), arr.ind=TRUE)
##  r.indices <- as.numeric(max.indices[1,])
  r.indices <- c()
  max.dst <- 0
  nr <- nrow(mat)
  nc <- ncol(mat)  
  for(i in 1:(nc-1)) {
    vi <- mat[,i]
    for(j in 2:i) {
      vj <- mat[,j]
      dst <- sum((vi-vj)^2)
      if(dst > max.dst) { max.dst <- dst; r.indices <- c(i,j); }
    }
  }
  q <- 2

  while(q < m) {

    ## Select the q+1 item that maximizes the minimum distance to any
    ## of the previously selected q items

    min.dsts <-
      unlist(llply(1:ncol(mat), .parallel = TRUE,
                   .fun = function(i) {
		     col <- as.numeric(mat[,i])
		     min.dst <-
		       min(unlist(lapply(r.indices,
		                  function(j) {
                                    vj <- as.numeric(mat[,j])
				    dist <- sum((col - vj)^2)
				    dist
                                  })))
		     min.dst
		   }))
    q <- q + 1
    r.indices[q] <- which.max(min.dsts)
  }
  if(any(duplicated(r.indices))) { stop("Unexpectedly got duplicated indices\n") }
  mat[, r.indices]
}

chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))

find.extremal.points.across.matrices <- function(mats, m) {

  ## Find the two mutually most distant points, where each point
  ## resides in a different matrix.
  ## NO: For efficiency, assume that the first and last matrix in the list
  ## NO: are most different and select points from them.
  ## For efficiency, take our first two points from the two smallest matrices.
  ## Since the complexity of the below n nested for loops is nc1 x nc2,
  ## where nci is the # of columns of the ith matrix.
  ## In the common case (where m >= # of mats), we are going to select
  ## points from these matrices anyway. This ordering will affect our
  ## greedy choice, but seems a fair tradeoff.
  num.cols <- unlist(llply(mats, ncol))
  o <- order(num.cols, decreasing=FALSE)
  mats <- mats[o]
  
  nr <- nrow(mats[[1]])
  mat.indices <- c()
  r.indices <- c()
  max.dst <- 0
  for(mat.i in c(1)) {
##    print(mat.i)
    nc.i <- ncol(mats[[mat.i]])
    indices <- 1:nc.i
    names(indices) <- indices
    ret <-
      ldply(indices, .parallel = TRUE,
            .fun = function(i) {           
                     vi <- mats[[mat.i]][,i]
		     mat.indices <- c(2)
		     names(mat.indices) <- mat.indices
                     rt <-
		       ldply(mat.indices, .parallel = FALSE,
		             .fun = function(mat.j) {
                                      print(c(i, mat.i, mat.j))
                                      nc.j <- ncol(mats[[mat.j]])
                                      dsts <- unlist(llply(1:nc.j, .parallel = FALSE,
                                                     .fun = function(j) {
                                                              vj <- mats[[mat.j]][,j]
                                                              dst <- sum((vi-vj)^2)
                                                              dst
                                                            }))
                                      df <- data.frame(dist = dsts, j = 1:nc.j, mat.j = mat.j)
				      df <- df[which.max(df$dist), ,drop=F]
                                    })
                     colnames(rt)[1] <- "mat.i"
		     rt <- rt[which.max(rt$dist), ,drop=F]
		     rt
		   })
    colnames(ret)[1] <- c("i")
    ret <- cbind(rep(mat.i,nrow(ret)), ret)
    colnames(ret)[1] <- c("mat.i")
    dst <- max(ret$dist)
    if(dst > max.dst) {
      indx <- which.max(ret$dist)
      mat.indices <- c(as.numeric(ret[indx, "mat.i"]),
                       as.numeric(ret[indx, "mat.j"]))
      r.indices <- c(as.numeric(ret[indx, "i"]),
                       as.numeric(ret[indx, "j"]))
		       
    }
  }

  q <- 2

  cat("Done selecting firste two points\n")

  ## Select points from each remaining matrix, recycling as necessary
  all.matrix.indices <- 1:length(mats)
  candidate.matrices <- all.matrix.indices[!(all.matrix.indices %in% mat.indices)]
  
  num.needed <- m - q
  if( num.needed > length(candidate.matrices) ) {
    candidate.matrices <- c(candidate.matrices, rep_len(all.matrix.indices, num.needed - length(candidate.matrices)))
  }
  candidate.matrices <- candidate.matrices[1:num.needed]

  ret.mat <- Reduce("cbind", lapply(1:length(r.indices), function(i) mats[[mat.indices[i]]][, r.indices[i]]))

  cat("About to iterate to find remaining points\n")

  for(next.matrix in candidate.matrices) {

    q <- q + 1

    ## Select the q+1 item that maximizes the minimum distance to any
    ## of the previously selected q items (in ret.mat)

    remaining.mats <- 1:length(mats)
    remaining.mats <- remaining.mats[!(remaining.mats %in% mat.indices)]
    remaining.mats <- q
    remaining.mats <- next.matrix
    names(remaining.mats) <- remaining.mats
    cat(paste0("Find extremal point ", q, " in matrix with # cols = ", ncol(mats[[remaining.mats[[1]]]]), "\n"))
    cat(paste0("ret.mat = ", nrow(ret.mat) , " x ", ncol(ret.mat), "\n"))
    dsts <- ldply(remaining.mats, .parallel = FALSE,
                  .fun = function(mat.i) {
		           cols <- 1:ncol(mats[[mat.i]])
			   cat("Start llply\n")
			   chunks <- chunk2(cols, num.cores-1)
                           ret <- unlist(llply(chunks, .parallel = TRUE,
			                .fun = function(is) {
					         unlist(llply(is,
						   .fun = function(i) {
                                                            col <- as.numeric(mats[[mat.i]][, i])
           		                                    min.dst <-
 		                                              min(unlist(apply(ret.mat, 2,
		                                                         function(vj) {
   				                                           dist <- sum((col - vj)^2)
				                                           dist
                                                                         })))
						            min.dst
							  }))
                                               }))
	                   cat("Done with llply\n")
	                   indx <- which.max(ret)
			   cat("Done with index\n")
			   dist <- ret[indx]
			   df <- data.frame(index = indx, dist = dist)
			   cat("Done with df\n")
			   df
                        })
    colnames(dsts) <- c("mat.index", "index", "dist")

    indx <- which.max(dsts$dist)[1]
    r.index <- as.numeric(dsts[indx, "index"])
    mat.index <- as.numeric(dsts[indx, "mat.index"])    
    col <- as.numeric(mats[[mat.index]][, r.index])
    r.indices[q] <- r.index
    mat.indices[q] <- mat.index
    ret.mat <- cbind(ret.mat, col)
  }
  tmp <- paste(mat.indices, r.indices)
  o <- order(mat.indices, r.indices)
  tmp <- tmp[o]
  if(any(duplicated(tmp))) { stop("Was not expecting duplicated indices\n") }
  print(tmp)
  colnames(ret.mat) <- NULL
  ret.mat
}

old.find.extremal.points.across.matrices <- function(mats, m) {

  ## Find the two mutually most distant points, where each point
  ## resides in a different matrix.
  ## NO: For efficiency, assume that the first and last matrix in the list
  ## NO: are most different and select points from them.
  ## For efficiency, take our first two points from the two smallest matrices.
  ## Since the complexity of the below n nested for loops is nc1 x nc2,
  ## where nci is the # of columns of the ith matrix.
  ## In the common case (where m >= # of mats), we are going to select
  ## points from these matrices anyway. This ordering will affect our
  ## greedy choice, but seems a fair tradeoff.
  num.cols <- unlist(llply(mats, ncol))
  o <- order(num.cols, decreasing=FALSE)
  mats <- mats[o]
  
  nr <- nrow(mats[[1]])
  mat.indices <- c()
  r.indices <- c()
  max.dst <- 0
  for(mat.i in c(1)) {
##    print(mat.i)
    nc.i <- ncol(mats[[mat.i]])
    for(i in 1:nc.i) {
      vi <- mats[[mat.i]][,i]
      for(mat.j in c(2)) {
        print(c(i, mat.i, mat.j))
        nc.j <- ncol(mats[[mat.j]])
	if(FALSE) {
          for(j in 1:nc.j) {
            vj <- mats[[mat.j]][,j]
            dst <- sum((vi-vj)^2)
            if(dst > max.dst) { max.dst <- dst; mat.indices <- c(mat.i, mat.j); r.indices <- c(i,j); }
          }
	}
	dsts <- unlist(llply(1:nc.j, .parallel = FALSE,
	              .fun = function(j) {
                               vj <- mats[[mat.j]][,j]
                               dst <- sum((vi-vj)^2)
                               dst
                             }))
        dst <- max(dsts)
	if(dst > max.dst) {
          j <- which.max(dsts)
	  mat.indices <- c(mat.i, mat.j)
	  r.indices <- c(i,j)
	}
      }
    }
  }
  q <- 2

  ## Select points from each remaining matrix, recycling as necessary
  all.matrix.indices <- 1:length(mats)
  candidate.matrices <- all.matrix.indices[!(all.matrix.indices %in% mat.indices)]
  
  num.needed <- m - q
  if( num.needed > length(candidate.matrices) ) {
    candidate.matrices <- c(candidate.matrices, rep_len(all.matrix.indices, num.needed - length(candidate.matrices)))
  }
  candidate.matrices <- candidate.matrices[1:num.needed]

  ret.mat <- Reduce("cbind", lapply(1:length(r.indices), function(i) mats[[mat.indices[i]]][, r.indices[i]]))

  for(next.matrix in candidate.matrices) {

    q <- q + 1

    ## Select the q+1 item that maximizes the minimum distance to any
    ## of the previously selected q items (in ret.mat)

    cat(paste0("Find extremal point ", q, "\n"))
    remaining.mats <- 1:length(mats)
    remaining.mats <- remaining.mats[!(remaining.mats %in% mat.indices)]
    remaining.mats <- q
    remaining.mats <- next.matrix
    names(remaining.mats) <- remaining.mats
    dsts <- ldply(remaining.mats, .parallel = FALSE,
                  .fun = function(mat.i) {
		           cols <- 1:ncol(mats[[mat.i]])
			   names(cols) <- cols
                           ret <- ldply(cols, .parallel = TRUE,
			                .fun = function(i) { 
                                                 col <- mats[[mat.i]][, i]
         		                         min.dst <-
 		                                   min(unlist(apply(ret.mat, 2,
		                                                    function(vj) {
   				                                      dist <- sum((col - vj)^2)
				                                      dist
                                                                    })))
						 min.dst
                                               })
                           colnames(ret) <- c("index", "dist")
			   ret
                        })
    colnames(dsts) <- c("mat.index", "index", "dist")

    indx <- which.max(dsts$dist)[1]
    r.index <- as.numeric(dsts[indx, "index"])
    mat.index <- as.numeric(dsts[indx, "mat.index"])    
    col <- as.numeric(mats[[mat.index]][, r.index])
    r.indices[q] <- r.index
    mat.indices[q] <- mat.index
    ret.mat <- cbind(ret.mat, col)
  }
  tmp <- paste(mat.indices, r.indices)
  o <- order(mat.indices, r.indices)
  tmp <- tmp[o]
  if(any(duplicated(tmp))) { stop("Was not expecting duplicated indices\n") }
##  print(tmp)
  colnames(ret.mat) <- NULL
  ret.mat
}

find.extremal.points.across.matrices.shalloway <- function(mats, m) {

  ## Find the two mutually most distant points, where each point
  ## resides in a different matrix.
  ## NO: For efficiency, assume that the first and last matrix in the list
  ## NO: are most different and select points from them.
  ## For efficiency, take our first two points from the two smallest matrices.
  ## Since the complexity of the below n nested for loops is nc1 x nc2,
  ## where nci is the # of columns of the ith matrix.
  ## In the common case (where m >= # of mats), we are going to select
  ## points from these matrices anyway. This ordering will affect our
  ## greedy choice, but seems a fair tradeoff.
  num.cols <- unlist(llply(mats, ncol))
  o <- order(num.cols, decreasing=FALSE)
  mats <- mats[o]
  
  nr <- nrow(mats[[1]])
  mat.indices <- c()
  r.indices <- c()
  max.dst <- 0
  for(mat.i in c(1)) {
    print(mat.i)
    nc.i <- ncol(mats[[mat.i]])
    for(i in 1:nc.i) {
      vi <- mats[[mat.i]][,i]
      for(mat.j in c(2)) {
        print(c(i, mat.i, mat.j))
        nc.j <- ncol(mats[[mat.j]])
        for(j in 1:nc.j) {
          vj <- mats[[mat.j]][,j]
          dst <- sum((vi-vj)^2)
          if(dst > max.dst) { max.dst <- dst; mat.indices <- c(mat.i, mat.j); r.indices <- c(i,j); }
        }
      }
    }
  }
  q <- 2

  ## Define projection matrix
  Pq <- diag(nr)
  r1 <- r.indices[1]
  m1 <- mat.indices[1]  
  rq <- r.indices[q]
  mq <- mat.indices[q]    
  diff.vec <- as.numeric(mats[[mq]][, rq]) - as.numeric(mats[[m1]][, r1])
##  denom <- sqrt(sum(diff.vec^2))^2
  denom <- sum(diff.vec^2)
  for(i in 1:nr) {
    for(j in 1:nr) {
      Pq[i,j] <- Pq[i,j] - ( diff.vec[i] * diff.vec[j] ) / denom
    }
  }

  ## Select points from each remaining matrix, recycling as necessary
  all.matrix.indices <- 1:length(mats)
  candidate.matrices <- all.matrix.indices[!(all.matrix.indices %in% mat.indices)]
  
  num.needed <- m - q
  if( num.needed > length(candidate.matrices) ) {
    candidate.matrices <- c(candidate.matrices, rep_len(all.matrix.indices, num.needed - length(candidate.matrices)))
  }
  candidate.matrices <- candidate.matrices[1:num.needed]

  for(next.matrix in candidate.matrices) {

    q <- q + 1
    ## Select the q+1 item as that which maximizes the
    ## projected distance
    remaining.mats <- 1:length(mats)
    remaining.mats <- remaining.mats[!(remaining.mats %in% mat.indices)]
    remaining.mats <- q
    remaining.mats <- next.matrix
    names(remaining.mats) <- remaining.mats
    dsts <- ldply(remaining.mats, .parallel = TRUE,
                  .fun = function(mat.i) {
		           cols <- 1:ncol(mats[[mat.i]])
			   names(cols) <- cols
                           ret <- ldply(cols, .parallel = FALSE,
			                .fun = function(i) { 
                                                 col <- mats[[mat.i]][, i]
						 diff.vec <- as.numeric(col) - as.numeric(mats[[m1]][, r1])
						 vec <- Pq %*% diff.vec
						 sum(vec^2)
                                               })
                           colnames(ret) <- c("index", "dist")
			   ret
                        })
    colnames(dsts) <- c("mat.index", "index", "dist")
    flag <- mat.indices == next.matrix
    if(any(flag)) { dsts$dist[r.indices[flag]] <- 0 }

    indx <- which.max(dsts$dist)[1]
    r.indices[q] <- as.numeric(dsts[indx, "index"])
    mat.indices[q] <- as.numeric(dsts[indx, "mat.index"])    

    ## Update the projection matrix
    rq <- r.indices[q]
    mq <- mat.indices[q]    
    diff.vec <- as.numeric(mats[[mq]][, rq]) - as.numeric(mats[[m1]][, r1])
##    denom <- sqrt(sum(diff.vec^2))^2
    denom <- sum(diff.vec^2)
    for(i in 1:nr) {
      for(j in 1:nr) {
        Pq[i,j] <- Pq[i,j] - ( diff.vec[i] * diff.vec[j] ) / denom
      }
    }
  }
  Reduce("cbind", lapply(1:length(r.indices), function(i) mats[[mat.indices[i]]][, r.indices[i]]))
}


find.extremal.points.across.matrices.old <- function(mats) {
  if(length(mats) < 2) { stop("Was expected > 2 matrices\n") }

  m <- length(mats)
  num.cols <- unlist(llply(mats, ncol))
  o <- order(num.cols, decreasing=FALSE)
  mats <- mats[o]
  ## Find the two mutually most distant points, where each point
  ## resides in a different matrix
  nr <- nrow(mats[[1]])
  mat.indices <- c()
  r.indices <- c()
  max.dst <- 0
  for(mat.i in c(1)) {
    print(mat.i)
    nc.i <- ncol(mats[[mat.i]])
    for(i in 1:nc.i) {
      vi <- mats[[mat.i]][,i]
      for(mat.j in c(2)) {
        print(c(i, mat.i, mat.j))
        nc.j <- ncol(mats[[mat.j]])
        for(j in 1:nc.j) {
          vj <- mats[[mat.j]][,j]
          dst <- sum((vi-vj)^2)
          if(dst > max.dst) { max.dst <- dst; mat.indices <- c(mat.i, mat.j); r.indices <- c(i,j); }
        }
      }
    }
  }
  q <- 2

  ## Define projection matrix
  Pq <- diag(nr)
  r1 <- r.indices[1]
  m1 <- mat.indices[1]  
  rq <- r.indices[q]
  mq <- mat.indices[q]    
  diff.vec <- as.numeric(mats[[mq]][, rq]) - as.numeric(mats[[m1]][, r1])
  denom <- sqrt(sum(diff.vec^2))^2
  for(i in 1:nr) {
    for(j in 1:nr) {
      Pq[i,j] <- Pq[i,j] - ( diff.vec[i] * diff.vec[j] ) / denom
    }
  }

  while(q < m) {

    q <- q + 1
    ## Select the q+1 item as that which maximizes the
    ## projected distance
    remaining.mats <- 1:length(mats)
    remaining.mats <- remaining.mats[!(remaining.mats %in% mat.indices)]
    remaining.mats <- q
    names(remaining.mats) <- remaining.mats
    dsts <- ldply(remaining.mats, .parallel = TRUE,
                  .fun = function(mat.i) {
		           cols <- 1:ncol(mats[[mat.i]])
			   names(cols) <- cols
                           ret <- ldply(cols, .parallel = FALSE,
			                .fun = function(i) { 
                                                 col <- mats[[mat.i]][, i]
						 vec <- Pq %*% col
						 sum(vec^2)
                                               })
                           colnames(ret) <- c("index", "dist")
			   ret
                        })
    colnames(dsts) <- c("mat.index", "index", "dist")
    print(head(dsts))
    indx <- which.max(dsts$dist)[1]
    r.indices[q] <- as.numeric(dsts[indx, "index"])
    mat.indices[q] <- as.numeric(dsts[indx, "mat.index"])    

    ## Update the projection matrix
    rq <- r.indices[q]
    mq <- mat.indices[q]    
    diff.vec <- as.numeric(mats[[mq]][, rq]) - as.numeric(mats[[m1]][, r1])
    denom <- sqrt(sum(diff.vec^2))^2
    for(i in 1:nr) {
      for(j in 1:nr) {
        Pq[i,j] <- Pq[i,j] - ( diff.vec[i] * diff.vec[j] ) / denom
      }
    }
  }
  Reduce("cbind", lapply(1:length(r.indices), function(i) mats[[mat.indices[i]]][, r.indices[i]]))
}

unif.broken.stick.proportion <- function(n, min.prop = 0, max.prop = 1, step = 0.01) {
    sq <- seq(from = min.prop, to = max.prop - step, by = step)
    sam <- sample(sq, size = n-1, replace = FALSE)
    sam <- c(0,sort(sam), max.prop)
    ret <- unlist(lapply(1:(length(sam)-1), function(i) sam[i+1]-sam[i]))
    eps <- 10^-4
    if(abs(sum(ret) - max.prop) > eps) { stop(paste0("unif.broken.stick.proportion sum is not 1, but = ", sum(ret), ": ", paste(ret, collapse=","), "\n")) }
    ret
}

generate.random.uniform.admixtures <- function(populations, n, min.prop = 0, max.prop = 1, lbs = NULL, ubs = NULL, step = NULL) {

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

plot.admixture.correlations <- function(admixtures, mar = c(1, 1, 1, 1), ...) {

    fc <- as.matrix(cor(admixtures, method = "spearman"))
    

    ##    corrplot(fc, method = "ellipse", type = "upper", order = "hclust", tl.cex = 0.6, mar = mar, ...)
    corrplot(fc, method = "ellipse", type = "upper", order = "original", tl.cex = 0.6, mar = mar, ...)
##    title(main, line=-1)
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

sample.flat.unif.model.exhaustive <- function(n, flat.model, steps, constraint.col = "tumor.fraction",
                                              min.step = 0.01, num.min.steps = 3) {
  o <- order(as.numeric(flat.model["max",]),decreasing=TRUE)
  flat.model <- flat.model[, o]
  x1 <- seq(from=as.numeric(flat.model["min",1]), to=as.numeric(flat.model["max",1]), by=steps[1])
  x2 <- seq(from=as.numeric(flat.model["min",2]), to=as.numeric(flat.model["max",2]), by=steps[2])
  x1 <- sort(unique(c(x1, seq(from=min.step, to=num.min.steps*min.step, by=min.step))))
  x2 <- sort(unique(c(x2, seq(from=min.step, to=num.min.steps*min.step, by=min.step))))
##  mat <- expand.grid(x1, x2)
  mat <- merge(as.data.frame(x1), as.data.frame(x2))
  colnames(mat) <- colnames(flat.model)[1:ncol(mat)]
  min.remaining <- sum(as.numeric(flat.model["min",(3):ncol(flat.model)]))
  flag <- unlist(apply(mat, 1, function(row) sum(row) <= (1-min.remaining)))
  mat <- mat[flag, ]
  for(i in 3:ncol(flat.model)) {
##  for(i in 3:8) {
    min.remaining <- 0
    if(i < ncol(flat.model)) { min.remaining <- sum(as.numeric(flat.model["min",(i+1):ncol(flat.model)])) }
    min.included <- min(unlist(apply(mat[,1:(i-1)], 1, sum)))
    mx <- as.numeric(flat.model["max",i])
    mx <- min(mx, 1-min.included-min.remaining)
    xi <- seq(from=as.numeric(flat.model["min",i]), to=mx, by=steps[i])
    xi <- sort(unique(c(xi, seq(from=min.step, to=num.min.steps*min.step, by=min.step))))
    mat <- merge(mat, as.data.frame(xi))
    colnames(mat) <- colnames(flat.model)[1:ncol(mat)]
    cat(paste0("Col: ", i, " min remaining = ", min.remaining, " mx = ", mx, " max = ", as.numeric(flat.model["max",i]), "\n"))
    flag <- unlist(apply(mat, 1, function(row) sum(row) <= (1-min.remaining)))
    mat <- mat[flag, ]
  }
  flag <- unlist(apply(mat, 1, function(row) sum(row) == 1))
  mat <- mat[flag, ]
  mat
}

sample.flat.unif.model.mixexp <- function(n, flat.model, min.prop = 0.01, constraint.col = "tumor.fraction") {
    library(mixexp)
##    save(file="fm.Rd", flat.model)
    nfac <- ncol(flat.model)
    ## Xvert only works in nfac < 12 dimensions
    df <- NULL
    nfac.max <- 10
    ncols <- ncol(flat.model)
    if(nfac < nfac.max) {
      df <- Xvert(n, uc = flat.model["max",], lc = flat.model["min", ], ndm = 0, plot = FALSE)
      df <- df[, !colnames(df) == "dimen"]
      colnames(df) <- colnames(flat.model)
    } else {
      uc <- as.numeric(flat.model["max",])
      lc <- as.numeric(flat.model["min",])
      uc1 <- uc[1:(nfac.max-1)]
      uc1[nfac.max-1] <- min(1,sum(uc[(nfac.max-1):ncols]))
      lc1 <- lc[1:(nfac.max-1)]
      lc1[nfac.max-1] <- sum(lc[(nfac.max-1):ncols])      
      df1 <- Xvert(nfac.max-1, uc = uc1, lc = lc1, ndm = 0, plot = FALSE)
      df1 <- df1[, !colnames(df1) == "dimen"]
      df1 <- round(100*df1)/100
      flag <- rowSums(df1) == 1
      df1 <- df1[flag,,drop=F]
      df1 <- df1[!(duplicated(df1)),,drop=F]
      library(plyr)
      library(dplyr)
      df2 <-
        ldply(1:5,
	      .fun = function(i) {
                       uc2 <- uc[c((nfac.max-1):ncols)]
                       lc2 <- lc[c((nfac.max-1):ncols)]
		       eps <- 10^-4
		       sm <- sum(df1[i, 1:(nfac.max-2)])
		       uc2 <- c(sm+eps/2, uc2)
		       lc2 <- c(sm-eps/2, lc2)
                       df2 <- Xvert(length(uc2), uc = uc2, lc = lc2, ndm = 0, plot = FALSE)
		       row <- as.numeric(df1[i, 1:(nfac.max-2)])
		       df2 <- df2[, 2:(length(uc2))]
		       cbind(as.data.frame(matrix(row, nrow=nrow(df2), ncol=length(row), byrow=TRUE)), df2)
                     })
      
##      uc2 <- uc[c(1, (nfac.max-1):ncols)]
##      lc2 <- lc[c(1, (nfac.max-1):ncols)]
        df <- df2
	colnames(df) <- colnames(flat.model)
    }
    df <- round(100*df)/100
    flag <- rowSums(df) == 1
    df <- df[flag,,drop=F]
    df <- df[!(duplicated(df)),,drop=F]
    df
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
