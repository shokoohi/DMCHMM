.combine.BSData <- function(obj1,obj2){
    if(any(is.element(colnames(obj1), colnames(obj2)))){
        stop("The BSData objects to combine should not have samples in common!")
    }
    colData.new <- rbind(colData(obj1), colData(obj2))
    rowRanges.new <- sort(unique(c(rowRanges(obj1), rowRanges(obj2))))
    ind.match.obj1<-findOverlaps(rowRanges.new, rowRanges(obj1), select="first")
    ind.match.obj2<-findOverlaps(rowRanges.new, rowRanges(obj2), select="first")
    nr <- length(rowRanges.new)
    nc <- nrow(colData.new)
    methReads.new <- matrix(integer(length = nr*nc), ncol=nc, nrow=nr,
    dimnames=list(names(rowRanges.new), rownames(colData.new)))
    methReads.new[,] <- cbind(methReads(obj1)[ind.match.obj1,],
    methReads(obj2)[ind.match.obj2,])
    methReads.new[is.na(methReads.new)] <- 0L
    totalReads.new <- matrix(integer(length = nr*nc), ncol=nc, nrow=nr,
    dimnames=list(names(rowRanges.new), rownames(colData.new)))
    totalReads.new[,] <- cbind(totalReads(obj1)[ind.match.obj1,],
    totalReads(obj2)[ind.match.obj2,])
    totalReads.new[is.na(totalReads.new)] <- 0L
    z <- cBSData(colData = colData.new, rowRanges = rowRanges.new,
    methReads = methReads.new, totalReads = totalReads.new
    )
    return(z)
}

.combine.BSDMCs <- function(obj1,obj2){
    if(any(is.element(colnames(obj1), colnames(obj2)))){
        stop("The BSDMCs objects to combine should not have samples in common!")
    }
    colData.new <- rbind(colData(obj1), colData(obj2))
    rowRanges.new <- sort(unique(c(rowRanges(obj1), rowRanges(obj2))))
    ind.match.obj1<-findOverlaps(rowRanges.new, rowRanges(obj1), select="first")
    ind.match.obj2<-findOverlaps(rowRanges.new, rowRanges(obj2), select="first")
    nr <- length(rowRanges.new)
    nc <- nrow(colData.new)
    methReads.new <- matrix(integer(length = nr*nc), ncol=nc, nrow=nr,
    dimnames=list(names(rowRanges.new), rownames(colData.new)))
    methReads.new[,] <- cbind(methReads(obj1)[ind.match.obj1,],
    methReads(obj2)[ind.match.obj2,])
    methReads.new[is.na(methReads.new)] <- 0L
    totalReads.new <- matrix(integer(length = nr*nc), ncol=nc, nrow=nr,
    dimnames=list(names(rowRanges.new), rownames(colData.new)))
    totalReads.new[,] <- cbind(totalReads(obj1)[ind.match.obj1,],
    totalReads(obj2)[ind.match.obj2,])
    totalReads.new[is.na(totalReads.new)] <- 0L
    methStates.new <- matrix(integer(length = nr*nc), ncol=nc, nrow=nr,
    dimnames=list(names(rowRanges.new), rownames(colData.new)))
    methStates.new[,] <- cbind(methStates(obj1)[ind.match.obj1,],
    methStates(obj2)[ind.match.obj2,])
    methStates.new[is.na(methStates.new)] <- 0L
    methLevels.new <- matrix(integer(length = nr*nc), ncol=nc, nrow=nr,
    dimnames=list(names(rowRanges.new), rownames(colData.new)))
    methLevels.new[,] <- cbind(methLevels(obj1)[ind.match.obj1,],
    methLevels(obj2)[ind.match.obj2,])
    methLevels.new[is.na(methLevels.new)] <- 0L
    z <- cBSDMCs(
    methReads = methReads.new, totalReads = totalReads.new,
    methLevels = methLevels.new, methStates = methStates.new,
    rowRanges = rowRanges.new, colData = colData.new
    )
    return(z)
}

#' @rdname combine-method
#' @aliases combine-method combine
setMethod("combine", signature(obj1 = "BSData", obj2 = "BSData"),
    .combine.BSData)

#' @rdname combine-method
#' @aliases combine-method combine
setMethod("combine", signature(obj1 = "BSDMCs", obj2 = "BSDMCs"),
    .combine.BSDMCs)
