.readBismark <- function(files, colData) {

    if (nrow(colData) != length(files)) {
        stop("Row number of colData must equal length of files.")
    }

    methData = list()

    for (i in seq_along(files)) {

        cat(paste("Processing sample ", rownames(colData)[i], " ... \n",
            sep = ""))

        bismark <- scan(files[i], skip = 0, sep = "\t", comment.char = "#",
            what = list("character", integer(), NULL, NULL, integer(),
                integer()))

        methData[[i]] = GRanges(seqnames = bismark[[1]],
            ranges = IRanges(start = bismark[[2]],
            width = 1), methylated = bismark[[5]], reads = bismark[[5]] +
            bismark[[6]])

        rm(bismark)
    }

    cat("Building BSData object.\n")

    fData <- methData[[1]]

    if (length(methData) > 1) {
        for (i in seq(along = methData)[-1]) {
            fData <- unique(c(fData, methData[[i]]))
        }
    }
    mcols(fData) <- NULL
    names(fData) <- as.character(1:length(fData))

    tReads <- matrix(integer(length = length(fData) * length(methData)),
        nrow = length(fData))
    mReads <- matrix(integer(length = length(fData) * length(methData)),
        nrow = length(fData))

    for (i in seq_along(methData)) {
        mtch <- findOverlaps(fData, methData[[i]])
        ind1 <- queryHits(mtch)
        ind2 <- subjectHits(mtch)
        tReads[ind1, i] <- mcols(methData[[i]])$reads[ind2]
        mReads[ind1, i] <- mcols(methData[[i]])$methylated[ind2]
    }

    colnames(tReads) <- rownames(colData)
    colnames(mReads) <- rownames(colData)
    rownames(tReads) <- names(fData)
    rownames(mReads) <- names(fData)

    BSDAT = cBSData(colData = colData, rowRanges = fData, methReads = mReads,
        totalReads = tReads)

    return(BSDAT)
}

#' @rdname readBismark-method
#' @aliases readBismark-method readBismark
setMethod("readBismark", signature(files = "character", colData = "DataFrame"),
    .readBismark)

#' @rdname readBismark-method
#' @aliases readBismark-method readBismark
setMethod("readBismark",
          signature = c(files = "character", colData = "data.frame"),
          function(files, colData) {
              colData = as(colData, "DataFrame")
              .readBismark(files, colData)
    })

#' @rdname readBismark-method
#' @aliases readBismark-method readBismark
setMethod("readBismark",
          signature = c(files = "character", colData = "character"),
          function(files, colData) {
              colData = DataFrame(row.names = colData)
              .readBismark(files, colData)
    })
