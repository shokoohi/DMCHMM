.readBismark <- function(files, colData, mc.cores) {

    if (nrow(colData) != length(files)) {
        stop("Row number of colData must equal length of files.")
    }

    methData = list()

    optbp <- MulticoreParam(workers = mc.cores, progressbar = TRUE)

    .runfile <- function(i){
        cat(paste("\nProcessing sample ", rownames(colData)[i], " ... \n",
            sep = ""))

        bismark <- scan(files[i], skip = 0, sep = "\t", comment.char = "#",
            what = list("character", integer(), NULL, NULL, integer(),
                integer()))

        Datafile = GRanges(seqnames = bismark[[1]],
            ranges = IRanges(start = bismark[[2]],
                width = 1), methylated = bismark[[5]], reads = bismark[[5]] +
                    bismark[[6]])

        rm(bismark)
        return(Datafile)
    }

    methData <- bplapply(seq_along(files), .runfile, BPPARAM = optbp)


    cat("Building BSData object.\n")

    fData <- methData[[1]]

    if (length(methData) > 1) {
        for (i in seq(along = methData)[-1]) {
            fData <- unique(c(fData, methData[[i]]))
        }
    }
    mcols(fData) <- NULL
    names(fData) <- as.character(seq_len(length(fData)))

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
setMethod("readBismark", signature(files = "character", colData = "DataFrame",
    mc.cores = "numeric"), .readBismark)

#' @rdname readBismark-method
#' @aliases readBismark-method readBismark
setMethod("readBismark", signature = c(files = "character",
        colData = "data.frame", mc.cores = "numeric"),
            function(files, colData, mc.cores) {
                colData = as(colData, "DataFrame")
                .readBismark(files, colData, mc.cores)
    })

#' @rdname readBismark-method
#' @aliases readBismark-method readBismark
setMethod("readBismark", signature = c(files = "character",
        colData = "character", mc.cores = "numeric"),
            function(files, colData, mc.cores) {
                colData = DataFrame(row.names = colData)
                .readBismark(files, colData, mc.cores)
    })
