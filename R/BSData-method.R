.BSData <- function(methReads,
    totalReads,
    rowRanges,
    colData = DataFrame(row.names=colnames(methReads)),
    metadata = list(),
    ...)
{ new("BSData", SummarizedExperiment( assays = SimpleList(
    totalReads = totalReads,
    methReads = methReads
    ),
    rowRanges = rowRanges,
    colData = colData,
    metadata = list())
    )
}

#' @rdname BSData-method
#' @aliases BSData-method
setMethod("BSData", signature(methReads = "matrix", totalReads = "matrix",
    rowRanges="GRanges"), .BSData)
