.cBSData <- function(methReads, totalReads, rowRanges,
    colData = DataFrame(row.names = colnames(methReads)), metadata = list(),
    ...) {
    new("BSData", SummarizedExperiment(assays =
        SimpleList(totalReads = totalReads,
        methReads = methReads), rowRanges = rowRanges, colData = colData,
        metadata = list()))
}

#' @rdname cBSData-method
#' @aliases cBSData-method cBSData
setMethod("cBSData", signature(methReads = "matrix", totalReads = "matrix",
    rowRanges = "GRanges"), .cBSData)
