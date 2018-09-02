.cBSDMCs <- function(methReads, totalReads, methLevels, methStates, methVars,
    rowRanges, colData = DataFrame(row.names = colnames(methReads)),
    metadata = list(), ...) {
    new("BSDMCs", SummarizedExperiment(assays =
        SimpleList(methReads = methReads, totalReads = totalReads,
        methLevels = methLevels, methStates = methStates, methVars = methVars),
        rowRanges = rowRanges, colData = colData, metadata = list()))
}

#' @rdname cBSDMCs-method
#' @aliases cBSDMCs-method cBSDMCs
setMethod("cBSDMCs", signature(methReads = "matrix", totalReads = "matrix",
    methLevels = "matrix", methStates = "matrix", methVars = "matrix",
    rowRanges = "GRanges"),
    .cBSDMCs)
