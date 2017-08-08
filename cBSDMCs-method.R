.cBSDMCs <- function(methReads,
    totalReads,
    methLevels,
    methStates,
    rowRanges,
    colData = DataFrame(row.names=colnames(methReads)),
    metadata = list(),
    ...)
{
    new("BSDMCs", SummarizedExperiment(
    assays = SimpleList(methReads = methReads,
    totalReads = totalReads,
    methLevels = methLevels,
    methStates = methStates),
    rowRanges = rowRanges,
    colData = colData,
    metadata = list()))
}

#' @rdname cBSDMCs-method
#' @aliases cBSDMCs-method cBSDMCs
setMethod("cBSDMCs", signature(methReads = "matrix", totalReads = "matrix",
    methLevels = "matrix", methStates = "matrix", rowRanges = "GRanges"),
    .cBSDMCs)
