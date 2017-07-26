.BSDMCs <- function(methReads,
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

#' @rdname BSDMCs-method
setMethod("BSDMCs", signature(methReads = "matrix", totalReads = "matrix",
    rowRanges="GRanges"), .BSDMCs)
