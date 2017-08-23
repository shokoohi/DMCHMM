.methReads <- function(object) {
    return(assays(object)$methReads)
}

.replace.methReads <- function(object, value) {
    assays(object)$methReads <- value
    return(object)
}

#' @rdname methReads-method
#' @aliases methReads-method methReads
setMethod("methReads", signature(object = "BSData"), .methReads)

#' @rdname methReads-method
#' @aliases methReads-method methReads<-
setReplaceMethod("methReads", signature(object = "BSData", value = "matrix"), 
    .replace.methReads)

#' @rdname methReads-method
#' @aliases methReads-method methReads
setMethod("methReads", signature(object = "BSDMCs"), .methReads)

#' @rdname methReads-method
#' @aliases methReads-method methReads<-
setReplaceMethod("methReads", signature(object = "BSDMCs", value = "matrix"), 
    .replace.methReads)

