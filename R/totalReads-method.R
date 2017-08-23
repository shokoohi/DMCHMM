.totalReads <- function(object) {
    return(assays(object)$totalReads)
}

.replace.totalReads <- function(object, value) {
    assays(object)$totalReads <- value
    return(object)
}

#' @rdname totalReads-method
#' @aliases totalReads-method totalReads
setMethod("totalReads", signature(object = "BSData"), .totalReads)

#' @rdname totalReads-method
#' @aliases totalReads-method totalReads<-
setReplaceMethod("totalReads", signature(object = "BSData", value = "matrix"), 
    .replace.totalReads)

#' @rdname totalReads-method
#' @aliases totalReads-method totalReads
setMethod("totalReads", signature(object = "BSDMCs"), .totalReads)

#' @rdname totalReads-method
#' @aliases totalReads-method totalReads<-
setReplaceMethod("totalReads", signature(object = "BSDMCs", value = "matrix"), 
    .replace.totalReads)
