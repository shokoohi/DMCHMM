.methLevels <- function(object) {
    return(assays(object)$methLevels)
}

.replace.methLevels <- function(object, value) {
    assays(object)$methLevels <- value
    return(object)
}

#' @rdname methLevels-method
#' @aliases methLevels-method methLevels
setMethod("methLevels", signature(object = "BSDMCs"), .methLevels)

#' @rdname methLevels-method
#' @aliases methLevels-method methLevels<-
setReplaceMethod("methLevels", signature(object = "BSDMCs", value = "matrix"), 
    .replace.methLevels)
