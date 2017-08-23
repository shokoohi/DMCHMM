.methStates <- function(object) {
    return(assays(object)$methStates)
}

.replace.methStates <- function(object, value) {
    assays(object)$methStates <- value
    return(object)
}

#' @rdname methStates-method
#' @aliases methStates-method methStates
setMethod("methStates", signature(object = "BSDMCs"), .methStates)

#' @rdname methStates-method
#' @aliases methStates-method methStates<-
setReplaceMethod("methStates", signature(object = "BSDMCs", value = "matrix"), 
    .replace.methStates)
