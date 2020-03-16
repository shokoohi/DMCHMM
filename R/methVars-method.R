.methVars <- function(object) {
    return(assays(object)$methVars)
}

.replace.methVars <- function(object, value) {
    assays(object, withDimnames = FALSE)$methVars <- value
    return(object)
}

#' @rdname methVars-method
#' @aliases methVars-method methVars
setMethod("methVars", signature(object = "BSDMCs"), .methVars)

#' @rdname methVars-method
#' @aliases methVars-method methVars<-
setReplaceMethod("methVars", signature(object = "BSDMCs", value = "matrix"),
    .replace.methVars)
