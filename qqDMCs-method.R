.qqDMCs <- function(object, ...) {

    if (is.null(metadata(object)$DMCHMM)) stop("Input must be an object
    obtianed from findDMCs method.")

    pvector <- metadata(object)$DMCHMM$pvalues
    # limit to not missing, not nan, not null, not infinite, between 0 and 1
    pvector <- pvector[!is.na(pvector) & !is.nan(pvector) & !is.null(pvector) &
    is.finite(pvector) & pvector<1 & pvector>0]

    # Observed and expected
    o = -log10(sort(pvector,decreasing=FALSE))
    e = -log10( ppoints(length(pvector) ))


    # The new way to initialize the plot.
    ## See http://stackoverflow.com/q/23922130/654296
    ## First, define your default arguments
    def_args <- list(pch=20, xlim=c(0, max(e)), ylim=c(0, max(o)),
    xlab=expression(Expected~~-log[10](italic(p))),
    ylab=expression(Observed~~-log[10](italic(p)))
    )
    ## Next, get a list of ... arguments
    #dotargs <- as.list(match.call())[-1L]
    dotargs <- list(...)
    ## And call the plot function passing NA, your ... arguments, and the
    ## default arguments that were not defined in the ... arguments.
    tryCatch(do.call("plot", c(list(x=e, y=o), def_args[!names(def_args) %in%
    names(dotargs)], dotargs)), warn=stop)

    # Add diagonal
    abline(0,1,col="red")

}

#' @rdname qqDMCs-method
#' @aliases qqDMCs-method qqDMCs
setMethod("qqDMCs", signature(object = "BSDMCs"), .qqDMCs)

