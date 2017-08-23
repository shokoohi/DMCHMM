.writeBED_BSData <- function(object, name, file) {

    if (ncol(object) != length(name) | ncol(object) != length(file)) {
        stop("Character vectors name and file must have the same length as the
    object has samples.")
    }

    colFunc <- colorRamp(colors = c("green", "black", "red"))
    strand(object) <- "*"
    object = sort(object)

    for (i in 1:ncol(object)) {
        object.i <- object[, i]
        ind.cov <- totalReads(object.i) != 0
        if (sum(ind.cov) > 1) {
            object.i <- object.i[ind.cov, ]
            bed <- rowRanges(object.i)
            mcols(bed)$score <- methReads(object.i)/totalReads(object.i)
            mcols(bed)$name <- totalReads(object.i)
            m <- colFunc(mcols(bed)$score)/255
            mcols(bed)$itemRgb <- rgb(m[, 1], m[, 2], m[, 3])

            bed <- as(bed, "UCSCData")
            bed@trackLine@name <- paste("\"", name[i], "\"", sep = "")
            export.ucsc(bed, con = file[i], subformat = "bed")

        } else {
            # no CpGs are covered
            warning(paste(name[i], ": No CpG is covered - no bed file has
    been written",
                sep = ""))
        }
    }
}

.writeBED_BSDMCs <- function(object, name, file) {

    if (ncol(object) != length(name) | ncol(object) != length(file)) {
        stop("Character vectors name and file must have the same length as
    the object has samples.")
    }

    colFunc <- colorRamp(colors = c("green", "black", "red"))
    object = sort(object)

    for (i in 1:ncol(object)) {
        object.i <- object[, i]
        ind.cov <- !is.na(methLevels(object.i))
        if (sum(ind.cov) > 1) {
            object.i <- object.i[ind.cov, ]
            bed <- rowRanges(object.i)
            mcols(bed)$score <- methLevels(object.i)
            mcols(bed)$name <- "\"\""
            m <- colFunc(mcols(bed)$score)/255
            mcols(bed)$itemRgb <- rgb(m[, 1], m[, 2], m[, 3])

            bed <- as(bed, "UCSCData")
            bed@trackLine@name <- paste("\"", name[i], "\"", sep = "")
            export.ucsc(bed, con = file[i], subformat = "bed")

        } else {
            # no CpGs are covered
            warning(paste(name[i], ": No CpG is covered - no bed file has been
    written",
                sep = ""))
        }
    }
}

#' @rdname writeBED-method
#' @aliases writeBED-method writeBED
setMethod("writeBED", signature = c(object = "BSData", name = "character",
    file = "character"), .writeBED_BSData)

#' @rdname writeBED-method
#' @aliases writeBED-method writeBED
setMethod("writeBED", signature = c(object = "BSData", name = "character",
    file = "missing"), function(object, name) {
    .writeBED_BSData(object, name = name, file = paste(colnames(object),
        ".bed", sep = ""))
})

#' @rdname writeBED-method
#' @aliases writeBED-method writeBED
setMethod("writeBED", signature = c(object = "BSData", name = "missing",
    file = "character"), function(object, file) {
    .writeBED_BSData(object, name = colnames(object), file = file)
})

#' @rdname writeBED-method
#' @aliases writeBED-method writeBED
setMethod("writeBED", signature = c(object = "BSData", name = "missing",
    file = "missing"), function(object) {
    .writeBED_BSData(object, name = colnames(object),
                     file = paste(colnames(object), ".bed", sep = ""))
})

#' @rdname writeBED-method
#' @aliases writeBED-method writeBED
setMethod("writeBED", signature = c(object = "BSDMCs", name = "character",
    file = "character"), .writeBED_BSDMCs)

#' @rdname writeBED-method
#' @aliases writeBED-method writeBED
setMethod("writeBED", signature = c(object = "BSDMCs", name = "character",
    file = "missing"), function(object, name) {
    .writeBED_BSDMCs(object, name = name, file = paste(colnames(object),
        ".bed", sep = ""))
})

#' @rdname writeBED-method
#' @aliases writeBED-method writeBED
setMethod("writeBED", signature = c(object = "BSDMCs", name = "missing",
    file = "character"), function(object, file) {
    .writeBED_BSDMCs(object, name = colnames(object), file = file)
})

#' @rdname writeBED-method
#' @aliases writeBED-method writeBED
setMethod("writeBED", signature = c(object = "BSDMCs", name = "missing",
    file = "missing"), function(object) {
    .writeBED_BSDMCs(object, name = colnames(object),
                     file = paste(colnames(object), ".bed", sep = ""))
})

