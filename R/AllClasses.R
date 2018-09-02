#' @title params
#' @description parameters name and their descriptions
#' @author Farhad Shokoohi <shokoohi@icloud.com>
#' @name params
#' @param methReads The matrix \code{methReads} contains the number of
#' methylated reads spanning a CpG-site. The rows represent the CpG sites in
#' \code{rowRanges} and the columns represent the samples in \code{colData}.
#' @param totalReads The matrix \code{totalReads} contains the number of reads
#' spanning a CpG-site. The rows represent the CpG sites in \code{rowRanges}
#' and the columns represent the samples in \code{colData}.
#' @param methLevels The matrix \code{methLevels} contains the predicted
#' methylation level spanning a CpG-site using Hidden Markov model. The rows
#' represent the CpG sites in \code{rowRanges} and the columns represent the
#' samples in \code{colData}.
#' @param methVars The matrix \code{methVars} contains the variances of the
#' corresponding \code{methLevels} obtianed from MCMC.
#' @param methStates The matrix \code{methStates} contains the state of
#' methylation obtained from Hidden Markov model spanning a CpG-site. The rows
#' represent the CpG sites in \code{rowRanges} and the columns represent the
#' samples in \code{colData}. The value of state is stored in \code{metadata},
#' named \code{Beta}.
#' @param rowRanges A \code{\link{GRanges}} or \code{\link{GRangesList}}
#' object describing the ranges of interest. Names, if present, become the row
#' names of the \code{\link{SummarizedExperiment}} object. The length of the
#' \code{\link{GRanges}} or \code{\link{GRangesList}} must equal the number of
#' rows of the matrices in \code{assays}. If \code{rowRanges} is missing, a
#' \code{\link{SummarizedExperiment}} instance is returned.
#' @param colData Object of class \code{"DataFrame"} containing information on
#' variable values of the samples
#' @param metadata An optional \code{list} of arbitrary content describing the
#' overall experiment
#' @param object A \code{\link{BSData-class}} or \code{\link{BSDMCs-class}}
#' object
#' @param value An integer matrix
#' @param obj1 A \code{\link{BSData-class}} or \code{\link{BSDMCs-class}}
#' @param obj2 A \code{\link{BSData-class}} or \code{\link{BSDMCs-class}}
#' @param files A character list
#' @param file A character
#' @param name A character list
#' @param MaxK An integer value
#' @param MaxEmiter An integer value
#' @param epsEM A positive numeric value
#' @param useweight A logical value
#' @param mc.cores An integer greater than 0
#' @param nburn An integer value
#' @param nthin An integer value
#' @param nsamp An integer value
#' @param formula A formula
#' @param FDRthreshold A numeric value
#' @param Methylthreshold A numeric value
#' @param weightfunction A function to create weights using variance obtained
#' form the MCMC algorithm
#' @param ... other possible parameters
#' @param col A character vector indicating which colors to alternate.
#' @param chrlabs A character vector equal to the number of chromosomes
#' specifying the chromosome labels (e.g., \code{c(1:22, "X", "Y", "MT")}).
#' @param suggestiveline Where to draw a "suggestive" line. Default
#' -log10(1e-5). Set to FALSE to disable.
#' @param genomewideline Where to draw a "genome-wide sigificant" line. Default
#' -log10(5e-8). Set to FALSE to disable.
#' @param highlight A character vector of SNPs in your dataset to highlight.
#' These SNPs should all be in your dataset.
#' @param logp If TRUE, the -log10 of the p-value is plotted. It isn't very
#' useful to plot raw p-values, but plotting the raw value could be useful for
#' other genome-wide plots, for example, peak heights, bayes factors, test
#' statistics, other "scores," etc.
#' @param annotatePval If set, SNPs below this p-value will be annotated on the
#' plot.
#' @param annotateTop If TRUE, only annotates the top hit on each chromosome
#' that is below the annotatePval threshold.
#' @docType NULL
NULL

#' @title BSData object
#' @description The \code{BSData} object is an S4 class that represents BS-Seq
#' Data.
#' @author Farhad Shokoohi <shokoohi@icloud.com>
#' @name BSData-class
#' @import methods
#' @import SummarizedExperiment
#' @import S4Vectors
#' @slot methReads An integer matrix
#' @slot totalReads An integer matrix
#' @param methReads The matrix \code{methReads} contains the number of
#' methylated reads spanning a CpG-site. The rows represent the CpG sites in
#' \code{rowRanges} and the columns represent the samples in \code{colData}.
#' @param totalReads The matrix \code{totalReads} contains the number of reads
#' spanning a CpG-site. The rows represent the CpG sites in \code{rowRanges}
#' and the columns represent the samples in \code{colData}.
#' @seealso \code{\link{SummarizedExperiment}} objects.
#' @docType class
#' @keywords object
#' @return A \code{\link{BSData-class}} object
#' @examples
#' nr <- 500; nc <- 16
#' metht<-matrix(as.integer(runif(nr * nc, 0, nr)), nr)
#' methc<-matrix(rbinom(n=nr*nc,c(metht),prob = runif(nr*nc)),nr,nc)
#' r1 <- GRanges(rep("chr1", nr), IRanges(1:nr, width=1), strand="*")
#' names(r1) <- 1:nr
#' cd1<-DataFrame(Group=rep(c("G1","G2"),each=nc/2),row.names=LETTERS[1:nc])
#' OBJ1<-cBSData(rowRanges=r1,methReads=methc,totalReads=metht,colData=cd1)
#' OBJ1
#' @exportClass BSData
BSData <- setClass("BSData",
    representation(methReads = "matrix", totalReads = "matrix"),
    prototype (methReads = matrix(), totalReads = matrix()),
    contains = "RangedSummarizedExperiment")

setValidity("BSData", function(object){
    if(length(assays(object)) != 2)
    return("The assays slot in BSData object must be of length two.")
    if(!(all( is.element(names(assays(object)), c("methReads", "totalReads")))))
    return("The assays slot in BSData object must contain methReads and
    totalReads")
    if(!all( vapply(assays(object), class, character(1)) == "matrix" ))
    return("The methReads and totalReads slots of an BSData object must be
    matrices.")
    if(!all(vapply(assays(object), typeof, character(1)) == "integer" ))
    return("The methReads and totalReads matrices of an BSData object must
    contain integer data.")})

#' @title BSDMCs object
#' @description The \code{BSDMCs} object is an S4 class that represents
#' differentially methylated CpG sites (DMCs) in BS-Seq Data.
#' @author Farhad Shokoohi <shokoohi@icloud.com>
#' @name BSDMCs-class
#' @import methods
#' @import SummarizedExperiment
#' @slot methReads An integer matrix
#' @slot totalReads An integer matrix
#' @slot methLevels A numeric matrix
#' @slot methStates An integer matrix
#' @slot methVars A double matrix
#' @param methReads The matrix \code{methReads} contains the number of
#' methylated reads spanning a CpG-site. The rows represent the CpG sites in
#' \code{rowRanges} and the columns represent the samples in \code{colData}.
#' @param totalReads The matrix \code{totalReads} contains the number of reads
#' spanning a CpG-site. The rows represent the CpG sites in \code{rowRanges}
#' and the columns represent the samples in \code{colData}.
#' @param methLevels The matrix \code{methLevels} contains the predicted
#' methylation level spanning a CpG-site using Hidden Markov model. The rows
#' represent the CpG sites in \code{rowRanges} and the columns represent the
#' samples in \code{colData}.
#' @param methStates The matrix \code{methStates} contains the state of
#' methylation obtained from Hidden Markov model spanning a CpG-site. The rows
#' represent the CpG sites in \code{rowRanges} and the columns represent the
#' samples in \code{colData}. The value of state is stored in \code{metadata},
#' named \code{Beta}.
#' @param methVars The matrix \code{methVars} contains the variances of the
#' corresponding \code{methLevels} obtianed from MCMC.
#' @docType class
#' @keywords object
#' @return A \code{\link{BSDMCs-class}} object
#' @examples
#' nr <- 500; nc <- 16
#' metht <- matrix(as.integer(runif(nr * nc, 0, nr)), nr)
#' methc <- matrix(rbinom(n=nr*nc,c(metht),prob = runif(nr*nc)),nr,nc)
#' meths <- matrix(as.integer(runif(nr * nc, 0, 10)), nr)
#' methl <- methc/metht
#' methv <- matrix((runif(nr * nc, 0.1, 0.5)), nr)
#' r1 <- GRanges(rep('chr1', nr), IRanges(1:nr, width=1), strand='*')
#' names(r1) <- 1:nr
#' cd1 <- DataFrame(Group=rep(c('G1','G2'),each=nc/2),row.names=LETTERS[1:nc])
#' OBJ2 <- cBSDMCs(rowRanges=r1,methReads=methc,totalReads=metht,
#' methLevels=methl,methStates=meths,methVars=methv,colData=cd1)
#' OBJ2
#' @exportClass BSDMCs
BSDMCs <- setClass("BSDMCs",
    representation(methReads = "matrix", totalReads = "matrix",
    methLevels = "matrix", methStates = "matrix", methVars = "matrix"),
    prototype (methReads = matrix(), totalReads = matrix(),
    methLevels = matrix(), methStates = matrix(), methVars = matrix()),
    contains = "RangedSummarizedExperiment")

setValidity("BSDMCs", function(object){
    if(length(assays(object)) != 5)
    return("The assays slot in BSDMCs object must be of length five.")
    if(!(all( is.element(names(assays(object)), c("methReads", "totalReads",
    "methLevels", "methStates", "methVars")))))
    return("The assays slot in BSDMCs object must contain totalReads,
    methReads, methStates, methLevels and methVars.")
    if(!all( vapply(assays(object), class, character(1)) == "matrix" ))
    return("The methReads, totalReads, methLevels, methStates and methVars slots
    of a BSDMCs object must be matrices.")
    if(!all(vapply(assays(object), typeof, character(1)) ==
    c("integer", "integer", "double", "integer", "double")))
    return("The methReads, totalReads, methLevels, methStates and methVars slots
    of a BSDMCs object must be integer, integer, double, integer and double,
    repectively.")}
    )
