#' @title cBSData method
#' @description Creates a \code{\link{BSData-class}} object
#' @author Farhad Shokoohi <shokoohi@icloud.com>
#' @name cBSData-method
#' @import S4Vectors
#' @inheritParams params
#' @details The rows of a \code{BSData} object represent ranges (in genomic
#' coordinates) of interest. The ranges of interest are described by a
#' \code{GRanges} or a \code{GRangesList} object,
#' accessible using the \code{rowRanges} function.
#' The \code{GRanges} and \code{GRangesList} classes
#' contains sequence (e.g., chromosome) name, genomic coordinates, and strand
#' information. Each range can be annotated with additional data; this data
#' might be used to describe the range or to
#' summarize results (e.g., statistics of differential abundance) relevant to
#' the range. Rows may or may not have row names; they often will not.
#' @return A \code{\link{BSData-class}} object
#' @examples
#' nr <- 150; nc <- 8
#' metht <- matrix(as.integer(runif(nr * nc, 0, 100)), nr)
#' methc <- matrix(rbinom(n=nr*nc,c(metht),prob = runif(nr*nc)),nr,nc)
#' r1 <- GRanges(rep('chr1', nr), IRanges(1:nr, width=1), strand='*')
#' names(r1) <- 1:nr
#' cd1 <- DataFrame(Group=rep(c('G1','G2'),each=nc/2),row.names=LETTERS[1:nc])
#' OBJ1 <- cBSData(rowRanges=r1,methReads=methc,totalReads=metht,colData=cd1)
#' OBJ1
#' @exportMethod cBSData
setGeneric("cBSData", function(methReads, totalReads, rowRanges, colData =
                            DataFrame(row.names = colnames(methReads)),
                            metadata = list(), ...) standardGeneric("cBSData"),
    signature = c("methReads", "totalReads", "rowRanges"))

#' @title methReads method
#' @description Returns \code{methReads} stored in \code{\link{BSData-class}}
#' @author Farhad Shokoohi <shokoohi@icloud.com>
#' @import S4Vectors
#' @name methReads-method
#' @inheritParams params
#' @return A matrix
#' @examples
#' nr <- 150; nc <- 8
#' metht <- matrix(as.integer(runif(nr * nc, 0, 100)), nr)
#' methc <- matrix(rbinom(n=nr*nc,c(metht),prob = runif(nr*nc)),nr,nc)
#' r1 <- GRanges(rep('chr1', nr), IRanges(1:nr, width=1), strand='*')
#' names(r1) <- 1:nr
#' cd1 <- DataFrame(Group=rep(c('G1','G2'),each=nc/2),row.names=LETTERS[1:nc])
#' OBJ1 <- cBSData(rowRanges=r1,methReads=methc,totalReads=metht,colData=cd1)
#' methReads(OBJ1)
#' @exportMethod methReads
setGeneric("methReads", function(object) standardGeneric("methReads"))

#' @title methReads method
#' @description Assigns \code{methReads} to \code{\link{BSData-class}}
#' @name methReads-method
#' @import S4Vectors
#' @inheritParams params
#' @return A \code{\link{BSData-class}} object
#' @examples
#' methReads(OBJ1) <- methc
#' @exportMethod methReads<-
setGeneric("methReads<-", function(object,value) standardGeneric("methReads<-"))

#' @title totalReads method
#' @description Returns \code{totalReads} stored in \code{\link{BSData-class}}
#' @author Farhad Shokoohi <shokoohi@icloud.com>
#' @name totalReads-method
#' @import S4Vectors
#' @inheritParams params
#' @return A matrix
#' @examples
#' nr <- 150; nc <- 8
#' metht <- matrix(as.integer(runif(nr * nc, 0, 100)), nr)
#' methc <- matrix(rbinom(n=nr*nc,c(metht),prob = runif(nr*nc)),nr,nc)
#' r1 <- GRanges(rep('chr1', nr), IRanges(1:nr, width=1), strand='*')
#' names(r1) <- 1:nr
#' cd1 <- DataFrame(Group=rep(c('G1','G2'),each=nc/2),row.names=LETTERS[1:nc])
#' OBJ1 <- cBSData(rowRanges=r1,methReads=methc,totalReads=metht,colData=cd1)
#' totalReads(OBJ1)
#' @exportMethod totalReads
setGeneric("totalReads", function(object) standardGeneric("totalReads"))

#' @title totalReads method
#' @description Assigns \code{totalReads} to \code{\link{BSData-class}}
#' @name totalReads-method
#' @import S4Vectors
#' @inheritParams params
#' @return A \code{\link{BSData-class}} object
#' @examples
#' totalReads(OBJ1) <- metht
#' @exportMethod totalReads<-
setGeneric("totalReads<-", function(object,value)
    standardGeneric("totalReads<-"))

#' @title cBSDMCs method
#' @description Creates a \code{\link{BSDMCs-class}} object
#' @author Farhad Shokoohi <shokoohi@icloud.com>
#' @name cBSDMCs-method
#' @import S4Vectors
#' @inheritParams params
#' @details The rows of a \code{BSDMCs} object represent ranges (in genomic
#' coordinates) of interest. The ranges of interest are described by a
#' \code{GRanges} or a \code{GRangesList} object,
#' accessible using the \code{rowRanges} function.
#' The \code{GRanges} and \code{GRangesList} classes
#' contains sequence (e.g., chromosome) name, genomic coordinates, and strand
#' information. Each range can be annotated with additional data; this data
#' might be used to describe the range or to
#' summarize results (e.g., statistics of differential abundance) relevant to
#' the range. Rows may or may not have row names; they often will not.
#' @return A \code{\link{BSDMCs-class}}
#' @examples
#' set.seed(1980)
#' nr <- 150; nc <- 8
#' metht <- matrix(as.integer(runif(nr * nc, 0, 100)), nr)
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
#' @exportMethod cBSDMCs
setGeneric("cBSDMCs", function(methReads, totalReads, methLevels, methStates,
        methVars, rowRanges,
        colData = DataFrame(row.names = colnames(methReads)),
        metadata = list(), ...) standardGeneric("cBSDMCs"), signature =
        c("methReads", "totalReads", "methLevels", "methStates", "methVars",
        "rowRanges"))

#' @title methReads method
#' @description Returns \code{methReads} stored in \code{\link{BSDMCs-class}}
#' @name methReads-method
#' @import S4Vectors
#' @inheritParams params
#' @return A matrix
#' @exportMethod methReads
setGeneric("methReads", function(object) standardGeneric("methReads"))

#' @title methReads method
#' @description Assigns \code{methReads} to \code{\link{BSDMCs-class}}
#' @name methReads-method
#' @import S4Vectors
#' @inheritParams params
#' @return A \code{\link{BSDMCs-class}} object
#' @exportMethod methReads<-
setGeneric("methReads<-", function(object, value)
    standardGeneric("methReads<-"))

#' @title totalReads method
#' @description Returns \code{totalReads} stored in \code{\link{BSDMCs-class}}
#' @name totalReads-method
#' @import S4Vectors
#' @inheritParams params
#' @return A matrix
#' @exportMethod totalReads
setGeneric("totalReads", function(object) standardGeneric("totalReads"))

#' @title totalReads method
#' @description Assigns \code{totalReads} to \code{\link{BSDMCs-class}}
#' @name totalReads-method
#' @import S4Vectors
#' @inheritParams params
#' @return A \code{\link{BSDMCs-class}} object
#' @exportMethod totalReads<-
setGeneric("totalReads<-", function(object, value)
    standardGeneric("totalReads<-"))

#' @title methLevels method
#' @description Returns \code{methLevels} stored in \code{\link{BSDMCs-class}}
#' @author Farhad Shokoohi <shokoohi@icloud.com>
#' @name methLevels-method
#' @import S4Vectors
#' @inheritParams params
#' @return A matrix
#' @examples
#' set.seed(1980)
#' nr <- 150; nc <- 8
#' metht <- matrix(as.integer(runif(nr * nc, 0, 100)), nr)
#' methc <- matrix(rbinom(n=nr*nc,c(metht),prob = runif(nr*nc)),nr,nc)
#' meths <- matrix(as.integer(runif(nr * nc, 0, 10)), nr)
#' methl <- methc/metht
#' methv <- matrix((runif(nr * nc, 0.1, 0.5)), nr)
#' r1 <- GRanges(rep('chr1', nr), IRanges(1:nr, width=1), strand='*')
#' names(r1) <- 1:nr
#' cd1 <- DataFrame(Group=rep(c('G1','G2'),each=nc/2),row.names=LETTERS[1:nc])
#' OBJ2 <- cBSDMCs(rowRanges=r1,methReads=methc,totalReads=metht,
#' methLevels=methl,methStates=meths,methVars=methv,colData=cd1)
#' methLevels(OBJ2)
#' @exportMethod methLevels
setGeneric("methLevels", function(object) standardGeneric("methLevels"))

#' @title methLevels method
#' @description Assigns \code{methLevels} to \code{\link{BSDMCs-class}}
#' @name methLevels-method
#' @import S4Vectors
#' @inheritParams params
#' @return A \code{\link{BSDMCs-class}} object
#' @examples
#' methLevels(OBJ2) <- methl
#' @exportMethod methLevels<-
setGeneric("methLevels<-", function(object, value)
    standardGeneric("methLevels<-"))

#' @title methStates method
#' @description Returns \code{methStates} stored in \code{\link{BSDMCs-class}}
#' @author Farhad Shokoohi <shokoohi@icloud.com>
#' @name methStates-method
#' @import S4Vectors
#' @inheritParams params
#' @return A matrix
#' @examples
#' set.seed(1980)
#' nr <- 150; nc <- 8
#' metht <- matrix(as.integer(runif(nr * nc, 0, 100)), nr)
#' methc <- matrix(rbinom(n=nr*nc,c(metht),prob = runif(nr*nc)),nr,nc)
#' meths <- matrix(as.integer(runif(nr * nc, 0, 10)), nr)
#' methl <- methc/metht
#' methv <- matrix((runif(nr * nc, 0.1, 0.5)), nr)
#' r1 <- GRanges(rep('chr1', nr), IRanges(1:nr, width=1), strand='*')
#' names(r1) <- 1:nr
#' cd1 <- DataFrame(Group=rep(c('G1','G2'),each=nc/2),row.names=LETTERS[1:nc])
#' OBJ2 <- cBSDMCs(rowRanges=r1,methReads=methc,totalReads=metht,
#' methLevels=methl,methStates=meths,methVars=methv,colData=cd1)
#' methStates(OBJ2)
#' @exportMethod methStates
setGeneric("methStates", function(object) standardGeneric("methStates"))

#' @title methStates method
#' @description Assigns \code{methStates} to \code{\link{BSDMCs-class}}
#' @name methStates-method
#' @import S4Vectors
#' @inheritParams params
#' @return A \code{\link{BSDMCs-class}} object
#' @examples
#' methStates(OBJ2)<- meths
#' @exportMethod methStates<-
setGeneric("methStates<-", function(object, value)
    standardGeneric("methStates<-"))

#' @title methVars method
#' @description Returns \code{methVars} stored in \code{\link{BSDMCs-class}}
#' @author Farhad Shokoohi <shokoohi@icloud.com>
#' @name methVars-method
#' @import S4Vectors
#' @inheritParams params
#' @return A matrix
#' @examples
#' set.seed(1980)
#' nr <- 150; nc <- 8
#' metht <- matrix(as.integer(runif(nr * nc, 0, 100)), nr)
#' methc <- matrix(rbinom(n=nr*nc,c(metht),prob = runif(nr*nc)),nr,nc)
#' meths <- matrix(as.integer(runif(nr * nc, 0, 10)), nr)
#' methl <- methc/metht
#' methv <- matrix((runif(nr * nc, 0.1, 0.5)), nr)
#' r1 <- GRanges(rep('chr1', nr), IRanges(1:nr, width=1), strand='*')
#' names(r1) <- 1:nr
#' cd1 <- DataFrame(Group=rep(c('G1','G2'),each=nc/2),row.names=LETTERS[1:nc])
#' OBJ2 <- cBSDMCs(rowRanges=r1,methReads=methc,totalReads=metht,
#' methLevels=methl,methStates=meths,methVars=methv,colData=cd1)
#' methVars(OBJ2)
#' @exportMethod methVars
setGeneric("methVars", function(object) standardGeneric("methVars"))

#' @title methVars method
#' @description Assigns \code{methVars} to \code{\link{BSDMCs-class}}
#' @name methVars-method
#' @import S4Vectors
#' @inheritParams params
#' @return A \code{\link{BSDMCs-class}} object
#' @examples
#' methVars(OBJ2)<- meths
#' @exportMethod methVars<-
setGeneric("methVars<-", function(object, value)
    standardGeneric("methVars<-"))

#' @title combine method
#' @description combine two \code{\link{BSData-class}} or
#' two \code{\link{BSDMCs-class}}
#' @author Farhad Shokoohi <shokoohi@icloud.com>
#' @name combine-method
#' @import S4Vectors
#' @inheritParams params
#' @return A \code{\link{BSData-class}} or \code{\link{BSDMCs-class}}
#' @examples
#' set.seed(1980)
#' nr <- 150; nc <- 8
#' metht <- matrix(as.integer(runif(nr * nc*2, 0, nr)), nr)
#' methc <- matrix(rbinom(n=nr*nc,c(metht),prob = runif(nr*nc*2)),nr,nc*2)
#' r1 <- GRanges(rep('chr1', nr), IRanges(1:nr, width=1), strand='*')
#' names(r1) <- 1:nr
#' cd1 <- DataFrame(Group=rep('G1',each=nc),row.names=LETTERS[1:nc])
#' OBJ1 <- cBSData(rowRanges=r1,methReads=methc[,1:nc],totalReads=metht[,1:nc],
#' colData=cd1)
#' cd2 <- DataFrame(Group=rep('G2',each=nc),row.names=LETTERS[nc+1:nc])
#' OBJ2 <- cBSData(rowRanges=r1,methReads=methc[,nc+1:nc],totalReads=
#' metht[,nc+1:nc],colData=cd2)
#' OBJ3 <- combine(OBJ1, OBJ2)
#' OBJ3
#' @exportMethod combine
setGeneric("combine", function(obj1, obj2) standardGeneric("combine"))

#' @title readBismark method
#' @description reads BS-Seq data
#' @author Farhad Shokoohi <shokoohi@icloud.com>
#' @name readBismark-method
#' @import GenomicRanges
#' @import IRanges
#' @import S4Vectors
#' @inheritParams params
#' @return A \code{\link{BSData-class}} object
#' @examples
#' fn <- list.files(system.file('extdata',package = 'DMCHMM'))
#' fn.f <- list.files(system.file('extdata',package='DMCHMM'), full.names=TRUE)
#' OBJ <- readBismark(fn.f, fn, mc.cores=2)
#' cdOBJ <- DataFrame(Cell = factor(c('BC', 'TC','Mono'),
#' labels = c('BC', 'TC', 'Mono')), row.names = c('BCU1568','BCU173','BCU551'))
#' colData(OBJ) <- cdOBJ
#' OBJ
#' @exportMethod readBismark
setGeneric("readBismark", function(files, colData, mc.cores)
    standardGeneric("readBismark"))

#' @title writeBED method
#' @description write BS-Seq data to BED files
#' @author Farhad Shokoohi <shokoohi@icloud.com>
#' @name writeBED-method
#' @inheritParams params
#' @import S4Vectors
#' @return BED files
#' @import grDevices
#' @import rtracklayer
#' @exportMethod writeBED
setGeneric("writeBED", function(object, name, file) standardGeneric("writeBED"))

#' @title methHMEM method
#' @description Estimates the HMM methylation paths and the HMM order for
#' each sample using the EM algorithm
#' @author Farhad Shokoohi <shokoohi@icloud.com>
#' @name methHMEM-method
#' @inheritParams params
#' @return \code{\link{BSDMCs-class}} object
#' @import BiocParallel
#' @import GenomicRanges
#' @import S4Vectors
#' @importFrom stats dbinom
#' @examples
#' set.seed(1980)
#' nr <- 150; nc <- 8
#' metht <- matrix(as.integer(runif(nr * nc, 0, 100)), nr)
#' methc <- matrix(rbinom(n=nr*nc,c(metht),prob = runif(nr*nc)),nr,nc)
#' r1 <- GRanges(rep('chr1', nr), IRanges(1:nr, width=1), strand='*')
#' names(r1) <- 1:nr
#' cd1 <- DataFrame(Group=rep(c('G1','G2'),each=nc/2),row.names=LETTERS[1:nc])
#' OBJ1 <- cBSData(rowRanges=r1,methReads=methc,totalReads=metht,colData=cd1)
#' OBJ2 <- methHMEM(OBJ1, MaxK=2, mc.cores=2)
#' OBJ2
#' @exportMethod methHMEM
setGeneric("methHMEM", function(object, MaxK, MaxEmiter, epsEM, useweight,
    mc.cores) standardGeneric("methHMEM"))

#' @title methHMMCMC method
#' @description Estimates the HMM methylation paths and the HMM order for
#' each sample using the MCMC algorithm
#' @author Farhad Shokoohi <shokoohi@icloud.com>
#' @name methHMMCMC-method
#' @inheritParams params
#' @return \code{\link{BSDMCs-class}} object
#' @import BiocParallel
#' @import GenomicRanges
#' @import S4Vectors
#' @importFrom stats dbinom
#' @importFrom stats lm
#' @importFrom stats pbeta
#' @importFrom stats qbeta
#' @importFrom stats runif
#' @importFrom utils combn
#' @examples
#' set.seed(1980)
#' nr <- 150; nc <- 8
#' metht <- matrix(as.integer(runif(nr * nc, 0, 100)), nr)
#' methc <- matrix(rbinom(n=nr*nc,c(metht),prob = runif(nr*nc)),nr,nc)
#' r1 <- GRanges(rep('chr1', nr), IRanges(1:nr, width=1), strand='*')
#' names(r1) <- 1:nr
#' cd1 <- DataFrame(Group=rep(c('G1','G2'),each=nc/2),row.names=LETTERS[1:nc])
#' OBJ1 <- cBSData(rowRanges=r1,methReads=methc,totalReads=metht,colData=cd1)
#' OBJ2 <- methHMEM(OBJ1, MaxK=2, mc.cores=2)
#' OBJ3 <- methHMMCMC(OBJ2, mc.cores=2)
#' OBJ3
#' @exportMethod methHMMCMC
setGeneric("methHMMCMC", function(object, useweight, nburn, nthin, nsamp,
    mc.cores) standardGeneric("methHMMCMC"))

#' @title findDMCs method
#' @description finds the DMCs after smoothing using HMM
#' @author Farhad Shokoohi <shokoohi@icloud.com>
#' @name findDMCs-method
#' @inheritParams params
#' @return \code{\link{BSDMCs-class}} object
#' @import BiocParallel
#' @import rtracklayer
#' @import GenomicRanges
#' @import fdrtool
#' @import multcomp
#' @import S4Vectors
#' @importFrom stats dbinom
#' @importFrom stats as.formula
#' @importFrom stats pf
#' @examples
#' set.seed(1980)
#' nr <- 150; nc <- 8
#' metht <- matrix(as.integer(runif(nr * nc, 0, 100)), nr)
#' methc <- matrix(rbinom(n=nr*nc,c(metht),prob = runif(nr*nc)),nr,nc)
#' r1 <- GRanges(rep('chr1', nr), IRanges(1:nr, width=1), strand='*')
#' names(r1) <- 1:nr
#' cd1 <- DataFrame(Group=rep(c('G1','G2'),each=nc/2),row.names=LETTERS[1:nc])
#' OBJ1 <- cBSData(rowRanges=r1,methReads=methc,totalReads=metht,colData=cd1)
#' OBJ2 <- methHMEM(OBJ1, MaxK=2, mc.cores=2)
#' OBJ3 <- methHMMCMC(OBJ2, mc.cores=2)
#' OBJ4 <- findDMCs(OBJ3, mc.cores=2)
#' head(metadata(OBJ4)$DMCHMM)
#' @exportMethod findDMCs
setGeneric("findDMCs", function(object, formula, FDRthreshold, Methylthreshold,
    mc.cores, windowsize, weightfunction) standardGeneric("findDMCs"))

#' @title qqDMCs method
#' @description Creates a Q-Q plot based on the p-values obtained
#' from \code{\link{findDMCs}} method
#' @author Farhad Shokoohi <shokoohi@icloud.com>
#' @name qqDMCs-method
#' @inheritParams params
#' @return A QQ plot
#' @importFrom stats ppoints
#' @import S4Vectors
#' @import graphics
#' @examples
#' set.seed(1980)
#' nr <- 150; nc <- 8
#' metht <- matrix(as.integer(runif(nr * nc, 0, 100)), nr)
#' methc <- matrix(rbinom(n=nr*nc,c(metht),prob = runif(nr*nc)),nr,nc)
#' r1 <- GRanges(rep('chr1', nr), IRanges(1:nr, width=1), strand='*')
#' names(r1) <- 1:nr
#' cd1 <- DataFrame(Group=rep(c('G1','G2'),each=nc/2),row.names=LETTERS[1:nc])
#' OBJ1 <- cBSData(rowRanges=r1,methReads=methc,totalReads=metht,colData=cd1)
#' OBJ2 <- methHMEM(OBJ1, MaxK=2, mc.cores=2)
#' OBJ3 <- methHMMCMC(OBJ2, mc.cores=2)
#' OBJ4 <- findDMCs(OBJ3, mc.cores=2)
#' qqDMCs(OBJ4)
#' @exportMethod qqDMCs
setGeneric("qqDMCs", function(object, ...) standardGeneric("qqDMCs"))

#' @title manhattanDMCs method
#' @description Creates a Manhattan plot based on the p-values obtained
#' from \code{\link{findDMCs}} method
#' @author Farhad Shokoohi <shokoohi@icloud.com>
#' @name manhattanDMCs-method
#' @inheritParams params
#' @return A Manhattan plot
#' @import S4Vectors
#' @import graphics
#' @importFrom calibrate textxy
#' @examples
#' set.seed(1980)
#' nr <- 150; nc <- 8
#' metht <- matrix(as.integer(runif(nr * nc, 0, 100)), nr)
#' methc <- matrix(rbinom(n=nr*nc,c(metht),prob = runif(nr*nc)),nr,nc)
#' r1 <- GRanges(rep('chr1', nr), IRanges(1:nr, width=1), strand='*')
#' names(r1) <- 1:nr
#' cd1 <- DataFrame(Group=rep(c('G1','G2'),each=nc/2),row.names=LETTERS[1:nc])
#' OBJ1 <- cBSData(rowRanges=r1,methReads=methc,totalReads=metht,colData=cd1)
#' OBJ2 <- methHMEM(OBJ1, MaxK=2, mc.cores=2)
#' OBJ3 <- methHMMCMC(OBJ2, mc.cores=2)
#' OBJ4 <- findDMCs(OBJ3, mc.cores=2)
#' manhattanDMCs(OBJ4)
#' @keywords visualization manhattan
#' @exportMethod manhattanDMCs
setGeneric("manhattanDMCs", function(object, col, chrlabs, suggestiveline,
    genomewideline, highlight, logp, annotatePval, annotateTop, ...)
    standardGeneric("manhattanDMCs"))
