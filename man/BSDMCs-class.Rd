% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllClasses.R
\docType{class}
\name{BSDMCs-class}
\alias{BSDMCs-class}
\alias{BSDMCs}
\title{BSDMCs object}
\arguments{
\item{methReads}{The matrix \code{methReads} contains the number of
methylated reads spanning a CpG-site. The rows represent the CpG sites in
\code{rowRanges} and the columns represent the samples in \code{colData}.}

\item{totalReads}{The matrix \code{totalReads} contains the number of reads
spanning a CpG-site. The rows represent the CpG sites in \code{rowRanges}
and the columns represent the samples in \code{colData}.}

\item{methLevels}{The matrix \code{methLevels} contains the predicted
methylation level spanning a CpG-site using Hidden Markov model. The rows
represent the CpG sites in \code{rowRanges} and the columns represent the
samples in \code{colData}.}

\item{methStates}{The matrix \code{methStates} contains the state of
methylation obtained from Hidden Markov model spanning a CpG-site. The rows
represent the CpG sites in \code{rowRanges} and the columns represent the
samples in \code{colData}. The value of state is stored in \code{metadata},
named \code{Beta}.}

\item{methVars}{The matrix \code{methVars} contains the variances of the
corresponding \code{methLevels} obtianed from MCMC.}
}
\value{
A \code{\link{BSDMCs-class}} object
}
\description{
The \code{BSDMCs} object is an S4 class that represents
differentially methylated CpG sites (DMCs) in BS-Seq Data.
}
\section{Slots}{

\describe{
\item{\code{methReads}}{An integer matrix}

\item{\code{totalReads}}{An integer matrix}

\item{\code{methLevels}}{A numeric matrix}

\item{\code{methStates}}{An integer matrix}

\item{\code{methVars}}{A double matrix}
}}

\examples{
nr <- 500; nc <- 16
metht <- matrix(as.integer(runif(nr * nc, 0, nr)), nr)
methc <- matrix(rbinom(n=nr*nc,c(metht),prob = runif(nr*nc)),nr,nc)
meths <- matrix(as.integer(runif(nr * nc, 0, 10)), nr)
methl <- methc/metht
methv <- matrix((runif(nr * nc, 0.1, 0.5)), nr)
r1 <- GRanges(rep('chr1', nr), IRanges(1:nr, width=1), strand='*')
names(r1) <- 1:nr
cd1 <- DataFrame(Group=rep(c('G1','G2'),each=nc/2),row.names=LETTERS[1:nc])
OBJ2 <- cBSDMCs(rowRanges=r1,methReads=methc,totalReads=metht,
methLevels=methl,methStates=meths,methVars=methv,colData=cd1)
OBJ2
}
\author{
Farhad Shokoohi <shokoohi@icloud.com>
}
\keyword{object}
