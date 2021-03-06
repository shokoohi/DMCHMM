% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllClasses.R
\docType{class}
\name{BSData-class}
\alias{BSData-class}
\alias{BSData}
\title{BSData object}
\arguments{
\item{methReads}{The matrix \code{methReads} contains the number of
methylated reads spanning a CpG-site. The rows represent the CpG sites in
\code{rowRanges} and the columns represent the samples in \code{colData}.}

\item{totalReads}{The matrix \code{totalReads} contains the number of reads
spanning a CpG-site. The rows represent the CpG sites in \code{rowRanges}
and the columns represent the samples in \code{colData}.}
}
\value{
A \code{\link{BSData-class}} object
}
\description{
The \code{BSData} object is an S4 class that represents BS-Seq
Data.
}
\section{Slots}{

\describe{
\item{\code{methReads}}{An integer matrix}

\item{\code{totalReads}}{An integer matrix}
}}

\examples{
nr <- 500; nc <- 16
metht<-matrix(as.integer(runif(nr * nc, 0, nr)), nr)
methc<-matrix(rbinom(n=nr*nc,c(metht),prob = runif(nr*nc)),nr,nc)
r1 <- GRanges(rep("chr1", nr), IRanges(1:nr, width=1), strand="*")
names(r1) <- 1:nr
cd1<-DataFrame(Group=rep(c("G1","G2"),each=nc/2),row.names=LETTERS[1:nc])
OBJ1<-cBSData(rowRanges=r1,methReads=methc,totalReads=metht,colData=cd1)
OBJ1
}
\seealso{
\code{\link{SummarizedExperiment}} objects.
}
\author{
Farhad Shokoohi <shokoohi@icloud.com>
}
\keyword{object}
