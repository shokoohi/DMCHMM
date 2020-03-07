## ---- eval=TRUE, message=FALSE------------------------------------------------
library(DMCHMM)
fn <- list.files(system.file("extdata",package = "DMCHMM"))
fn.f <- list.files(system.file("extdata",package="DMCHMM"), full.names=TRUE)
OBJ <- readBismark(fn.f, fn)
cdOBJ <- DataFrame(Cell = factor(c("BC", "TC","Mono"),
labels = c("BC", "TC", "Mono")), row.names = c("BCU1568","BCU173","BCU551"))
colData(OBJ) <- cdOBJ
OBJ

## ---- eval=TRUE, message=FALSE------------------------------------------------
nr <- 150; nc <- 8
metht <- matrix(as.integer(runif(nr * nc, 0, 20)), nr)
methc <- matrix(rbinom(n=nr*nc,c(metht),prob = runif(nr*nc)),nr,nc)
r1 <- GRanges(rep("chr1", nr), IRanges(1:nr, width=1), strand="*")
names(r1) <- 1:nr
cd1 <- DataFrame(Group=rep(c("G1","G2"),each=nc/2),row.names=LETTERS[1:nc])
OBJ1 <- cBSData(rowRanges=r1,methReads=methc,totalReads=metht,colData=cd1)
OBJ1

## ---- eval=TRUE, message=FALSE------------------------------------------------
OBJ2 <- methHMEM(OBJ1, MaxK=2)
OBJ2

## ---- eval=TRUE, message=FALSE------------------------------------------------
OBJ3 <- methHMMCMC(OBJ2)
OBJ3

## ---- eval=TRUE, message=FALSE------------------------------------------------
OBJ4 <- findDMCs(OBJ3)
head(metadata(OBJ4)$DMCHMM)

