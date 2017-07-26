## ---- eval=FALSE, results='hide', message=FALSE--------------------------
#  library(DMCHMM)
#  fn <- list.files(system.file("extdata",package = "DMCHMM"))
#  fn.f <- list.files(system.file("extdata",package="DMCHMM"), full.names=TRUE)
#  OBJ <- readBismark(fn.f, fn)
#  cdOBJ <- DataFrame(Cell = factor(c("BC", "TC","Mono"),
#  labels = c("BC", "TC", "Mono")), row.names = c("BCU1568","BCU173","BCU551"))
#  colData(OBJ) <- cdOBJ
#  OBJ

## ---- eval=FALSE, results='hide', message=FALSE--------------------------
#  nr <- 500; nc <- 8
#  metht <- matrix(as.integer(runif(nr * nc, 0, 100)), nr)
#  methc <- matrix(rbinom(n=nr*nc,c(metht),prob = runif(nr*nc)),nr,nc)
#  r1 <- GRanges(rep("chr1", nr), IRanges(1:nr, width=1), strand="*")
#  names(r1) <- 1:nr
#  cd1 <- DataFrame(Group=rep(c("G1","G2"),each=nc/2),row.names=LETTERS[1:nc])
#  OBJ1 <- BSData(rowRanges=r1,methReads=methc,totalReads=metht,colData=cd1)
#  OBJ1

## ---- eval=FALSE, results='hide', message=FALSE--------------------------
#  OBJ2 <- methHMEM(OBJ1, mc.cores=detectCores())
#  OBJ2

## ---- eval=FALSE, results='hide', message=FALSE--------------------------
#  OBJ3 <- methHMMCMC(OBJ2, mc.cores=detectCores())
#  OBJ3

## ---- eval=FALSE, results='hide', message=FALSE--------------------------
#  OBJ4 <- findDMCs(OBJ3, mc.cores=detectCores())
#  head(metadata(OBJ4)$DMCHMM)

