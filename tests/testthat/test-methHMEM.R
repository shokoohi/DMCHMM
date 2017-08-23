test_that("test methHMEM", {
    set.seed(1980)
    nr <- 150; nc <- 8
    metht <- matrix(as.integer(runif(nr * nc, 0, 100)), nr)
    methc <- matrix(rbinom(n=nr*nc,c(metht),prob = runif(nr*nc)),nr,nc)
    r1 <- GRanges(rep('chr1', nr), IRanges(1:nr, width=1), strand='*')
    names(r1) <- 1:nr
    cd1 <- DataFrame(Group=rep(c('G1','G2'),each=nc/2),row.names=LETTERS[1:nc])
    OBJ1 <- cBSData(rowRanges=r1,methReads=methc,totalReads=metht,colData=cd1)
    OBJ2 <- methHMEM(OBJ1, MaxK=2, mc.cores=2)
    expect_equal(metadata(OBJ2)$K, c(3, 3, 3, 3, 3, 3, 3, 3))
})
