.marglike <- function(beta.est, Pm.est, Kv, Yv, nv, useweight) {
    npos = length(Yv)

    pX.filter = array(0, c(npos, Kv + 2))
    dimnames(pX.filter) = list(index = seq_len(npos),
                            states = c(seq_len(Kv + 2)))

    pX.smooth = array(0, c(npos, Kv + 2))
    dimnames(pX.smooth) = list(index = seq_len(npos),
                            states = c(seq_len(Kv + 2)))

    pX.predict = array(0, c(npos, Kv + 2))
    dimnames(pX.predict) = list(index = seq_len(npos),
                            states = c(seq_len(Kv + 2)))

    Xfilter = rep(0, npos)
    Xsmooth = rep(0, npos)

    pv <- beta.est

    Pm <- Pm.est

    if (useweight) {
        nw <- log(nv + 2)
    } else {
        nw <- rep(1, length(Yv))
    }

    Ptmp <- Pm
    for (i in seq_len(100)) {
        Ptmp <- Ptmp %*% Pm
    }

    pX.predict[1, ] <- Ptmp[1, ]  #  Pi
    pmarg <- sum(dbinom(Yv[1], nv[1], pv) * pX.predict[1, ])
    pX.filter[1, ] <- c(dbinom(Yv[1], nv[1], pv) * pX.predict[1, ])/pmarg
    # filtering

    pX.predict[1, ] <- t(t(Pm) %*% pX.filter[1, ])
    # prediction
    pmarg <- sum(dbinom(Yv[1], nv[1], pv) * pX.predict[1, ])
    # marginal
    pX.filter[1, ] <- c(dbinom(Yv[1], nv[1], pv) * pX.predict[1, ])/pmarg

    mlike <- 0
    mlike <- mlike + nw[1] * log(pmarg)

    for (i in 2:npos) {
        pX.predict[i, ] <- t(t(Pm) %*% pX.filter[i - 1, ])
        pmarg <- sum(dbinom(Yv[i], nv[i], pv) * pX.predict[i, ])
        pX.filter[i, ] <- c(dbinom(Yv[i], nv[i], pv) * pX.predict[i, ])/pmarg
        mlike <- mlike + nw[i] * log(pmarg)
    }
    return(mlike)
}

.marglike.M <- function(X.estimate, Yv, nv, useweight) {
    npos = length(Yv)
    if (useweight) {
        nw <- log(nv + 2)
    } else {
        nw <- rep(1, length(Yv))
    }
    mlike <- sum(nw * dbinom(Yv, nv, X.estimate, log=TRUE), na.rm = TRUE)
    return(mlike)
}

.HMMmcmc.M <- function(Yv, nv, Kv, old.prob, old.Uvec, old.P.U, nburn,
    nthin, nsamp, useweight) {
    nits <- nburn + nthin * nsamp
    ico <- 0
    X.estimate <- matrix(0, nrow = nsamp, ncol = length(Yv))
    p.samp <- matrix(0, nrow = nsamp, ncol = Kv + 2)
    P.U.samp <- array(0, c(nsamp, Kv + 2, Kv + 2))
    Prob.U.vec <- matrix(0, nrow = length(Yv), ncol = Kv + 2)
    x.samp <- matrix(0, nrow = nsamp, ncol = length(Yv))
    index <- seq_len(length(Yv))
    n = length(Yv)
    if (useweight) {
        nw <- log(nv + 2)
    } else {
        nw <- rep(1, length(Yv))
    }
    wYv = nw * Yv
    wnv = nw * nv

    for (Itr in seq_len(nits)) {

        if (Itr == nburn) {
            Prob.U.vec <- Prob.U.vec * 0
        }
        dtmp1 <- dbinom(Yv[1], nv[1], old.prob, log = TRUE)
        if(any(is.infinite(dtmp1))){
        dtmp1 = dbinom(Yv[1], nv[1], rep(1/length(dtmp1), length(dtmp1)),
            log = TRUE)
        }
        dtmp3 <- log(old.P.U[, old.Uvec[2] + 1])
        if(any(is.infinite(dtmp3))){
            dtmp3 = log(rep(1/length(dtmp3), length(dtmp3)))
        }

        dtmp <- dtmp1 + dtmp3
        ptmp <- exp(dtmp - max(dtmp))

        old.Uvec[1] <- sample(0:(Kv + 1), size = 1, prob = ptmp)
        Prob.U.vec[1, ] <- Prob.U.vec[1, ] + ptmp/sum(ptmp)

        for (i in 2:(n - 1)) {
            dtmp1 <- dbinom(Yv[i], nv[i], old.prob, log = TRUE)
            if(any(is.infinite(dtmp1))){
            dtmp1 = dbinom(Yv[i], nv[i], rep(1/length(old.prob),
                length(old.prob)), log = TRUE)
            }

            dtmp3 <- log(old.P.U[old.Uvec[i - 1] + 1, ] * old.P.U[, old.Uvec[i +
                1] + 1])
            if(any(is.infinite(dtmp3))){
                dtmp3= log(rep(1/length(dtmp3), length(dtmp3)))
            }

            dtmp <- dtmp1 + dtmp3
            ptmp <- exp(dtmp - max(dtmp))
            old.Uvec[i] <- sample(0:(Kv + 1), size = 1, prob = ptmp)
            Prob.U.vec[i, ] <- Prob.U.vec[i, ] + ptmp/sum(ptmp)
        }
        dtmp1 <- dbinom(Yv[n], nv[n], old.prob, log = TRUE)
        if(any(is.infinite(dtmp1))){
            dtmp1 = dbinom(Yv[n], nv[n], rep(1/length(old.prob),
                length(old.prob)), log = TRUE)
        }

        dtmp3 <- log(old.P.U[, old.Uvec[n - 1] + 1])
        if(any(is.infinite(dtmp3))){
            dtmp3 = log(rep(1/length(dtmp3), length(dtmp3)))
        }

        dtmp <- dtmp1 + dtmp3
        ptmp <- exp(dtmp - max(dtmp))

        old.Uvec[n] <- sample(0:(Kv + 1), size = 1, prob = ptmp)
        Prob.U.vec[n, ] <- Prob.U.vec[n, ] + ptmp/sum(ptmp)

        for (k in 0:(Kv + 1)) {
            kpos <- index[old.Uvec == k]
            tvec <- table(old.Uvec[kpos + 1])
            states <- as.numeric(names(tvec))
            tvec <- as.numeric(tvec)
            old.P.U[k + 1, states + 1] <- tvec/sum(tvec)
        }

        # Now do the old.probs

        for (k in seq_len(Kv)) {
            kpos1 <- index[old.Uvec == k]
            ytot <- sum(wYv[kpos1])
            ntot <- sum(wnv[kpos1])
            FL <- pbeta(old.prob[k], ytot + 1/2, ntot - ytot + 1/2)
            FU <- pbeta(old.prob[k + 2], ytot + 1/2, ntot - ytot + 1/2)
            u <- runif(1)
            old.prob[k + 1] <- qbeta(FL + u * (FU - FL), ytot + 1/2, ntot -
                ytot + 1/2)
        }
        if (Itr > nburn & Itr%%nthin == 0) {
            ico <- ico + 1
            p.samp[ico, ] <- old.prob
            P.U.samp[ico, , ] <- old.P.U
            x.samp[ico, ] <- old.Uvec
            X.estimate[ico, ] = c((Prob.U.vec/rowSums(Prob.U.vec)) %*%
                t(old.prob))
        }
    }

    return(list(p.samp = p.samp, P.U.samp = P.U.samp, Prob.U.vec = Prob.U.vec,
        x.samp = x.samp, X.estimate = X.estimate))
}

.runMCMC.M <- function(Y.tot, n.tot, K, EM.old.prob, EM.old.Uvec, EM.old.P.U,
    nburn, nthin, nsamp, useweight) {

    xMCMC <- .HMMmcmc.M(Yv = Y.tot, nv = n.tot, Kv =K, EM.old.prob, EM.old.Uvec,
        EM.old.P.U, nburn, nthin, nsamp, useweight)

    p.ests <- apply(xMCMC$p.samp, 2, mean)
    p.var.ests <- apply(xMCMC$p.samp, 2, var)

    p.U.ests <- apply(xMCMC$P.U.samp, 2:3, mean)
    p.U.var.ests <- apply(xMCMC$P.U.samp, 2:3, var)

    logMlik <- .marglike(p.ests, p.U.ests, K, Y.tot, n.tot, useweight)
    BIC <- -2 * logMlik + log(length(Y.tot)) * (K + (K + 2) * (K + 1))

    X.estimate <- apply(xMCMC$X.estimate, 2, mean)
    X.estimate.var <- apply(xMCMC$X.estimate, 2, var)

    logMlik.M <- .marglike.M(X.estimate, Y.tot, n.tot, useweight)
    BIC.M <- -2 * logMlik.M + log(length(Y.tot)) * (K + (K + 2) * (K + 1))

    return(list(xMCMC = xMCMC, p.ests = p.ests, p.var.ests = p.var.ests,
        p.U.ests = p.U.ests, p.U.var.ests = p.U.var.ests, logMlik = logMlik,
        BIC = BIC, K =K, X.estimate=X.estimate, X.estimate.var = X.estimate.var,
        logMlik.M = logMlik.M, BIC.M = BIC.M))
}

.methHMMCMC <- function(object, useweight, nburn, nthin, nsamp, mc.cores) {

    if (missing(object) | is(object)[1] != "BSDMCs")
        stop("A BSDMCs object must be provided.")

    if (missing(useweight)) {
        useweight = TRUE
    } else if (!is.logical(useweight)) {
        stop("A logical value must be provided for useweight.")
    }

    if (missing(nburn)) {
        nburn = 500
    } else if (!is.numeric(nburn) | nburn < 10) {
        stop("An integer value greated than 10 must be provided for nburn.")
    }
    nburn = as.integer(nburn)

    if (missing(nthin)) {
        nthin = 1
    } else if (!is.numeric(nthin) | nthin < 1) {
        stop("An integer value greated than 1 must be provided for nthin.")
    }
    nthin = as.integer(nthin)

    if (missing(nsamp)) {
        nsamp = 500
    } else if (!is.numeric(nsamp) | nsamp < 300) {
        stop("An integer value greated than 300 must be provided for nsamp.")
    }
    nsamp = as.integer(nsamp)

    if (missing(mc.cores)) {
        mc.cores = multicoreWorkers()
    } else if (!is.numeric(mc.cores) | mc.cores <= 0)
    {
        stop("An integer value greater than 0 must be provided for mc.cores")
    }
    mc.cores = as.integer(mc.cores)

    strand(object) <- "*"
    object <- sort(object)

    nGroup = length(unique(colData(object))[[1]])
    if (nGroup <= 1)
        stop("Assign groups using colData.")

    nPos = nrow(object)
    nSam = length(colnames(object))

    old.object = list(methReads = assays(object)$methReads,
                    totalReads = assays(object)$totalReads,
                    methStates = assays(object)$methStates,
                    Beta = metadata(object)$Beta,
                    K = metadata(object)$K,
                    Pm = metadata(object)$Pm)

    optbp <- MulticoreParam(workers = mc.cores, progressbar = TRUE)
    .mfun2 <- function(j){
        .runMCMC.M(c(old.object$methReads[, j]),
                c(old.object$totalReads[, j]), c(old.object$K[j]),
                old.object$Beta[[j]], c(old.object$methStates[,j]),
                old.object$Pm[[j]], nburn, nthin, nsamp, useweight)
    }
    totres2 <- bplapply(seq_len(nSam), .mfun2, BPPARAM = optbp)

    Khat <- vapply(totres2, function(elt) elt$K, numeric(1))
    Betahat <- lapply(totres2, function(elt) elt$p.ests)
    Phat <- lapply(totres2, function(elt) elt$p.U.ests)
    methLevels.new <- vapply(totres2, function(elt) elt$X.estimate,
                                numeric(nPos))
    methVars.new <- vapply(totres2, function(elt) elt$X.estimate.var,
                            numeric(nPos))
    methStates.new <- vapply(totres2, function(elt) apply(elt$xMCMC$Prob.U.vec,
                                1, which.max) - 1, numeric(nPos))
    storage.mode(methStates.new) <- "integer"

    colnames(methLevels.new) <- colnames(object)
    rownames(methLevels.new) <- NULL
    colnames(methVars.new) <- colnames(object)
    rownames(methVars.new) <- NULL
    colnames(methStates.new) <- colnames(object)
    rownames(methStates.new) <- NULL

    predictedMeth <- cBSDMCs(methReads = methReads(object),
                            totalReads = totalReads(object),
                            methLevels = methLevels.new,
                            methVars = methVars.new,
                            methStates = apply(methStates.new, 2, as.integer),
                            rowRanges = rowRanges(object),
                            colData = colData(object))
    metadata(predictedMeth)$K = Khat
    metadata(predictedMeth)$Beta = Betahat
    metadata(predictedMeth)$Pm = Phat

    return(predictedMeth)
}

#' @rdname methHMMCMC-method
#' @aliases methHMMCMC-method methHMMCMC
setMethod("methHMMCMC", signature = c(object = "BSDMCs", useweight = "ANY",
    nburn = "ANY", nthin = "ANY", nsamp = "ANY", mc.cores = "ANY"), .methHMMCMC)

