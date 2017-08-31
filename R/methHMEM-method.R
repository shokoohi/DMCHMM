.EMHMM.M <- function(xv, Kv, Yv, nv, MaxEmiter, epsEM, useweight) {
    if (useweight) {
        nw <- log(nv + 2)
    } else {
        nw <- rep(1, length(Yv))
    }
    nPos = length(Yv)

    pX.filter = array(0, c(nPos, Kv + 2))
    dimnames(pX.filter) = list(index = 1:nPos, states = c(1:(Kv + 2)))

    X.filter = array(0, c(nPos))
    dimnames(X.filter) = list(index = 1:nPos)

    pX.predict = array(0, c(nPos, Kv + 2))
    dimnames(pX.predict) = list(index = 1:nPos, states = c(1:(Kv + 2)))

    pX.smooth = array(0, c(nPos, Kv + 2))
    dimnames(pX.smooth) = list(index = 1:nPos, states = c(1:(Kv + 2)))

    X.smooth = array(0, c(nPos))
    dimnames(X.smooth) = list(index = 1:nPos)

    W01k.Estep = array(0, c(Kv + 2))
    dimnames(W01k.Estep) = list(states = c(1:(Kv + 2)))

    Wik.Estep = array(0, c(nPos, Kv + 2))
    dimnames(Wik.Estep) = list(index = 1:nPos, states = c(1:(Kv + 2)))

    Wijk.Estep = array(0, c(nPos - 1, Kv + 2, Kv + 2))
    dimnames(Wijk.Estep) = list(index = 2:nPos, states = c(1:(Kv + 2)),
        states = c(1:(Kv + 2)))

    pv <- exp(c(xv[1:Kv], 0))
    pv <- pv/sum(pv)
    pv <- c(0, cumsum(pv))

    Pm <- matrix(0, nrow = Kv + 2, ncol = Kv + 2)
    Pm[, 1:(Kv + 1)] <- matrix(xv[c((Kv + 1):(Kv + (Kv + 2) * (Kv + 1)))],
        Kv + 2, Kv + 1)
    Pm <- t(apply(exp(Pm), 1, function(x) {
        return(x/sum(x))
    }))

    Pi0 <- Pm
    for (i in 1:100) Pi0 <- Pi0 %*% Pm
    Pi0 <- Pi0[1, ]

    diffpar = 1
    niter = 1

    while ((diffpar > epsEM) & (niter < MaxEmiter)) {
        # E-step
        pX.predict[1, ] <- t(t(Pm) %*% Pi0)
        # prediction
        pmarg <- sum(dbinom(Yv[1], nv[1], pv) * pX.predict[1, ])
        # marginal
        pX.filter[1, ] <- c(dbinom(Yv[1], nv[1], pv) * pX.predict[1, ])/pmarg

        for (i in 2:nPos) {
            pX.predict[i, ] <- t(t(Pm) %*% pX.filter[i - 1, ])
            pmarg <- sum(dbinom(Yv[i], nv[i], pv) * pX.predict[i, ])
            pmarg <- ifelse(pmarg < 1e-200, 1e-200, pmarg)
            pX.filter[i, ] <- c(dbinom(Yv[i], nv[i], pv) * pX.predict[i,
                ])/pmarg
        }
        pX.smooth[nPos, ] <- pX.filter[nPos, ]
        for (i in (nPos - 1):1) {
            hh <- (pX.smooth[i + 1, ]/pX.predict[i + 1, ])
            hh[is.na(hh)] <- 0
            pX.smooth[i, ] <- pX.filter[i, ] * (Pm %*% hh)
        }
        W01k.Estep = pX.smooth[1, ]
        Wik.Estep = pX.smooth

        HWijk.Estep = Wijk.Estep

        for (i in 2:nPos) {
            Wijk.Estep[i - 1, , ] = ((Pm * matrix((pX.filter[i - 1, ]),
                Kv + 2, Kv + 2)) * matrix(pX.smooth[i, ], Kv + 2, Kv +
                2, byrow = TRUE))/matrix(pX.predict[i, ], Kv + 2, Kv +
                2, byrow = TRUE)
        }
        Wijk.Estep[is.na(Wijk.Estep)] <- 0
        new.pv = ((nw * Yv) %*% Wik.Estep)/((nw * nv) %*% Wik.Estep)

        new.Pm <- apply(Wijk.Estep, 2:3, function(x) sum(x)) /
            apply(Wijk.Estep, 3, function(x) sum(x))

        # folowing lines will be activated if we are seeking an alternative
        # estimator for initial probabilities.
        tmppi <- new.Pm
        for (i in 1:100) tmppi = tmppi %*% new.Pm
        new.Pi0 = tmppi[1, ]

        diffpar = 0
        diffpar = diffpar + sum((Pi0 - new.Pi0)^2)
        diffpar = diffpar + sum((Pm - new.Pm)^2)
        diffpar = diffpar + sum((pv - new.pv)^2)

        niter = niter + 1

        Pi0 = new.Pi0
        Pm = new.Pm
        pv = new.pv
    }

    mlike <- 0
    pX.predict[1, ] <- t(t(Pm) %*% Pi0)
    pmarg <- sum(dbinom(Yv[1], nv[1], pv) * pX.predict[1, ])
    pX.filter[1, ] <- c(dbinom(Yv[1], nv[1], pv) * pX.predict[1, ])/pmarg
    X.filter[1] <- which.max(pX.filter[1, ])
    mlike <- mlike + nw[1] * log(pmarg)

    for (i in 2:nPos) {
        pX.predict[i, ] <- t(t(Pm) %*% pX.filter[i - 1, ])
        pmarg <- sum(dbinom(Yv[i], nv[i], pv) * pX.predict[i, ])
        pX.filter[i, ] <- c(dbinom(Yv[i], nv[i], pv) * pX.predict[i, ])/pmarg
        X.filter[i] <- which.max(pX.filter[i, ])
        mlike <- mlike + nw[i] * log(pmarg)
    }
    pX.smooth[nPos, ] <- pX.filter[nPos, ]
    X.smooth[nPos] <- which.max(pX.smooth[nPos, ])
    for (i in (nPos - 1):1) {
        hh = (pX.smooth[i + 1, ]/pX.predict[i + 1, ])
        hh[is.na(hh)] = 0
        pX.smooth[i, ] <- pX.filter[i, ] * (Pm %*% hh)
        X.smooth[i] <- which.max(pX.smooth[i, ])
    }
    X.estimate = rowSums(pX.smooth * matrix(pv, ncol = Kv + 2, nrow = nPos,
        byrow = TRUE))
    X.estimate[X.estimate > 1] = 1
    X.estimate[X.estimate < 0] = 0

    mlike <- sum(nw * dbinom(Yv, nv, X.estimate, log=TRUE))

    return(list(mlike = mlike, beta = pv, Pm = Pm, Pi0 = Pi0,
                X.smooth = X.smooth, X.filter = X.filter,
                X.estimate = X.estimate))
}

.runEM.M <- function(yy, nn, MaxK, MaxEmiter, epsEM, useweight) {
    K = 1
    xstart <- rep(0, K + (K + 2) * (K + 1) + (K + 1))
    xEM <- .EMHMM.M(xstart, K, yy, nn, MaxEmiter, epsEM, useweight = useweight)
    logMlik <- xEM$mlike
    BIC <- -2 * logMlik + log(length(yy)) * (K + (K + 2) * (K + 1) + (K +
        1))
    BICmax = 10^20
    if (!is.na(BIC)) {
        BICmax = BIC
        EstHMMSam1 = list(xEM = xEM, K = K, logMlik = logMlik, BIC = BIC)
    } else {
        BIC = 10^20
    }
    while ((K <= MaxK) & (BICmax >= BIC)) {
        K = K + 1
        xstart <- rep(0, K + (K + 2) * (K + 1) + (K + 1))
        xEM <- .EMHMM.M(xstart, K, yy, nn, MaxEmiter, epsEM, useweight)
        logMlik <- xEM$mlike
        BIC <- -2 * logMlik + log(length(yy)) * (K + (K + 2) * (K + 1) +
            (K + 1))
        if ((BIC < BICmax) & (!is.na(BIC))) {
            BICmax = BIC
            EstHMMSam1 = list(xEM = xEM, K = K, logMlik = logMlik, BIC = BIC)
        }
    }
    return(EstHMMSam1)
}

.methHMEM <- function(object, MaxK, MaxEmiter, epsEM, useweight, mc.cores) {

    if (missing(object) | class(object) != "BSData")
        stop("A BSData object must be provided.")

    if (missing(MaxK)) {
        MaxK = 10
    } else if (!is.numeric(MaxK) | MaxK < 1) {
        stop("An integer value greated than 0 must be provided for MaxK")
    }
    MaxK = as.integer(MaxK)

    if (missing(MaxEmiter)) {
        MaxEmiter = 300
    } else if (!is.numeric(MaxEmiter) | MaxEmiter < 10) {
        stop("An integer value greated than 10 must be provided for MaxEmiter.")
    }
    MaxEmiter = as.integer(MaxEmiter)

    if (missing(epsEM)) {
        epsEM = 1e-05
    } else if (!is.numeric(epsEM) | epsEM <= 0 | epsEM >= 0.1) {
        stop("A numeric value between 0 and 0.1 must be provided for epsEM.")
    }

    if (missing(useweight)) {
        useweight = TRUE
    } else if (!is.logical(useweight)) {
        stop("A logical value must be provided for useweight.")
    }
    useweight = as.logical(useweight)

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

    optbp <- MulticoreParam(workers = mc.cores, progressbar = TRUE)
    .mfun <- function(itr){
        .runEM.M(c(assays(object)$methReads[, itr]),
                c(assays(object)$totalReads[, itr]), MaxK, MaxEmiter, epsEM,
                useweight)
    }
    totres <- bplapply(seq_len(nSam), .mfun, BPPARAM = optbp)

    Khat <- vapply(totres, function(elt) elt$K, numeric(1))
    Betahat <- lapply(totres, function(elt) elt$xEM$beta)
    Phat <- lapply(totres, function(elt) elt$xEM$Pm)
    methLevels.new <- vapply(totres, function(elt) elt$xEM$X.estimate,
                            numeric(nPos))
    methStates.new <- vapply(totres, function(elt) elt$xEM$X.smooth - 1,
                            numeric(nPos))
    storage.mode(methStates.new) <- "integer"

    colnames(methLevels.new) <- colnames(object)
    rownames(methLevels.new) <- NULL
    colnames(methStates.new) <- colnames(object)
    rownames(methStates.new) <- NULL

    predictedMeth <- cBSDMCs(methReads = methReads(object),
                            totalReads = totalReads(object),
                            methStates = apply(methStates.new, 2, as.integer),
                            methLevels = methLevels.new,
                            colData = colData(object),
                            rowRanges = rowRanges(object))
    metadata(predictedMeth)$K = Khat
    metadata(predictedMeth)$Beta = Betahat
    metadata(predictedMeth)$Pm = Phat

    return(predictedMeth)
}

#' @rdname methHMEM-method
#' @aliases methHMEM-method methHMEM
setMethod("methHMEM", signature = c(object = "BSData", MaxK = "ANY",
    MaxEmiter = "ANY", epsEM = "ANY", useweight = "ANY", mc.cores = "ANY"),
    .methHMEM)

