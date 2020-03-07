.findDMCs <- function(object, formula, FDRthreshold, Methylthreshold, mc.cores,
    weightfunction){

    if (missing(object) | is(object)[1] != "BSDMCs")
        stop("A BSDMCs object must be provided.")

    if (missing(formula)) {
        formula = as.formula(paste("Methylation",
            paste(c(names(colData(object))), collapse = " + "), sep = " ~ "))
    } else if (is(formula)[1] == "formula") {
        stop("Either provide a formula or leave it blank.")
    } else {
        if (!all(all.vars(formula) %in% names(colData(object)))) {
            stop("The variables that are provided in the formula do no match
            any variable in", paste(names(colData(object)), collpase=", "), ".")
        } else {
            formula = as.formula(paste("Methylation ~ ",
                paste(all.vars(formula)[all.vars(formula) %in%
    names(colData(object))], collapse = " + ")))
        }
    }

    if (missing(FDRthreshold)) {
        FDRthreshold = 0.05
    } else if (!is.numeric(FDRthreshold) | FDRthreshold <=
        0 | FDRthreshold > 0.1) {
        stop("A numeric value between 0 and 0.1 must be
    provided for FDRthreshold.")
    }
    FDRthreshold = as.numeric(FDRthreshold)

    if (missing(Methylthreshold)) {
        Methylthreshold = 0.1
    } else if (!is.numeric(Methylthreshold) | Methylthreshold <
        0 | Methylthreshold > 0.9) {
        stop("A numeric value between 0 and 0.9 must be provided for
    Methylthreshold.")
    }
    Methylthreshold = as.numeric(Methylthreshold)

    if (missing(mc.cores)) {
        mc.cores = multicoreWorkers()
    } else if (!is.numeric(mc.cores) | mc.cores <= 0)
    {
        stop("An integer value greater than 0 must be provided for mc.cores.")
    }
    mc.cores = as.integer(mc.cores)

    if (missing(weightfunction)) {
        weightfunction = NULL
    } else if (!is.function(weightfunction)){
        stop("An function must be provided to create weights based on variances
        obtained form the MCMC algorithm.")
    }
    strand(object) <- "*"
    object <- sort(object)

    object.df <- as.data.frame(colData(object))
    object.df <- object.df[all.vars(formula)[all.vars(formula) %in%
                                            names(colData(object))]]
    ind1 <- vapply(object.df, is.character, logical(1))
    object.df[ind1] <- lapply(object.df[ind1], as.factor)
    ind1 <- vapply(object.df, is.factor, logical(1))
    object.df.factor <- as.data.frame(object.df[, ind1])
    names(object.df.factor) <- names(object.df)[ind1]
    ind1 <- vapply(object.df.factor, nlevels, numeric(1))>1
    fnames <- names(object.df.factor)[ind1]
    object.df.factor <- as.data.frame(object.df.factor[,ind1])
    names(object.df.factor) <- fnames

    if (ncol(object.df) < 1)
        stop("There is no covarite in the data.")

    nPos = nrow(object)
    nSam = length(colnames(object))
    methlist <- list(methLevels = assays(object)$methLevels,
                    methVars = assays(object)$methVars)
    Weights = c()

    optbp <- MulticoreParam(workers = mc.cores, progressbar = TRUE)
    .mfun3 <- function(itr){
        dat = data.frame(Methylation = c(methlist$methLevels[itr, ]),
                        Weights = c(methlist$methVars[itr, ]))
        dat = cbind(dat, object.df)
        dat$Methylation[dat$Methylation < 1e-10] <- 1e-10
        dat$Methylation[dat$Methylation > 0.9999999999] <- 0.9999999999
        dat$Methylation[is.na(dat$Methylation)] <- mean(dat$Methylation,
                                                        na.rm = TRUE)
        dat$Methylation <- log(dat$Methylation/(1 - dat$Methylation))
        dat$Methylation[is.na(dat$Methylation)] <- mean(dat$Methylation,
                                                        na.rm = TRUE)
        if(!is.null(weightfunction)){
            dat$Weights = weightfunction(c(methlist$methVars[itr, ]))
            options(show.error.messages = FALSE)
            suppressWarnings(lmodel <- try(lm(formula = formula,
                                            weights = Weights,
                                            data = dat), silent = TRUE))
            options(show.error.messages = TRUE)
        }else{
            options(show.error.messages = FALSE)
            suppressWarnings(lmodel <- try(lm(formula = formula,
                                            data = dat), silent = TRUE))
            options(show.error.messages = TRUE)
        }
        if ((is(lmodel)[1] == "try-error")) {
            p.val <- NA
            dmcs <- NA
            methDir <- NA
            return(list(p.val = p.val, dmcs = dmcs, methDir = methDir))
        } else {
            fst <- summary(lmodel)$fstatistic
            p.val <- pf(fst[1], fst[2], fst[3], lower.tail = FALSE)
            attributes(p.val) <- NULL

            if (dim(object.df.factor)[2] == 0) {
                return(list(p.val = p.val))
            } else {
                dmcs <- c()
                methDir <- c()
                if ("(Intercept)" %in% names(lmodel$coefficients)) {
                    baseL <- lmodel$coefficients[1]
                } else {
                    baseL <- lmodel$coefficients[c(grep(names(
                        object.df.factor)[1],
                        names(lmodel$coefficients), value = TRUE))][1]
                }

                ITR1 <- 0
                for (g in names(object.df.factor)) {
                    nGroup = length(levels(object.df.factor[, g]))
                    if (nGroup > 1) {
                        if ("(Intercept)" %in% names(lmodel$coefficients) &
                            ITR1 == 0) {
                            lc <- lmodel$coefficients[c("(Intercept)",
                                grep(g, names(lmodel$coefficients),
                                    value = TRUE))]
                            lc[-1] = lc[-1] + baseL
                        } else if (!("(Intercept)" %in%
                                names(lmodel$coefficients)) & ITR1 == 0) {
                            lc <- lmodel$coefficients[c(grep(g,
                                    names(lmodel$coefficients), value = TRUE))]
                        } else {
                            lc <- lmodel$coefficients[c(grep(g,
                                    names(lmodel$coefficients), value = TRUE))]
                            lc <- c(baseL, baseL + lc)
                        }
                        d <- (outer(1/(1 + exp(-lc)),
                                    1/(1 + exp(-lc)), "-"))
                        d = d[upper.tri(d)]
                        temp1 = d
                        temp1[d > Methylthreshold] = "hyper"
                        temp1[d < Methylthreshold] = "hypo"
                        temp1[d >= (-Methylthreshold) &
                                d <= Methylthreshold] = "equal"
                        methDir <- c(methDir, temp1)
                        bb <- NULL
                        bb$g <- "Tukey"
                        names(bb) <- g
                        if (!("(Intercept)" %in% names(lmodel$coefficients)) &
                            ITR1 == 0) {
                            attr(bb, "interaction_average") <- TRUE
                        } else {
                            attr(bb, "interaction_average") <- FALSE
                        }
                        attr(bb, "covariate_average") <- FALSE
                        attr(bb, "class") <- "mcp"

                        ITR1 <- 1

                        options(show.error.messages = FALSE)
                        suppressWarnings(lmodelc <- try(glht(lmodel,
                                                linfct = bb, silent = TRUE)))
                        options(show.error.messages = TRUE)
                        if ((is(lmodelc)[1] == "try-error")) {
                            dmcs <- c(dmcs, rep(NA, nGroup))
                        } else {
                            hhh <- (summary(lmodelc)$test$pvalues)
                            attributes(hhh) <- NULL
                            dmcs <- c(dmcs, hhh)
                        }
                    }
                }
                dmcs <- (dmcs < 0.05) * 1
                return(list(p.val = p.val, dmcs = dmcs,
                            methDir = methDir))
            }
        }
    }
    totres3 <- bplapply(seq_len(nPos), .mfun3, BPPARAM = optbp)

    pval <- vapply(totres3, function(elt) elt$p.val, numeric(1))

    fdr.fun <- function(x) {
        (fdrtool(x, statistic = "pvalue", plot = FALSE,
            verbose = FALSE))$qval
    }
    if (nPos < 1000) {
        qval <- fdr.fun(pval)
    } else {
        new.pval <- split(pval, ceiling(seq_along(pval)/500))
        aa = length(new.pval)
        if(length(new.pval[[aa]])<300){
            new.pval[[aa-1]] <- c(new.pval[[aa-1]] , new.pval[[aa]])
            new.pval = new.pval[-aa]
        }
        qval <- as.vector(unlist(lapply(new.pval, fdr.fun)))
    }

    DMCs = (qval < FDRthreshold) * 1

    nfactor <- dim(object.df.factor)[2]
    namesPairDMC <- namesPairDIR <- NULL
    if (nfactor > 0) {
        nGroupS <- unlist(vapply(object.df.factor,
                                function (x) {length(levels(x))}, numeric(1)))
        nPairs <- vapply(nGroupS, function (x) {ncol(combn(seq_len(x), 2))},
                        numeric(1))
        np <- unlist(vapply(object.df.factor,
                function(x) {paste(combn(levels(x), 2)[1, ],
                    combn(levels(x), 2)[2, ], sep = "vs")}, character(1)))
        names(np) <- NULL
        namesPairDMC <- c(paste("DMCs", rep(names(object.df.factor),
                                            nPairs), np, sep = ""))
        namesPairDIR <- c(paste("methDir", rep(names(object.df.factor),
                                            nPairs), np, sep = ""))
    }
    NCOL <- length(namesPairDIR)
    if (NCOL == 0) {
        outputEMMCMC.EM.W <- data.frame(DMCs = DMCs,
            pvalues = pval, qvalues = qval)
        names(outputEMMCMC.EM.W) <- c("DMCs", "pvalues",
            "qvalues")
    } else {
        DMCsPair <- as.matrix(vapply(totres3, function(elt) elt$dmcs,
            numeric(NCOL)))
        methDirPair <- as.matrix(vapply(totres3, function(elt) elt$methDir,
            character(NCOL)))
        if(NCOL>1){
            DMCsPair <- t(DMCsPair)
            methDirPair <- t(methDirPair)
        }
        DMCsPair[DMCs == 0, ] <- 0

        outputEMMCMC.EM.W <- data.frame(DMCs = DMCs,
            pvalues = pval, qvalues = qval, DMCsPair = DMCsPair,
            methDirPair = methDirPair)
        names(outputEMMCMC.EM.W) <- c("DMCs", "pvalues",
            "qvalues", namesPairDMC, namesPairDIR)
    }

    metadata(object)$DMCHMM <- outputEMMCMC.EM.W
    metadata(object)$formula <- formula

    return(object)
}

#' @rdname findDMCs-method
#' @aliases findDMCs-method findDMCs
setMethod("findDMCs", signature = c(object = "BSDMCs",
    formula = "ANY", FDRthreshold = "ANY", Methylthreshold = "ANY",
    mc.cores = "ANY", weightfunction =  "ANY"), .findDMCs)

