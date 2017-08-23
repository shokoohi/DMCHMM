.findDMCs <- function(object, formula, FDRthreshold, Methylthreshold, mc.cores){

    if (missing(object) | class(object) != "BSDMCs")
        stop("A BSDMCs object must be provided.")

    if (missing(formula)) {
        formula = as.formula(paste("Methylation",
            paste(c(names(colData(object))), collapse = " + "), sep = " ~ "))
    } else if (!class(formula) == "formula") {
        stop("Either provide a formula or leave it blank.")
    } else {
        if (!all(all.vars(formula) %in% names(colData(object)))) {
            stop("The variables that are provided in the formula do no match
            any variable in", paste(names(colData(object)), collpase=", "), ".")
        } else {
            formula = as.formula(paste("Methylation ~ ",
                paste(all.vars(formula)[all.vars(formula) %in%
    names(colData(object))], collapse = " + ")))
            # formula = as.formula(paste('Methylation ~ ',
            # paste(c(names(colData(object)),
            # all.vars(formula)), collapse = ' + ')))
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
        mc.cores = max(1, detectCores() - 1)
    } else if (!is.numeric(mc.cores) | mc.cores <= 0 |
        mc.cores > detectCores()) {
        stop("An integer value greater than 0 and less than or equal to
    detectCores() must be provided for mc.cores")
    }
    mc.cores = as.integer(mc.cores)

    strand(object) <- "*"
    object <- sort(object)

    object.df <- as.data.frame(colData(object))
    ind1 <- sapply(object.df, is.character)
    ind2 <- sapply(object.df, is.factor)
    object.df[ind1] <- lapply(object.df[ind1], as.factor)
    object.df <- object.df[all.vars(formula)[all.vars(formula) %in%
        names(colData(object))]]
    object.df.factor <- object.df[, ind2]
    object.df.factor = as.data.frame(object.df.factor)
    names(object.df.factor) <- names(object.df)[ind2]

    if (ncol(object.df) < 1)
        stop("There is no covarite information in the data.")

    # nGroup = length(unique(colData(object))[[1]])
    # if(nGroup<=1) stop('Assign groups using
    # colData.')

    nPos = nrow(object)
    nSam = length(colnames(object))

    optbp <- MulticoreParam(workers = mc.cores, progressbar = TRUE)
    .mfun3 <- function(itr){
        dat = data.frame(Methylation = c(assays(object)$methLevels[itr, ]))
        dat = cbind(dat, object.df)
        dat$Methylation[dat$Methylation < 1e-10] <- 1e-10
        dat$Methylation[dat$Methylation > 0.9999999999] <- 0.9999999999
        dat$Methylation[is.na(dat$Methylation)] <- mean(dat$Methylation,
                                                        na.rm = TRUE)
        dat$Methylation <- log(dat$Methylation/(1 - dat$Methylation))
        dat$Methylation[is.na(dat$Methylation)] <- mean(dat$Methylation,
                                                        na.rm = TRUE)

        options(show.error.messages = FALSE)
        suppressWarnings(lmodel <- try(lm(formula = formula, data = dat),
            silent = TRUE))
        options(show.error.messages = TRUE)
        if ((class(lmodel) == "try-error")) {
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
                 baseL <- lmodel$coefficients[c(grep(names(object.df.factor)[1],
                    names(lmodel$coefficients), value = TRUE))][1]
                }

                ITR1 = 0
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

                        ITR1 = 1

                        options(show.error.messages = FALSE)
                        suppressWarnings(lmodelc <- try(glht(lmodel,
                                                linfct = bb, silent = TRUE)))
                        options(show.error.messages = TRUE)
                        if ((class(lmodelc) == "try-error")) {
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

    # itr = 1
    # cl <- makeCluster(mc.cores)
    # registerDoSNOW(cl)
    # pb <- txtProgressBar(max = nPos, style = 3)
    # progress <- function(n) setTxtProgressBar(pb, n)
    # opts <- list(progress = progress)
    # totres3 <- foreach(itr = 1:nPos, .options.snow = opts,
    #     .packages = c("stats", "multcomp", "SummarizedExperiment")) %dopar%
    #     {
    #         dat = data.frame(Methylation = c(assays(object)$methLevels[itr,
    #             ]))
    #         dat = cbind(dat, object.df)
    #         dat$Methylation[dat$Methylation < 1e-10] <- 1e-10
    #         dat$Methylation[dat$Methylation > 0.9999999999] <- 0.9999999999
    #         dat$Methylation[is.na(dat$Methylation)] <- mean(dat$Methylation,
    #             na.rm = TRUE)
    #         dat$Methylation <- log(dat$Methylation/(1 -
    #             dat$Methylation))
    #         dat$Methylation[is.na(dat$Methylation)] <- mean(dat$Methylation,
    #             na.rm = TRUE)
    #
    #         options(show.error.messages = FALSE)
    #         suppressWarnings(lmodel <- try(lm(formula = formula,
    #             data = dat), silent = TRUE))
    #         options(show.error.messages = TRUE)
    #         if ((class(lmodel) == "try-error")) {
    #             p.val <- NA
    #             dmcs <- NA
    #             methDir <- NA
    #             return(list(p.val = p.val, dmcs = dmcs, methDir = methDir))
    #         } else {
    #             fst <- summary(lmodel)$fstatistic
    #             p.val <- pf(fst[1], fst[2], fst[3], lower.tail = FALSE)
    #             attributes(p.val) <- NULL
    #
    #             if (dim(object.df.factor)[2] == 0) {
    #               return(list(p.val = p.val))
    #             } else {
    #               dmcs <- c()
    #               methDir <- c()
    #               if ("(Intercept)" %in% names(lmodel$coefficients)) {
    #                 baseL <- lmodel$coefficients[1]
    #               } else {
    #            baseL <- lmodel$coefficients[c(grep(names(object.df.factor)[1],
    #                   names(lmodel$coefficients), value = TRUE))][1]
    #               }
    #
    #               ITR1 = 0
    #               for (g in names(object.df.factor)) {
    #                 nGroup = length(levels(object.df.factor[, g]))
    #                 if (nGroup > 1) {
    #                   if ("(Intercept)" %in% names(lmodel$coefficients) &
    #                     ITR1 == 0) {
    #                     lc <- lmodel$coefficients[c("(Intercept)",
    #                       grep(g, names(lmodel$coefficients), value = TRUE))]
    #                     lc[-1] = lc[-1] + baseL
    #                   } else if (!("(Intercept)" %in%
    #                     names(lmodel$coefficients)) &
    #                     ITR1 == 0) {
    #                     lc <- lmodel$coefficients[c(grep(g,
    #                       names(lmodel$coefficients),
    #                       value = TRUE))]
    #                   } else {
    #                     lc <- lmodel$coefficients[c(grep(g,
    #                       names(lmodel$coefficients),
    #                       value = TRUE))]
    #                     lc <- c(baseL, baseL + lc)
    #                   }
    #                   d <- (outer(1/(1 + exp(-lc)),
    #                     1/(1 + exp(-lc)), "-"))
    #                   d = d[upper.tri(d)]
    #                   temp1 = d
    #                   temp1[d > Methylthreshold] = "hyper"
    #                   temp1[d < Methylthreshold] = "hypo"
    #                   temp1[d >= (-Methylthreshold) &
    #                     d <= Methylthreshold] = "equal"
    #                   methDir <- c(methDir, temp1)
    #                   bb <- NULL
    #                   bb$g <- "Tukey"
    #                   names(bb) <- g
    #                   if (!("(Intercept)" %in% names(lmodel$coefficients)) &
    #                     ITR1 == 0) {
    #                     attr(bb, "interaction_average") <- TRUE
    #                   } else {
    #                     attr(bb, "interaction_average") <- FALSE
    #                   }
    #                   attr(bb, "covariate_average") <- FALSE
    #                   attr(bb, "class") <- "mcp"
    #
    #                   ITR1 = 1
    #
    #                   options(show.error.messages = FALSE)
    #                   suppressWarnings(lmodelc <- try(glht(lmodel,
    #                     linfct = bb, silent = TRUE)))
    #                   options(show.error.messages = TRUE)
    #                   if ((class(lmodelc) == "try-error")) {
    #                     dmcs <- c(dmcs, rep(NA, nGroup))
    #                   } else {
    #                     hhh <- (summary(lmodelc)$test$pvalues)
    #                     attributes(hhh) <- NULL
    #                     dmcs <- c(dmcs, hhh)
    #                   }
    #                 }
    #               }
    #               dmcs <- (dmcs < 0.05) * 1
    #               return(list(p.val = p.val, dmcs = dmcs,
    #                 methDir = methDir))
    #             }
    #         }
    #     }
    # close(pb)
    # stopCluster(cl)

    pval <- qval <- DMCs <- namesPair <- c()
    pval <- vapply(totres3, function(elt) elt$p.val, numeric(1))

    fdr.fun <- function(x) {
        (fdrtool(x, statistic = "pvalue", plot = FALSE,
            verbose = FALSE))$qval
    }
    if (nPos < 1000) {
        qval <- fdr.fun(pval)
    } else {
        new.pval <- split(pval, ceiling(seq_along(pval)/500))
        qval <- as.vector(unlist(lapply(new.pval, fdr.fun)))
    }

    DMCs = (qval < FDRthreshold) * 1

    nfactor <- dim(object.df.factor)[2]
    NCOL <- 0
    namesPairDMC <- namesPairDIR <- NULL
    if (nfactor > 0) {
        for (i in 1:nfactor) {
            nGroup = length(levels(object.df.factor[,
                i]))
            if (nGroup > 1) {
                NCOL = NCOL + ncol(combn(1:nGroup,
                  2))
                np <- paste(combn(levels(object.df.factor[,
                  i]), 2)[1, ], combn(levels(object.df.factor[,
                  i]), 2)[2, ], sep = "vs")
                namesPairDMC <- c(namesPairDMC, paste("DMCs",
                  names(object.df.factor[i]), np, sep = ""))
                namesPairDIR <- c(namesPairDIR, paste("methDir",
                  names(object.df.factor[i]), np, sep = ""))
            }
        }
    }

    if (NCOL == 0) {
        outputEMMCMC.EM.W <- data.frame(DMCs = DMCs,
            pvalues = pval, qvalues = qval)
        names(outputEMMCMC.EM.W) <- c("DMCs", "pvalues",
            "qvalues")
    } else {
        DMCsPair = matrix(NA, nrow = nPos, ncol = NCOL)
        methDirPair = matrix(NA, nrow = nPos, ncol = NCOL)

        DMCsPair <- vapply(totres3, function(elt) elt$dmcs,
            numeric(NCOL))
        methDirPair <- vapply(totres3, function(elt) elt$methDir,
            character(NCOL))
        # for (i in 1:nPos) { DMCsPair[i, ] =
        # totres3[[i]]$dmcs methDirPair[i, ] =
        # totres3[[i]]$methDir }
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
    mc.cores = "ANY"), .findDMCs)

