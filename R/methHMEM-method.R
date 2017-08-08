.EMHMM.M <- function(xv,Kv,Yv,nv, MaxEmiter, epsEM, useweight){
    if(useweight){
        nw <- log(nv+2)
    }else{
        nw <- rep(1,length(Yv))
    }

    nPos = length(Yv)

    pX.filter = array(0, c(nPos, Kv+2))
    dimnames(pX.filter) = list(index = 1:nPos, states = c(1:(Kv+2)))

    X.filter = array(0, c(nPos))
    dimnames(X.filter) = list(index = 1:nPos)

    pX.predict = array(0, c(nPos, Kv+2))
    dimnames(pX.predict) = list(index = 1:nPos, states = c(1:(Kv+2)))

    pX.smooth = array(0, c(nPos, Kv+2))
    dimnames(pX.smooth) = list(index = 1:nPos, states = c(1:(Kv+2)))

    X.smooth = array(0, c(nPos))
    dimnames(X.smooth) = list(index = 1:nPos)

    W01k.Estep = array(0, c(Kv+2))
    dimnames(W01k.Estep) = list(states = c(1:(Kv+2)))

    Wik.Estep = array(0, c(nPos, Kv+2))
    dimnames(Wik.Estep) = list(index = 1:nPos, states = c(1:(Kv+2)))

    Wijk.Estep = array(0, c(nPos-1, Kv+2, Kv+2))
    dimnames(Wijk.Estep) = list(index = 2:nPos, states = c(1:(Kv+2)),
                                states = c(1:(Kv+2)))

    pv<-exp(c(xv[1:Kv],0))
    pv<-pv/sum(pv)
    pv <-c(0,cumsum(pv))

    Pm<-matrix(0,nrow=Kv+2,ncol=Kv+2)
    Pm[,1:(Kv+1)] <- matrix(xv[c((Kv+1):(Kv+(Kv+2)*(Kv+1)))],Kv+2,Kv+1)
    Pm <- t(apply(exp(Pm),1,function(x){return(x/sum(x))}))

    #Pi0 <- exp(c(xv[((Kv+(Kv+2)*(Kv+1)+(1))):(Kv+(Kv+2)*(Kv+1)+(Kv+1))],0))
    #Pi0 <- Pi0 / sum(Pi0)

    # folowing lines will be activated if we are seeking an alternative
    # estimator for initial probabilities.
    Pi0 <- Pm
    for(i in 1:100)
        Pi0 <- Pi0 %*% Pm
    Pi0 <- Pi0[1,]

    diffpar = 1
    niter = 1

    while((diffpar > epsEM ) & (niter < MaxEmiter)){
        # E-step
        pX.predict[1,]<-t(t(Pm) %*% Pi0  )
        # p(X_{t}=k|y_1^{t-1}?)  prediction
        pmarg<-sum(dbinom(Yv[1],nv[1],pv)*pX.predict[1,])
        # p(y_t)=\sum_k p(y_t,X_{t-1}=k) marginal
        pX.filter[1,] <- c(dbinom(Yv[1],nv[1],pv)*pX.predict[1,])/pmarg
        #filtering
        #Xfilter[1]<-which.max(pX.filter[1,])

        for(i in 2:nPos){
            pX.predict[i,] <- t(t(Pm) %*% pX.filter[i-1,])
            # p(X_{t}=k|y_1^{t-1}?)  prediction
            pmarg <- sum(dbinom(Yv[i],nv[i],pv)*pX.predict[i,])
            pmarg <- ifelse(pmarg<1e-200,1e-200,pmarg)
            # p(y_t)=\sum_k p(y_t,X_{t-1}=k) marginal
            pX.filter[i,] <- c(dbinom(Yv[i],nv[i],pv)*pX.predict[i,])/pmarg
            #filtering
            #Xfilter[i]<-which.max(pX.filter[i,])
            # mlike<- mlike+nw[i]*log(pmarg)
        }
        #pX.smooth[is.na(pX.smooth)]=1/(Kv+2)
        pX.smooth[nPos,]<- pX.filter[nPos,]
        #Xsmooth[nPos]<-which.max(pX.smooth[nPos,])
        for(i in (nPos-1):1){
            #if(is.na(sum(pX.smooth[i,]))) stop(print(i))
            hh <- (pX.smooth[i+1,]/pX.predict[i+1,])
            hh[is.na(hh)]<-0
            pX.smooth[i,]<-pX.filter[i,] * (Pm %*% hh)
            #  Xsmooth[i]<-which.max(pX.smooth[i,])
        }
        #pen <- - sum(log( abs(pv[-c(1)]-pv[-c(Kv+2)])))
        #mlike<- mlike - lamb * pen

        W01k.Estep =  pX.smooth[1,]
        Wik.Estep = pX.smooth

        HWijk.Estep = Wijk.Estep

        for(i in 2:nPos){
            Wijk.Estep [i-1,,] =
            ((Pm * matrix((pX.filter[i-1,]), Kv+2, Kv+2)) *
    matrix(pX.smooth[i,], Kv+2, Kv+2, byrow = TRUE))/
            matrix(pX.predict[i,], Kv+2, Kv+2, byrow = TRUE)}
        Wijk.Estep[is.na(Wijk.Estep)] <- 0
        #new.Pi0 = W01k.Estep / sum(W01k.Estep)
        new.pv = ((nw*Yv) %*% Wik.Estep) / ((nw*nv) %*% Wik.Estep)

        new.Pm<-matrix(0, Kv+2, Kv+2)

        for(j in 1:(Kv+2)){
            for(k in 1:(Kv+2))
                new.Pm[j,k] = sum( Wijk.Estep[,j,k])/sum( Wijk.Estep[,j,])}

        # folowing lines will be activated if we are seeking an alternative
        # estimator for initial probabilities.
        tmppi <- new.Pm
        for(i in 1:100)
            tmppi = tmppi %*% new.Pm
        new.Pi0 = tmppi[1,]

        diffpar = 0
        diffpar = diffpar + sum((Pi0-new.Pi0)^2)
        diffpar = diffpar + sum((Pm - new.Pm)^2)
        diffpar = diffpar + sum((pv - new.pv)^2)
        #if(is.na(diffpar)) stop(paste(niter))

        niter = niter + 1

        Pi0 = new.Pi0
        Pm = new.Pm
        pv = new.pv
    }

    mlike <-0
    pX.predict[1,]<-t(t(Pm) %*% Pi0  )  # p(X_{t}=k|y_1^{t-1}?)  prediction
    pmarg<-sum(dbinom(Yv[1],nv[1],pv)*pX.predict[1,])
    # p(y_t)=\sum_k p(y_t,X_{t-1}=k) marginal
    pX.filter[1,] <- c(dbinom(Yv[1],nv[1],pv)*pX.predict[1,])/pmarg
    #filtering
    X.filter[1]<-which.max(pX.filter[1,])
    mlike<-mlike+nw[1]*log(pmarg)

    for(i in 2:nPos){
        pX.predict[i,]<-t(t(Pm) %*% pX.filter[i-1,]  )
        # p(X_{t}=k|y_1^{t-1}?)  prediction
        pmarg<-  sum(dbinom(Yv[i],nv[i],pv)*pX.predict[i,])
        # p(y_t)=\sum_k p(y_t,X_{t-1}=k) marginal
        pX.filter[i,]<-c(dbinom(Yv[i],nv[i],pv)*pX.predict[i,])/pmarg
        #filtering
        X.filter[i]<-which.max(pX.filter[i,])
        mlike<- mlike+nw[i]*log(pmarg)
    }
    pX.smooth[nPos,]<-pX.filter[nPos,]
    X.smooth[nPos]<-which.max(pX.smooth[nPos,])
    for(i in (nPos-1):1){
        hh=(pX.smooth[i+1,]/pX.predict[i+1,])
        hh[is.na(hh)]=0
        pX.smooth[i,]<-pX.filter[i,] *
            (Pm %*%hh)
        X.smooth[i]<-which.max(pX.smooth[i,])
    }
    X.estimate = rowSums(pX.smooth * matrix(pv, ncol = Kv+2, nrow = nPos,
                                            byrow = TRUE))
    X.estimate[X.estimate>1]=1
    X.estimate[X.estimate<0]=0

    #pen <- - sum(log( abs(pv[-c(1)]-pv[-c(Kv+2)])))
    mlike = 0.0
    for(i in 1:nPos){
        mlike = mlike + nw[i] * dbinom(Yv[i],nv[i],X.estimate[i],log = TRUE)
    }

    return(list(mlike=mlike, beta=pv, Pm=Pm, Pi0=Pi0, X.smooth=X.smooth,
                X.filter=X.filter, X.estimate=X.estimate))
}

.runEM.M <- function(yy, nn, MaxK, MaxEmiter, epsEM, useweight){
    K=1
    xstart<-rep(0,K+(K+2)*(K+1)+(K+1))
    xEM <- .EMHMM.M(xstart, K, yy, nn, MaxEmiter, epsEM, useweight = useweight)
    logMlik <- xEM$mlike
    BIC <--2*logMlik+log(length(yy))*(K+(K+2)*(K+1)+(K+1))
    BICmax = 10^20
    if(!is.na(BIC)){
        BICmax = BIC
        EstHMMSam1 = list(xEM=xEM, K=K, logMlik=logMlik, BIC=BIC)
    }else{
        BIC=10^20
    }
    while ((K<=MaxK) & (BICmax>=BIC) ) {
        K=K+1
        xstart<-rep(0,K+(K+2)*(K+1)+(K+1))
        xEM <-.EMHMM.M(xstart, K, yy, nn, MaxEmiter, epsEM, useweight)
        logMlik <- xEM$mlike
        BIC <--2*logMlik+log(length(yy))*(K+(K+2)*(K+1)+(K+1))
        if((BIC<BICmax) & (!is.na(BIC))){
            BICmax=BIC
            EstHMMSam1 = list(xEM=xEM, K=K, logMlik=logMlik, BIC=BIC)
        }
    }
    return(EstHMMSam1)
}

.methHMEM <- function(object, MaxK, MaxEmiter, epsEM, useweight,
    mc.cores){

    if(missing(object) | class(object)!="BSData")
        stop("A BSData object must be provided.")

    if(missing(MaxK)){
        MaxK = 10
    }else if(!is.numeric(MaxK) | MaxK<1){
        stop("An integer value greated than 0 must be provided for MaxK")
    }
    MaxK = as.integer(MaxK)

    if(missing(MaxEmiter)){
        MaxEmiter = 300
    }else if(!is.numeric(MaxEmiter) | MaxEmiter<10){
        stop("An integer value greated than 10 must be provided for MaxEmiter.")
    }
    MaxEmiter = as.integer(MaxEmiter)

    if(missing(epsEM)){
        epsEM = 1e-5
    }else if(!is.numeric(epsEM) | epsEM<=0 | epsEM>=0.1){
        stop("A numeric value between 0 and 0.1 must be provided for epsEM.")
    }

    if(missing(useweight)){
        useweight = TRUE
    }else if(!is.logical(useweight)){
        stop("A logical value must be provided for useweight.")
    }
    useweight = as.logical(useweight)

    if(missing(mc.cores)){
        mc.cores = max(1, detectCores()-1)
    }else if(!is.numeric(mc.cores) | mc.cores<=0 | mc.cores>detectCores()){
        stop("An integer value greater than 0 and less than or equal to
    detectCores() must be provided for mc.cores")
    }
    mc.cores = as.integer(mc.cores)

    strand(object) <- "*"
    object <- sort(object)

    nGroup = length(unique(colData(object))[[1]])
    if(nGroup<=1)
        stop("Assign groups using colData.")

    nPos = nrow(object)
    nSam = length(colnames(object))

    itr = 1
    cl <- makeCluster(mc.cores)
    registerDoSNOW(cl)
    pb <- txtProgressBar(max = nSam, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    totres <- foreach(itr=1:nSam, .options.snow = opts) %dopar%{
        EstHMMSam <- .runEM.M(c(assays(object)$methReads[,itr]),
    c(assays(object)$totalReads[,itr]), MaxK, MaxEmiter,epsEM, useweight)
    return(EstHMMSam)
    }
    close(pb)
    stopCluster(cl)

    Khat = c()
    Betahat = list()
    Phat = list()
    methStates.new = methLevels.new = matrix(0, nPos, nSam)

    for(i in 1:nSam){
        Khat[i] = totres[[i]]$K
        Betahat[[i]] = totres[[i]]$xEM$beta
        Phat[[i]] = totres[[i]]$xEM$Pm
        methLevels.new[,i] = c(totres[[i]]$xEM$X.estimate)
        methStates.new[,i] = as.integer(c(totres[[i]]$xEM$X.smooth-1))
    }

    colnames(methLevels.new) <- colnames(object)
    rownames(methLevels.new) <- NULL
    colnames(methStates.new) <- colnames(object)
    rownames(methStates.new) <- NULL

    predictedMeth <- cBSDMCs(colData = colData(object),
                            rowRanges = rowRanges(object),
                            methReads = methReads(object),
                            totalReads = totalReads(object),
                            methStates = apply(methStates.new, 2, as.integer),
                            methLevels = methLevels.new
    )
    metadata(predictedMeth)$K = Khat
    metadata(predictedMeth)$Beta = Betahat
    metadata(predictedMeth)$Pm = Phat

    return(predictedMeth)
    }

#' @rdname methHMEM-method
#' @aliases methHMEM-method methHMEM
setMethod("methHMEM",
    signature=c(object = "BSData", MaxK = "ANY", MaxEmiter = "ANY",
    epsEM = "ANY", useweight = "ANY", mc.cores = "ANY"),
    .methHMEM)

