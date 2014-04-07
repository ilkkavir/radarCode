optimizeR <- function(nCode,nBit,nPhase,nIter=1000,maxLag=1,verb=FALSE){
# random search for optimal code groups using the FFT method

    code1 <- randomCodeSequence(nCode=nCode,nBits=nBit,nPhases=nPhase)
    temp <- evalR(code=code1,lags=seq(1,maxLag),fracs=c(0),bitLen=1,resolution=1)
    eval1 <- 0
    for (j in seq(1,maxLag)) eval1 <- eval1 + temp[[j]]
    eval1 <- eval1/maxLag

    for (k in seq(1,nIter)){
        nf<-ceiling(5*runif(1))
        code2 <- code1
        for(j in seq(1,nf)){
            nc <- ceiling(nCode*runif(1))
            nb <- ceiling(nBit*runif(1))
            code2[[nc]][nb] <- exp(pi*2i*ceiling(runif(1)*nPhase)/nPhase)*code1[[nc]][nb]
        }
        temp <- evalR(code=code2,lags=seq(1,maxLag),fracs=c(0),bitLen=1,resolution=1)
        eval2 <- 0
        for (j in seq(1,maxLag)) eval2 <- eval2 + temp[[j]]
        eval2 <- eval2/maxLag
        if (eval2 > eval1){
            eval1 <- eval2
            code1 <- code2
        }

        if (k%%10 == 0){
            code2 <- randomCodeSequence(nCode=nCode,nBits=nBit,nPhases=nPhase)
            temp <- evalR(code=code2,lags=seq(1,maxLag),fracs=c(0),bitLen=1,resolution=1)
            eval2 <- 0
            for (j in seq(1,maxLag)) eval2 <- eval2 + temp[[j]]
            eval2 <- eval2/maxLag
            if (eval2 > eval1){
                eval1 <- eval2
                code1 <- code2
            }
        }
        if (verb){
            if(k%%100 == 0) cat(sprintf("\r %5.2f %7.5f",k/nIter*100,eval1))
                                        #       print(c(k,eval1))
        }
        if (eval1==1) break
    }

    result <- list()
    result[['code']] <- code1
    result[['eval']] <- eval1

    return(result)

} #optimizeR
