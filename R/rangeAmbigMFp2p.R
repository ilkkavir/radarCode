rangeAmbigMFp2p <- function(targetRange,lags,radarExp,maxRange,nClutter,bitLen,edgeCoeff=1){
    ## calculate range ambiguity functions for targets at a given range
    ## for the given radar experiment when the data are first decoded by means of matched
    ## filtering in amplitude domain, and pulse-to-pulse lags are then calculated.
    ## The ambiguities depend on range of the target because the IPPs are irregular.
    ##
    ## IV 2023, 2025

    nCode <- length(radarExp$code)
    
    # create the continuos TX envelope
    TXenv <- rep(makeTXenv(radarExp,bitLen),2)

    # The echo signal is a shifted copy of the envelope, but we do not have echoes during transmissions and for nClutter samples after each transmission. Zero-pad to make sure that we will not cut off anything from the end
    RXsig <- c(TXenv,numeric(length(TXenv)))*0
    RXsig[(targetRange+1):(length(TXenv)+targetRange)] <- TXenv
    if(edgeCoeff!=1){
        phaseShifts <- abs(diff(RXsig))!=0
        nShifts <- c(phaseShifts,FALSE) + c(FALSE,phaseShifts)
        phaseShifts <- c(phaseShifts,FALSE) | c(FALSE,phaseShifts)
        RXsig[phaseShifts] <- RXsig[phaseShifts]*(1-nShifts[phaseShifts]*(1-edgeCoeff))
    }
    
    IPP  <- rep(radarExp$IPP,length.out=length(radarExp$code))*bitLen
    icur <- 1
    for(k in seq(nCode)){
        RXsig[icur:(icur+length(radarExp$code[[k]])*bitLen-1+nClutter)] <- 0
        icur <- icur + IPP[k]
    }

    ## matched filter "decoding" that produces point spread functions of individual pulses
    pointspread <- 0*TXenv
    pstarts <- cumsum(IPP)
    icur  <- 1
    for(ip in seq(nCode)){
        for(r in seq(nClutter+1,maxRange)){
            pointspread[ icur + r ] <- sum( RXsig[(icur+r):(pstarts[ip])] * Conj(TXenv[icur:(pstarts[ip]-r)]) )
        }
        icur <- icur + IPP[ip]
    }
    
    # range ambiguities for each lag at each range from nClutter to maxRange
    nLags <- length(lags)
    rAmbMF <- matrix(0,(maxRange-nClutter),nLags)
    nn <- length(pointspread)/2
    ii <- 1
    icur <- 1
    for(l in lags){
        lProd <- pointspread[1:nn] * Conj(pointspread[(l+1):(nn+l)])
        icur <- 1
        for(ip in seq(nCode)){
            rAmbMF[,ii] <- rAmbMF[,ii] + lProd[((nClutter+1):maxRange)+icur]
            icur <- icur + IPP[ip]
        }
        if(any(abs(rAmbMF[,ii])>0)){
            rAmbMF[,ii] <- rAmbMF[,ii]/max(abs(rAmbMF[,ii]))
        }
        ii <- ii+1
    }

    return(list(ramb=rAmbMF,pointspread=pointspread))
    
    }
