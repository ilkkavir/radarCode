rangeAmbigMF <- function(targetRange,lags,radarExp,maxRange,nClutter,bitLen){
    # calculate range ambiguity functions for targets at a given range
    # for the given radar experiment in matched filter decoding. The ambiguities depend on range of the target,
    # because the IPPs are irregular.
    # IV 2023

    # create the continuos TX envelope
    TXenv <- makeTXenv(radarExp,bitLen)

    # The echo signal is a shifted copy of the envelope, but we do not have echoes during transmissions and for nClutter samples after each transmission. Zero-pad to make sure that we will not cut off anything from the end
    RXsig <- c(TXenv,numeric(length(TXenv)))
    RXsig[(targetRange+1):(length(TXenv)+targetRange)] <- TXenv

    IPP  <- rep(radarExp$IPP,length.out=length(radarExp$code))
    icur <- 1
    for(k in seq(length(radarExp$code))){
        RXsig[icur:(icur+length(radarExp$code[[k]])-1+nClutter)] <- 0
        icur <- icur + IPP[k]
    }
    
    # range ambiguities for each lag at each range from nClutter to maxRange
    nLags <- length(lags)
    rAmbMF <- matrix(0,(maxRange-nClutter+1),nLags)
    nTXenv <- length(TXenv)
    nRXsig <- length(RXsig)
    ii <- 1
    for(l in lags){
        rAmb <- TXenv[1:(nTXenv-l)] * Conj(TXenv[(l+1):nTXenv])
        lProd <- RXsig[1:(nRXsig-l)] * Conj(RXsig[(l+1):nRXsig])
        for(r in seq(nClutter,maxRange)){
            rAmbMF[r-nClutter+1,ii] <- sum( rAmb * Conj(lProd[(r+1):(r+nTXenv-l)]) )
        }
        if(any(abs(rAmbMF[,ii])>0)){
            rAmbMF[,ii] <- rAmbMF[,ii]/max(abs(rAmbMF[,ii]))
        }
        ii <- ii+1
    }

    return(rAmbMF)
    
    }
