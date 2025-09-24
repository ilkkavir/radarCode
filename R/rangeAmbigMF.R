rangeAmbigMF <- function(targetRange,lags,radarExp,maxRange,nClutter,bitLen,monostatic=TRUE){
    # calculate range ambiguity functions for targets at a given range
    # for the given radar experiment in matched filter decoding. The ambiguities depend on range of the target,
    # because the IPPs are irregular.
    # IV 2023, 2025

    ## seems to work now, although very slow..

    
    # create the continuos TX envelope
    TXenv <- rep(makeTXenv(radarExp,bitLen),2)
    nTXenv <- length(TXenv)
    ## filtering...
    if(bitLen>1){
        
        for(k in seq(nTXenv,bitLen,by=-1)){
            TXenv[k] <- sum(TXenv[(k-bitLen+1):k])/bitLen
        }
        for(k in seq((bitLen-1),1,by=-1)){
            TXenv[k] <- sum(TXenv[1:k])/k
        }
        
    }
    
    ## The echo signal is a shifted copy of the envelope, but we do not have echoes during transmissions and for nClutter samples after each transmission. Zero-pad to make sure that we will not cut off anything from the end
    RXsig <- rep(0,nTXenv*2)#c(TXenv,numeric(length(TXenv)))
    RXsig[(targetRange+1):(nTXenv+targetRange)] <- TXenv

    IPP  <- rep(radarExp$IPP*bitLen,length.out=length(radarExp$code))
    icur <- 1
    if(monostatic){
        for(k in seq(length(radarExp$code))){
            RXsig[icur:(icur+length(radarExp$code[[k]])*bitLen-1+nClutter)] <- 0
            icur <- icur + IPP[k]
        }
    }

    ## the decoding filter contains only the samples at full integer multiples of bit length
    decoEnv <- TXenv*0
    decoii <- seq(bitLen,nTXenv,by=bitLen)
    decoEnv[decoii] <- TXenv[decoii]
    
    ## range ambiguities for each lag at each range from nClutter to maxRange
    nLags <- length(lags)
    rAmbMF <- matrix(0,(maxRange-nClutter+1+1),nLags)
    nRXsig <- length(RXsig)
    ii <- 1
    for(l in lags){
        rAmb <- ( decoEnv[1:(nTXenv-l)] * Conj(decoEnv[(l+1):nTXenv]) )[1:(nTXenv/2)]
        lProd <- ( RXsig[1:(nRXsig-l)] * Conj(RXsig[(l+1):nRXsig]) )
        for(r in seq(nClutter,maxRange)){
            rAmbMF[r-nClutter+1,ii] <- sum( rAmb * Conj(lProd[(r+1):(r+nTXenv/2)]) )
        }
        if(any(abs(rAmbMF[,ii])>0)){
            rAmbMF[,ii] <- rAmbMF[,ii]/max(abs(rAmbMF[,ii]))
        }
        ii <- ii+1
    }

    return(rAmbMF)
    
    }
