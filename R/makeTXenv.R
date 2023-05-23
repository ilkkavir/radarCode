makeTXenv <- function(radarExp,bitLen){

    IPP <- rep(radarExp$IPP,length.out=length(radarExp$code))

    NS <- sum(unlist(IPP))

    TXenv <- numeric(NS)


    icur <- 1
    for(k in seq(length(radarExp$code))){
        nbit  <- length(radarExp$code[[k]])
        TXenv[icur:(icur+nbit-1)] <- radarExp$code[[k]]
        icur <- icur + IPP[[k]]
    }

    TXenv <- rep(TXenv,each=bitLen)

    return(TXenv)

}



