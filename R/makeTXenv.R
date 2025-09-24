makeTXenv <- function(radarExp,bitLen){
    ##
    ## make TX envelope for a given experiment and bit length.
    ##
    ## Notice that radarExp$IPP is given as number of one bit long samples in this version
    ## (scales with bitLen so that duty cycle is independent of bitLen)
    ##
    ## IV 2025
    ##

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



