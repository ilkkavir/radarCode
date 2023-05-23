randomCodedExp <- function(nCode,nBit,nPhase,nAmp,IPP){

    experiment <- list()
    experiment$IPP <- IPP
    experiment$code  <- randomCodeSequence(nCode=nCode,nBits=nBit,nPhases=nPhase,nAmplitudes=nAmp)

    return(experiment)

    }
