testPPATC  <- function(rTarget=seq(80,1200,by=100),maxLag=20,nClutter=60,maxRange=1300){
#
# test a PPATC version that fills the radar duty-cycle
#
#
# IV 2023
#


    IPP <- c(1,2)*1710/5
    
    # random coded experiment with the given IPPs
    rexp <- randomCodedExp(64,128,2,1,IPP)

    # range ambiguities for targets at different ranges...
    rAmbs <- list()
    for(k in seq(length(rTarget))){
        print(k)
        rAmbs[[k]] <- rangeAmbigMF(rexp,rTarget[k],maxRange,nClutter,maxLag,1)
    }

    return(rAmbs)
}


