testCPsuggestion  <- function(rTarget=seq(80,1200,by=100),maxLag=20,nClutter=60,maxRange=1300){
#
# test a CP mode suggested for E3D by Ingemar
#
#
# IV 2023
#

    
    TXtimes <- c(0, 1, 4, 6, 7, 17, 22, 25, 29, 31, 32, 37, 41, 54, 57, 65, 72, 84, 86, 88, 98, 105, 111, 116, 119, 120, 125, 148, 158, 161, 169, 170, 188, 195, 212, 227, 241, 243, 247, 250, 261, 263, 267, 288, 295, 300, 310, 318, 319, 322, 328, 336, 338, 351, 358, 362, 363, 364 ,367 ,373 ,375 ,380 ,382 ,383)

    IPP  <-  c(diff(TXtimes),5) # not sure what the last IPP should be...

    IPP <- IPP * 3125/5 # 5 us bit length, IPPs are multiples of 3125 us
    
    # random coded experiment with the given IPPs
#    rexp <- randomCodedExp(64,128,2,1,IPP)
    AC128 <- searchStrongBinaryAlterCodes(7)
    rexp$code <- AC128[[1]]#[1:64]

    # range ambiguities for targets at different ranges...
    rAmbs <- list()
    for(k in seq(length(rTarget))){
        print(k)
        rAmbs[[k]] <- rangeAmbigMF(rexp,rTarget[k],maxRange,nClutter,maxLag,1)
    }

    return(rAmbs)
}


