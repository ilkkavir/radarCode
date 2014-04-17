codeSearchV <- function(nCode,nBit,nPhase=2,nIter=1000,limit=Inf,verb=FALSE,nCores=1,nRun=Inf){
    # code search with the FFT-method
    #
    # Similar with codeSearchR but uses the routines in file quadSearch.R from Juha
    # Vierinen. This version minimises V whereas codeSearchR maximizes R, but V=1/R...
    #
    codes <- list()
    sfile <- paste('codesV_',as.character(nCode),'_',as.character(nBit),'_',as.character(nPhase),'.Rdata',sep='')
    i <- 1
    nn <- 0
    repeat{
        if(nCores==1){
            code <- quadphasesearch(ncodes=nCode,codelen=nBit,nphases=nPhase,niter=nIter,verbose=verb)
            if (code$var < limit){
                print(code$var)
                if (any(dir()==sfile)){
                    load(sfile)
                    i <- length(codes)+1
                }
                codes[[i]] <- code
                save(codes,file=sfile)
            }
        }else{
            pList <- list()
            for(k in seq(nCores)){
                pList[[k]] <- parallel(quadphasesearch(ncodes=nCode,codelen=nBit,nphases=nPhase,niter=nIter,verbose=FALSE),mc.set.seed=TRUE)
            }
            for(k in seq(nCores)){
                code <- collect( pList[[k]] )[[1]]
                if (code$var < limit){
                    print(code$var)
                    if (any(dir()==sfile)){
                        load(sfile)
                        i <- length(codes)+1
                    }
                    codes[[i]] <- code
                    save(codes,file=sfile)
                }
            }
        }
        nn <- nn + nCores
        if( nn >= nRun ) return(nn)
    }
}
