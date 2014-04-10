codeSearchR <- function(nCode,nBit,nPhase=2,nIter=1000,maxLag=1,limit=0,verb=FALSE,nCores=1){
# code search with the FFT-method
  codes <- list()
  sfile <- paste('codes_',as.character(nCode),'_',as.character(nBit),'_',as.character(nPhase),'.Rdata',sep='')
  i <- 1
  while(TRUE){
      if(nCores==1){
          code <- optimizeR(nCode=nCode,nBit=nBit,nPhase=nPhase,nIter=nIter,maxLag=maxLag,verb=verb)
          if (code$eval > limit){
              print(code$eval)
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
              pList[[k]] <- parallel(optimizeR(nCode=nCode,nBit=nBit,nPhase=nPhase,nIter=nIter,maxLag=maxLag,verb=FALSE),mc.set.seed=TRUE)
          }
          for(k in seq(nCores)){
              code <- collect( pList[[k]] )[[1]]
              if (code$eval > limit){
                  print(code$eval)
                  if (any(dir()==sfile)){
                      load(sfile)
                      i <- length(codes)+1
                  }
                  codes[[i]] <- code
                  save(codes,file=sfile)
              }
          }
      }
  }
} #runRsearch
