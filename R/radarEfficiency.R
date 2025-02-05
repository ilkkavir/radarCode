radarEfficiency <- function(experiment,rangeLimits=seq(1,sum(unique(experiment$IPP/experiment$baudLength))),lagLimits=seq(0,sum(unique(experiment$IPP/experiment$baudLength))),remote=FALSE,maxr=Inf){
    ## 
    ## Calculate absolute radar efficiencies of modulation for given range-gates and time lags
    ##  
    ## INPUT:
    ##   experiment   an "experiment list" of phase-codes and IPPs
    ##   rangeLimits  range gate limits (in sample units)
    ##   lagLimits    lag integration limits (in sample units)
    ##   
    ## IV 2025
    ## 


    ## a logical vector that true when TX is on, repeat twice to make the evaluation simpler
    expTmp <- experiment
    expTmp$IPP <- expTmp$IPP/experiment$baudLength
    TXbit <- makeTXenv(expTmp,1)!=0
    ns <- length(TXbit)
    nl <- length(lagLimits)-1
    nr <- length(rangeLimits)-1
    
    ## a matrix for the efficiencies (first without integration in range/lag)
    minLag <- min(lagLimits)
    maxLag <- max(lagLimits)-1
    minRange <- min(rangeLimits)
    maxRange <- max(rangeLimits)-1
    nl1 <- maxLag - minLag + 1
    nr1 <- maxRange - minRange + 1
    eff1 <- matrix(0,nrow=nr1,ncol=nl1)

    nrep <- ceiling((max(maxRange,maxLag)+ns)/length(TXbit)) + 1
    TXbit <- rep(TXbit,nrep)

    if(length(maxr)<nl){
        maxr <- c(maxr,rep(maxr[length(maxr)],nl-length(maxr)))
    }
    
    ## calculate efficiencies at each range and lag (the lowest range will be at the last row...)
    for(il in seq(minLag,maxLag)){
        icol <- il - minLag + 1
        Rbit <- TXbit[1:(nrep*ns-il)] & TXbit[(il+1):(nrep*ns)]
        OKbit <- !TXbit[1:(nrep*ns-il)] & !TXbit[(il+1):(nrep*ns)]

        ## all data samples are usable at remotes
        if(remote){
            OKbit[] <- TRUE
        }
        
        lowCount <- 0
        for(rl in seq(1,maxRange)){
            if(Rbit[rl]){
                lowCount <- 0
            }else{
                lowCount <- lowCount + 1
            }
        }
        for(rl in seq(maxRange+1,maxRange+ns)){
            if(Rbit[rl]){
                lowCount <- 0
            }else{
                lowCount <- lowCount + 1
            }
            if(remote) lowCount <- Inf
            if( OKbit[rl]  & (lowCount>=minRange)){
                eff1[,icol] <- eff1[,icol] + Rbit[(rl-maxRange):(rl-minRange)]
            }
        }
    }
    eff1 <- eff1/ns

    ## Integration in range and lag, plus flip the range axis
    nl <- length(lagLimits)-1
    nr <- length(rangeLimits)-1
    eff <- matrix(NA,ncol=nl,nrow=nr)
    for(il in seq(nl)){
        lagInds <- c(lagLimits[il],lagLimits[il+1]-1) - minLag + 1
        for(ir in seq(nr)){
            rangeInds <- c(rangeLimits[ir],rangeLimits[ir+1]-1) - minRange + 1
            eff[ir,il] <- mean(eff1[(nr1-rangeInds[2]+1):(nr1-rangeInds[1]+1),lagInds[1]:lagInds[2]])
        }
    }

    out <- list()
    out$efficiency <- eff
    out$rangeLimits <- rangeLimits
    out$lagLimits <- lagLimits
    out$range <- (rangeLimits[1:nr]+rangeLimits[2:(nr+1)]-1)/2
    out$lag <- (lagLimits[1:nl]+lagLimits[2:(nl+1)]-1)/2

    ## cut off unsolved columns like in LPI.gdf
    for(il in seq(nl)){
        out$efficiency[out$range>maxr[il],il] <- NA
    }

    rmcols <- maxr<=0
    out$efficiency <- out$efficiency[,!rmcols]
    out$lag <- out$lag[!rmcols]
    
    # reverse the colums of eff1
    out$eff1 <- apply(eff1,FUN=rev,MARGIN=2)
    out$lag1 <- seq(minLag,maxLag)
    out$range1 <- seq(minRange,maxRange)


    return(out)
    
}
