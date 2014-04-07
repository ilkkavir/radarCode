evalR <- function(code=list(),lags=c(1),fracs=c(0),bitLen=c(1),resolution=1)
    {
        # calculate the R-value for a given set of codes at given lags

        R <- list()
        nCodes <- length(code)
        nLags <- length(lags)
        nFracs <- length(fracs)
        nBauds <- length(code[[1]])

        for (k in seq(1,length(code)))
            {
                if (length(code[[k]])>nBauds) nBauds <- length(code[[k]])
            }
        fftLength <- nextn(10*nBauds*bitLen)


        code2 <- lengthenBits(code,bitLen)

        convVec <- rep(0,fftLength)

        if (bitLen==1){
                                        # boxcar impulse response
            convVec[1:resolution] <- 1
        }else{
                                        # triangle impulse response
            convVec[1:resolution] <- seq(1,resolution)/resolution
            convVec[(resolution+1):(2*resolution)] <- (1.0 - seq(1,resolution)/resolution)
        }
        convVec <- convVec/resolution


        fftConv <- rep(0,fftLength)
        fftConv <- (abs(fft(convVec)))^2

        if(nLags==0){
            fftSum <- rep(0,fftLength)
            nTot <- 0
            for (j in seq(1,nCodes))
                {
                    fftVec <- rep(0.0,fftLength)
                    fftVec[1:(nBauds*bitLen)] <- code2[[j]][1:(nBauds*bitLen)]
                    nTot  <- nTot + sum(abs(fftVec))#sum(fftVec!=0)
                    fftSum <- fftSum + ((Mod(fft(fftVec)))^2)
                }
            fftSum[fftConv==0] <- 1.0
            R[[1]] <- 1.0 / ( min(bitLen) * resolution * nTot * sum(fftConv/fftSum) / fftLength )
            if (is.na(R[[1]])|(R[[1]]==Inf)) R[[k]] <- 0

            return(R)
        }



        for (k in seq(1,nLags))
            {
                fftSum <- rep(0,fftLength)
                nTot <- 0
                for (m in seq(1,nFracs))
                    {
                        for (j in seq(1,nCodes))
                            {
                                fftVec <- rep(0.0,fftLength)
                                fftVec[1:(nBauds*bitLen-lags[k]-fracs[m])] <- code2[[j]][1:(nBauds*bitLen-lags[k]-fracs[m])]*Conj(code2[[j]][(lags[k]+fracs[m]+1):(nBauds*bitLen)])
                                nTot  <- nTot + sum(abs(fftVec))#sum(fftVec!=0)
                                fftSum <- fftSum + ((Mod(fft(fftVec)))^2)
                            }
                    }
                fftSum[fftConv==0] <- 1.0
                R[[k]] <- 1.0 / ( min(bitLen) * resolution * nTot * sum(fftConv/fftSum) / fftLength )
                if (is.na(R[[k]])|(R[[k]]==Inf)) R[[k]] <- 0
            }

        return(R)

    } #evalR

