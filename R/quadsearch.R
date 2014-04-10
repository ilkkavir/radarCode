

#
#
# Quadriphase code search from Juha Vierinen (2014)
#
#
#
#
#
#

# The variance of a single lag
#
lagVariance <- function(codeGroup=randomCodeGroup(codeLen=13,ncodes=4),
                        lag=1,fftlen=length(codeGroup[[1]])*10)
{
  cl <- length(codeGroup[[1]])
  nc <- length(codeGroup)
  fsum <- rep(0.0,fftlen)
  pwr <- 0

  for(i in 1:nc)
  {
    a <- codeGroup[[i]][1:(cl-lag)]
    b <- Conj(codeGroup[[i]][(lag+1):cl])

    pwr <- pwr + sum(abs(a*b)**2)
    zpcode <- c(a*b,rep(0,fftlen-cl+lag))

    fsum <- fsum + abs(fft(zpcode))**2

#    cat(sprintf("%d %1.2f %1.2f\n",i,pwr,pwr*sum(1/fsum)/fftlen))
  }
  lagvar <- pwr*sum(1/fsum)/fftlen
#  nvar2 <- (nc*(cl-lag))*sum(1/fsum)/fftlen
  npwr <- pwr/(nc*(cl-lag))
  res <- list(npwr=npwr,pwr=pwr,lagvar=lagvar)
}

calcAll <- function(codeGroup=randomCodeGroup(codeLen=13,ncodes=4))
{
  nc <- length(codeGroup)
  cl <- length(codeGroup[[1]])

  pwrs <- rep(0,cl-1)
  npwrs <- rep(0,cl-1)
  lagvars <- rep(0,cl-1)

  for(l in 1:(cl-1))
  {
    res <- lagVariance(codeGroup=codeGroup,lag=l)
    pwrs[l] <- res$pwr
    npwrs[l] <- res$npwr
    lagvars[l] <- res$lagvar
  }
  return(list(pwrs=pwrs,npwrs=npwrs,lagvars=lagvars))
}

phmat2cl <- function(phmat)
{
  cl <- list()
  for(i in 1:dim(phmat)[1])
  {
    cl[[i]] <- phmat[i,]
  }
  return(cl)
}

quadphasesearch <- function(ncodes=2,codelen=20,nphases=4, niter=2000,verbose=FALSE)
{
  phase <- exp(1i*2*pi/nphases)
  phmat <- array(dim=c(ncodes,codelen))
  for(i in 1:ncodes)
  {
    phmat[i,1:codelen] <- phase**round(2*rnorm(codelen)*nphases)
  }
  res <- calcAll(phmat2cl(phmat))
  bestvar <- mean(res$lagvars)
  bestcode <- phmat

  for(i in 1:niter)
  {
    nchanges <- floor(abs(rnorm(1)*5))+1
    tryph <- phmat
    for(i in 1:nchanges)
    {
      codeidx <- floor(runif(1)*ncodes)+1
      baudidx <- floor(runif(1)*codelen)+1
      tryph[codeidx,baudidx] <- tryph[codeidx,baudidx]**(round(rnorm(1)*nphases))
    }
    tryvar <- mean(calcAll(phmat2cl(tryph))$lagvars)
    if(tryvar < bestvar){
        if(verbose) cat("new best ",tryvar,"\n")
        bestvar <- tryvar
        bestcode <- tryph
    }
  }
  res <- calcAll(phmat2cl(bestcode))
  res$phmat <- bestcode
  res$var <- bestvar
  return(res)
}

