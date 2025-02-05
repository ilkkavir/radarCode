#
# A collection to make strong and weak type 1 alternating codes
#
#
# I. Virtanen 2010
#

searchStrongBinaryAlterCodes <- function(m){
# return all sequences of 2^m baud strong binary alternating codes
  P <- makePseq(m)
  k <- find.k(P)
  codePseqs <- list()
  nCode <- 0

  for (n in seq(1,length(k))){
    Pt <- list()
    for (l in seq(2,length(P))){
      Pt[[l-1]] <- trans.p(P[[l]],k[n])
    }
    tmpPseq <- test.k(P,Pt)
    if (!is.null(tmpPseq)){
      nCode <- nCode + 1
      codePseqs[[nCode]] <- tmpPseq
    }
  }


  codeSeqs <- makeCodeSeqs(codePseqs)

  codes <- makeCodes(codeSeqs,2^m)

  return(codes)

} # searchStrongBinaryAlterCodes

searchWeakBinaryAlterCodes <- function(m){
# return all sequences of 2^m baud weak binary alternating codes
  P <- makePseq(m)
  k <- find.k(P)
  codePseqs <- list()
  nCode <- 0

  for (n in seq(1,length(k))){
    Pt <- list()
    for (l in seq(2,length(P))){
      Pt[[l-1]] <- trans.p(P[[l]],k[n])
    }
    tmpPseq <- test.k(P,Pt)
    if (!is.null(tmpPseq)){
      nCode <- nCode + 1
      codePseqs[[nCode]] <- tmpPseq
    }
  }

  codeSeqs <- makeCodeSeqs(codePseqs)

  # conversion to weak codes
  for(k in seq(length(codeSeqs))){
    codeSeqs[[k]] <- floor(codeSeqs[[k]]/2)
  }

  codes <- makeCodes(codeSeqs,2^m,strong=F)

  return(codes)

} # searchWeakBinaryAlterCodes


makePseq <- function(m){
# make the p-sequences needed in search of
# strong binary alternating codes

  Pseq <- list()
  aMax <- 2^m - 1
  for (a in seq(0,aMax)){
    if (isEven(a)){
      Pseq <- c(Pseq,list(make.p(a,m)))
    }
  }

  return(Pseq)

} #makePseq

find.k <- function(P){
# find the possible k values for the search of strong alternating codes

  np <- length(P)
  k <- c()
  for (m in seq(1,np)){
    if (length(P[[m]])==2){
      k <- c(k,binxor(P[[m]]))
    }
  }

  return(k)

} #find.k

binxor <- function(p){
# exclusive OR operation of the binary presentations of
# the two decimal values of p

  charp <- intToBin(p)
  nC <- nchar(charp[1])
  reschar <- ''
  for (k in seq(1,nC)){
    if (substr(charp[1],k,k)==substr(charp[2],k,k)){
      reschar <- paste(reschar,'0',sep='')
    }else{
      reschar <- paste(reschar,'1',sep='')
    }
  }
  return(binToDec(reschar))

} #binxor

trans.p <- function(p,k){
# exclusive or operation for a p-sequence and k
  pTrans <- c()
  len <- length(p)
  for (m in seq(1,len)){
    pTrans[m] <- binxor(c(p[m],k))
  }

  return(pTrans)

} #trans.p

test.k <- function(P.orig,Pt.orig){

  P <- P.orig
  Pt <- Pt.orig
  pseq <- list()
  pseq[[1]] <- P[[1]]
  P[[1]] <- c()
  trans <- TRUE
  pNum <- 1

  while(TRUE){
    if(trans){
      Pcur <- Pt
    }else{
      Pcur <- P
    }

    found <- FALSE
    for (m in seq(1,length(Pcur))){
      if(Pcur[[m]][1]==pseq[[pNum]][length(pseq[[pNum]])]){
        pNum <- pNum + 1
        pseq[[pNum]] <- Pcur[[m]]
        Pt[[m]] <- c()
        P[[m]] <- c()
        trans <- !trans
        found <- TRUE
      }
    }
    if (!found){
      pseq <- NULL
      break
    }
    if (length(P)==1){
     break
    }
  }

  return(pseq)

} #test.k

makeCodeSeqs <- function(codePseqs){

  codeSeqs <- list()
  for (k in seq(1,length(codePseqs))){
    codeSeqs[[k]] <- codePseqs[[k]][[1]]
    for(m in seq(2,length(codePseqs[[k]]))){
      codeSeqs[[k]] <- c(codeSeqs[[k]],codePseqs[[k]][[m]][2:length(codePseqs[[k]][[m]])])
    }
  }

  return(codeSeqs)

} #makeCodeSeqs

makeCodes <- function(codeSeqs,nBaud,strong=T){

  codes <- list()
  for (k in seq(1,length(codeSeqs))){
    codes[[k]] <- alterCode(codeSeqs[[k]],nBaud=nBaud,strong=strong,octal=FALSE)
  }

  return(codes)

} #makeCodes

isEven <- function(a){
# returns TRUE if the number of 1s in the
# binary representation of a is even

  binChar <- intToBin(a)
  nC <- nchar(binChar)
  chars <- c()
  for(k in seq(1,nC)){
    chars[k]<-substr(binChar,k,k)
  }
  binNum <- as.integer(chars)

  result <- sum(binNum)%%2 == 0

  return(result)

} #isEven

make.p <- function(a,m){
# make one p-sequence for search of strong
# alternating codes

  p <- c(a)
  k <- 1
  pMax <- 2^m - 1
  while(p[k] <= pMax){
    k <- k+1
    p[k] <- ofunc(p[k-1])
  }

  return(p)

} #make.p

ofunc <- function(a){
# o-function for search of strong alternating codes
  binChar <- intToBin(a)
  if (isEven(a)){
    binChar <- paste(binChar,'1',sep='')
  }else{
    binChar <- paste(binChar,'0',sep='')
  }

  return(binToDec(binChar))

} #ofunc

binToDec <- function(binChar){
  chars <- c()
  nC <- nchar(binChar)
  for(k in seq(1,nC)){
    chars[k] <- substr(binChar,k,k)
  }
  binNum <- as.integer(chars)

  decNum <- 0
  for(k in seq(1,nC)){
    decNum <- decNum + binNum[nC-k+1]*2^(k-1)
  }

  return(decNum)

} #binToDec

alterCode <- function(indices,nBaud=0,strong=TRUE,octal=TRUE){
# alternating code correspoding the given Walsh indices

  codeList <- list()

  if (nBaud==0) {nBaud = length(indices)}

  if (octal){
    decInds <- octToDec(indices)
  }else{
    decInds <- indices
  }

  nB <- length(indices)
  nCode <- 2*nB

  walsh <- matrix(nrow=nCode,ncol=nCode)

  walsh <- walshMatrix(nCode)

  codes <- matrix(nrow=nCode,ncol=nB)

  for (k in seq(1,nB)){
    codes[,k] <- walsh[,decInds[k]+1]
    }

  if (strong){
    groupLen <- nCode
  }else{
    groupLen <- nB
    }

  for (k in seq(1,groupLen)){
    codeList[[k]] <- codes[k,1:nBaud]
    }

  return(codeList)

  } #alterCode

walshMatrix <- function(n){
# n x n Walsh matrix
#
#
  walsh <- matrix(ncol=2*n,nrow=2*n)
  walsh[1,1] <- 1

  for (k in seq(1,log2(2*n))){
    walsh[(2^(k-1)+1):(2^k),1:(2^(k-1))] <- walsh[1:(2^(k-1)),1:(2^(k-1))]
    walsh[(2^(k-1)+1):(2^k),(2^(k-1)+1):(2^k)] <- -walsh[1:(2^(k-1)),1:(2^(k-1))]
    walsh[1:(2^(k-1)),(2^(k-1)+1):(2^k)] <- walsh[1:(2^(k-1)),1:(2^(k-1))]
    }

  return(walsh[1:n,1:n])

  } #walshMatrix

octToDec <- function(octal){
# conversion from octal to decimal values

  binary <- c()

  if (length(octal)>1){

    for (k in seq(1,length(octal))){
      binary[k] <- octToDec(octal[k])
      }

  }else{

    binary <- 0
    if (is.numeric(octal)){
      charNum <- as.character(octal)
    }else{
      charNum <- octal
      }
    numLen <- nchar(charNum)
    for (k in seq(0,(numLen-1))){
       binary <- binary + as.integer(substr(charNum,(numLen-k),(numLen-k)))*8^k
      }

    }

  return(binary)

  } # octToDec

ACtest <- function(cg,maxLag=0){
# test if the given code is a type 1 alternating code
  nb <- length(cg[[1]])
  nc <- length(cg)
  if (maxLag==0)  maxLag <- nb-1
  eval <- TRUE
  for(lag in seq(1,maxLag)){
    lagAmb <- matrix(ncol=(nb-lag),nrow=nc)
    for (k in seq(1,nc)){
      lagAmb[k,] <- cg[[k]][1:(nb-lag)]*Conj(cg[[k]][(lag+1):nb])
    }
    tmp <- lagAmb
    for(k in seq(1,(nb-lag))){
      for(j in seq(1,(nb-lag))){
        tmp[,j] <- lagAmb[,j]*lagAmb[,k]
      }
      testsum <- rep(0,(nb-lag))
      for(j in seq(1,(nb-lag))){
        testsum[j] <- sum(tmp[,j])
      }
      eval <- eval && (sum(testsum!=0)==1)
    }
  }

  if(eval){
    cat('The code group is an alternating code!\n')
  }else{
    cat('The code group is not an alternating code!\n')
  }

  return(eval)

} # ACtest

ACtest.strong <- function(cg){
# test if the given code is a strong type 1 alternating code

    if(!is.list(cg)) stop("the code sequence must be a list")
    ncode <- length(cg)
    nbit <- length(cg[[1]])

    codemat <- matrix(unlist(cg),nrow=ncode,byrow=T)

    test <- 0
    for(l1 in seq(nbit-1)){
        l2max <- 1
        if(l1==(nbit-1)) l2max <- 0
        for(l2 in seq(-1,l2max)){
            lmat1 <- codemat[,1:(nbit-l1)]*Conj(codemat[,(l1+1):nbit])
            lmat2 <- codemat[,1:(nbit-l2-l1)]*Conj(codemat[,(l2+l1+1):nbit])
            for(m in seq(dim(lmat1)[2])){
                for(n in seq(dim(lmat2)[2])){
                    if(!((l2==0)&(m==n))){
                        test <- test + abs(sum(lmat1[,m]*Conj(lmat2[,n])))
                    }
                }
            }
        }
    }

    if(test){
        cat("The code is not a strong alternating code.\n")
    }else{
        cat("The code is a strong alternating code.\n")
    }

  invisible(test)

} # ACtest.strong


#makeStrong <- function(codeOrig){
## convert a set of weak alternating codes into a
## strong one
#  code <- codeOrig
#  nc <- length(code)
#  nb <- length(code[[1]])
#  for (k in seq(1,nc)){
#    code[[k+nc]] <- code[[k]]
#    for (j in seq(2,nb,by=2)){
#      code[[k+nc]][j] <- -code[[k]][j]
#    }
#  }
#
#  return(code)
#
#} #makeStrong

multiplycg <- function(cg,coefs){
# multiply the bauds of the codes with the given  coefficients
   mcg <- list()
   nc <- length(cg)
   for (k in seq(1,nc)){
     mcg[[k]] <- coefs*cg[[k]]
   }

  return(mcg)

}

randomizeAC <- function(codeOrig){
# randomize a set of alternating codes
  nCode <- length(codeOrig)
  nBit <- length(codeOrig[[1]])
  code <- codeOrig
  rSeq <- ceiling(2*(runif(nBit)))*2 - 3
  for (k in seq(1,nCode)){
    for (m in seq(1,nBit)){
      code[[k]][m] <- code[[k]][m]*rSeq[m]
    }
  }

  return(code)

} #randomizeAC

## intToBin from R.utils..
intToBin <- function (x) 
{
    y <- as.integer(x)
    y <- format.binmode(y)
    y <- as.character(y)
    dim(y) <- dim(x)
    return(y)
}

format.binmode <- function (x, ...) 
{
    isna <- is.na(x)
    y <- x[!isna]
    ans0 <- character(length = length(y))
    neg <- which(y < 0)
    if (length(neg) > 0) {
        y[neg] <- y[neg] + 1L + .Machine$integer.max
    }
    z <- NULL
    while (any(y > 0) || is.null(z)) {
        z <- y%%2
        y <- floor(y/2)
        ans0 <- paste(z, ans0, sep = "")
    }
    ans <- rep(NA_character_, times = length(x))
    ans[!isna] <- ans0
    ans
}
