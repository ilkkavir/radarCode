makeStrongCode <- function(code)
  {
    # convert a weak (alternating) code sequence into a strong one
    #
    #
    # INPUT:
    #   code list( c(code1[1],code1[2],...,code1[N]) , c(code2[1],...,code2[N]) , ... , c( codeM[1],...,codeM[N]) )
    #        a list containing code phases and amplitudes as complex values
    #
    # OUTPUT:
    #   strongCode  a strong version of the code. The original setis repeated twice, and every second bit in the copied part is multiplied with -1

  nc <- length(code)
  nb <- length(code[[1]])
  for (k in seq(1,nc)){
    code[[k+nc]] <- code[[k]]
    for (j in seq(2,nb,by=2)){
      code[[k+nc]][j] <- -code[[k]][j]
    }
  }
  
  return(code)

} #makeStrongCode
