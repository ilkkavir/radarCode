lengthenBits <- function(code,N){
# change the given code into the form where
# we have N samples from each bit

newCode <- list()
nCodes <- length(code)

for (k in seq(1,nCodes)){
newCode[[k]] <- rep(0,(N*length(code[[k]])))
  for (j in seq(1,length(code[[k]]))){
    newCode[[k]][((j-1)*N+1):(j*N)] <- code[[k]][j]
  }

}

return(newCode)

}

