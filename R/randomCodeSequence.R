randomCodeSequence <- function( nCode=1 , nBits=1 , nPhases=2 , nAmplitudes=1)
  {
    # a random code sequence with nCode codes, each code with nBits bits, phases randomly chosen from 2pi/nPhases, 2*2pi/nPhases, ... , 2pi, and amplitudes randomly chosen from
    # 1/nAmplitudes, ... , 1

    code <- list()


    for ( k in seq( nCode ) )
      {
        code[[k]] <- exp( 2i * pi * floor( runif( nBits ) * nPhases ) / nPhases ) * ceiling( runif( nBits ) *  nAmplitudes  ) / nAmplitudes
      }
    

    return(code)

  } #randomCodeSequence


