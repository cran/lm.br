
# this function is internal, not meant for the user

# R code to process user input, then call
# the corresponding function of the C++ object


.mle <- function( z, output ="T" )  {

  value <- verbose <- logical(1)
  output <- toupper(output)
  if( output=="T" )  {
    value <- FALSE
    verbose <- TRUE  
  } else {
    if( output=="V" )  {
      value <- TRUE
      verbose <- FALSE
    } else {
      if( output=="B" )
        value <- verbose <- TRUE
      else  stop( "'output' must be \"T\", \"V\" or \"B\"" )
    }
  }

  if(verbose) {
    (z$CppObj)$mle( )
  }
  
  if(value) {
    par <- (z$CppObj)$param()
    mles <- c( par[1], par[2], par[3], par[4], par[5] )  
    return( mles )
  }
  
}

