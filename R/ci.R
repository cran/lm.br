
# this function is internal, not meant for the user

# R code to process user input, then call
# the corresponding function of the C++ object


.ci <- function( z, CL =0.95, method ="clr", output ="T" )  {

  method <- toupper(method)
  met <- integer(1)
  if( method=="CLR" )  met <- 1  else  {
    if( method=="AF" )  met <- 2  else
      stop( "'method' must be \"CLR\" or \"AF\"" )
  }  

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

  if(value) {
    bounds <- (z$CppObj)$ci2( CL, met, as.integer(verbose) )
    return( bounds )
  }  else  {
    (z$CppObj)$ci( CL, met )
  }
  
}

