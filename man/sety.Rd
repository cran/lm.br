\name{sety}
\alias{sety,Cpp_Clmbr-method}
\alias{sety}

\title{Set y-Values}

\description{
  Reset the response values in the C++ object.  
}

\usage{
\S4method{sety}{Cpp_Clmbr}( rWy )
}

\arguments{
\item{rWy}{vector of 'y' values, pre-multiplied by the square-root of 'weights'}
}

\details{
The pre-multiplied vector is more convenient as input during simulation tests.  The 'rWy' vector is simply the y-vector if the model does not include 'weights'.  

The square-root of a vector 'W' is the vector 'rW' of the square-roots of the elements of 'W'.  The square-root of a matrix 'W' is the matrix 'rW' such that rW*rW = W (a stricter definition than  rW*transpose(rW) = W).
}

\note{
This method changes the y-values for the accessor functions 'sl', 'ci', 'cr' and 'mle'.  But the output list of 'coefficients', 'fitted.values' and so on retains its values from the original call.
}


\examples{
#  Data for Patient B from Smith and Cook (1980)
y <- c( 36.0, 45.5, 50.0, 60.0, 73.3, 71.0, 66.7, 60.0, 30.5, 18.3 )
x <- 1:10
sc <- lm.br( y ~ x )
sc$ci()
newy <- c( 37.3, 47.1, 51.5, 67.6, 75.9, 73.3, 69.4, 61.5, 31.8, 19.4 )
sc$sety( newy )
sc$ci()



#  An example simulation test

simtest  <-  function( x, theta, alpha, B, Bp, var, N =10000 )  {
  cat("\n\n   Monte Carlo simulation test:\n")
  cat("for model   theta=",theta," alpha=",alpha," B=",B," B'=",Bp," var=",var,"\n")
  cat("with observations at  x=",x,"\n\n")
  cat("     no. of     coverage frequency of the 0.95-confidence interval by\n")
  cat("   iterations                CLR                    AF\n")
  flush.console()

  y <- x
  LLmodel <- lm.br( y ~ x )    # construct an object with arbitrary 'y'
  n <- length(x)
  countCLR <- countAF <- 0
  for( i in 1:N )  {
    err <- rnorm( n, 0, sqrt(var) )
    y <- alpha + B*pmin(x-theta,0) + Bp*pmax(x-theta,0) + err
    LLmodel$sety( y )
    stest <- LLmodel$sl( theta, 'clr', .0001, FALSE )
    if( stest > 0.05 )  countCLR <- countCLR + 1
    stest <- LLmodel$sl( theta, 'af', .0001, FALSE )
    if( stest > 0.05 )  countAF <- countAF + 1

    if( floor(i/1000) - i/1000 == 0) {
      cat( format( i, width=10 ),
        format( countCLR/i, digits=4, nsmall=4, width=22 ),
        format(  countAF/i, digits=4, nsmall=4, width=22 ), "\n" )
      flush.console()
    }
  }
}

x <- c( 1.0, 1.1, 1.3, 1.7, 2.4, 3.9, 5.7, 7.6, 8.4, 8.6 )
\donttest{simtest( x, 3, 0, -1, 0.5, 1 )}
}

