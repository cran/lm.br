\name{cr}
\alias{cr,Cpp_Clmbr-method}
\alias{cr}

\title{Confidence Region for the Changepoint}

\description{
  Joint confidence region for ( theta, alpha ),  the changepoint's (x,y)-coordinates.
}

\usage{
\S4method{cr}{Cpp_Clmbr}(  CL =0.95 ,  method ="CLR",  incr,  output ="G"  )
}

\arguments{
  \item{CL}{confidence level, between 0 and 1}
  \item{method}{"CLR" or "AF" which stand for conditional likelihood-ratio  or  approximate-F (rapid), see \code{\link{sl}} for details}
  \item{incr}{increment of theta values in the confidence region boundary-points}
  \item{output}{"G", "T" or "V" which stand for graph, text printout  or  value}
}

\details{
This subroutine scans to determine the postulate values of (theta,alpha) that have significance level greater than 1-CL.  It scans first along the (theta, alpha-MLE) ridge to determine the 'theta' limit-boundaries.
}

\value{
If 'output' is "V", then 'cr' returns an N x 3 matrix of boundary points (theta, min-alpha, max-alpha).  Otherwise, 'cr' graphs or prints-out the confidence region, but does not return a value.
}


\examples{
##  Data for Patient B from Smith and Cook (1980)
library(lm.br)
y <- c(37.3, 47.1, 51.5, 67.6, 75.9, 73.3, 69.4, 61.5, 31.8, 19.4)
x <- 1:10
sc <- lm.br(y~x)
\donttest{sc$cr()}
\donttest{sc $ cr( .9, incr = 0.1 )}
sc $ cr( .99, 'af', out='t' )
}

