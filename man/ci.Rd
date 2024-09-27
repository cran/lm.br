\name{ci}
\alias{ci,Cpp_Clmbr-method}
\alias{ci}

\title{Confidence Interval for the Changepoint}

\description{
  Confidence interval for 'theta', the changepoint's x-coordinate.
}

\usage{
\S4method{ci}{Cpp_Clmbr}(  CL =0.95, method ="CLR", output ="T" )
}

\arguments{
  \item{CL}{confidence level, between 0 and 1.}
  \item{method}{"CLR" or "AF" which stand for conditional likelihood-ratio or approximate-F, see \code{\link{sl}} for details.}
  \item{output}{"T", "V" or "B" which stand for text, value or both.}
}

\details{
This subroutine scans to determine the postulate values of 'theta' that have significance level greater than 1-CL.
}

\value{
'ci' prints-out the confidence interval for 'theta' but does not return a value if 'output' is "T".  'sl' returns a numeric vector of boudaries for the contiguous segments of the confidence interval if 'output' is "V" or "B".
}

\examples{
#  Data for Patient B from Smith and Cook (1980)
y <- c(37.3, 47.1, 51.5, 67.6, 75.9, 73.3, 69.4, 61.5, 31.8, 19.4)
x <- 1:10
sc <- lm.br( y ~ x )
sc$ci()
sc $ ci( 0.90 )
sc $ ci( .99, 'af' )
sc $ ci( out= 'v' )
}

