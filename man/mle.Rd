\name{mle}
\alias{mle,Cpp_Clmbr-method}
\alias{mle}

\title{Maximum Likelihood Estimates}

\description{
  Maximum-likelihood estimates of parameters.  Estimates are without bias correction except for the variance.
}

\usage{
\S4method{mle}{Cpp_Clmbr}( output ="T" )
}

\arguments{
  \item{output}{"T", "V" or "B" which stand for text, value or both.}
}

\value{
'mle' prints-out the maximum-likelihood estimates but does not return a value if 'output' is "T".  'mle' returns a numeric vector of maximum-likelihood estimates if 'output' is "V" or "B".
}


\examples{
#  Data for Patient B from Smith and Cook (1980)
y <- c(37.3, 47.1, 51.5, 67.6, 75.9, 73.3, 69.4, 61.5, 31.8, 19.4)
x <- 1:10
sc <- lm.br(y~x)
sc$mle()
estimates <- sc$mle( 'v' )
estimates
}

