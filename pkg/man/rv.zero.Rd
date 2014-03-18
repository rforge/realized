\name{rv.zero}
\alias{rv.zero}
\title{ Calculates the percentage of zero returns at a specified sampling period }
\description{
Calculates the percentage of zero returns at a specified sampling period.
}
\usage{
rv.zero(x, period, align.by="seconds", align.period = 1, cts = TRUE, makeReturns = FALSE, ...)
}

\arguments{
    \item{x}{Tick data in xts object.}
  \item{period}{ Sampling period }
  \item{align.by}{Align the tick data to seconds|minutes|hours}
  \item{align.period}{ Align the returns to this period first}
  \item{cts}{ Create calendar time sampling if a non realizedObject is passed}
  \item{makeReturns}{ Prices are passed make them into log returns}
    \item{...}{ ...}
}
\value{
Percentage of zero returns.
}
\references{
S. W. Payseur.  A One Day Comparison of Realized Variance and Covariance Estimators. \emph{Working Paper: University 
of Washington}, 2007
}
\author{ Scott Payseur <scott.payseur@gmail.com>}

\seealso{ \code{\link{rSignature}}, \code{\link{rc.zero}} }


\examples{
data(sbux.xts)
rv.zero(x = sbux.xts, period = 60, align.by ="seconds", align.period=1)
}
\keyword{methods}