\name{rc.zero}
\alias{rc.zero}
\title{ Calculates the percentage of co-zero returns at a specified sampling period }
\description{
Calculates the percentage of co-zero returns at a specified sampling period.
}
\usage{
rc.zero(x, y, period, align.by="seconds", align.period = 1, cts = TRUE, makeReturns = FALSE, ...)
}

\arguments{
   \item{x}{Tick data in xts object.}
  \item{y}{Tick data in xts object.}
  \item{period}{ Sampling period }
  \item{align.by}{Align the tick data to seconds|minutes|hours}
  \item{align.period}{ Align the returns to this period first}
  \item{cts}{ Create calendar time sampling if a non realizedObject is passed}
  \item{makeReturns}{ Prices are passed make them into log returns}
   \item{...}{ ...}
}
\value{
Percentage of co-zero returns.
}
\references{
S. W. Payseur.  A One Day Comparison of Realized Variance and Covariance Estimators. \emph{Working Paper: University 
of Washington}, 2007
}
\author{ Scott Payseur <scott.payseur@gmail.com>}

\seealso{ \code{\link{rSignature}}, \code{\link{rv.zero}} }


\examples{
data(sbux.xts)
data(lltc.xts)

rc.zero(x = sbux.xts, y=lltc.xts, period = 60, align.by ="seconds", align.period=1)
}
\keyword{methods}
