\name{rc.naive}
\alias{rc.naive}
\title{ Realized Covariance }
\description{
  Traditional realized covariance estimator.
}
\usage{
rc.naive(x, y, period, align.by = "seconds", align.period = 1, cts = TRUE, makeReturns = FALSE, ...)
}

\arguments{
  \item{x}{Tick data in xts object.}
  \item{y}{Tick data in xts object.}
  \item{period}{ Sampling period }
  \item{align.by}{Align the tick data to seconds|minutes|hours}
  \item{align.period}{Align the tick data to this many [seconds|minutes|hours]}
  \item{cts}{Calendar Time Sampling is used}
  \item{makeReturns}{Convert to Returns}
  \item{...}{...}
  }

\value{
Sum of cross products of high frequency returns.
}
\references{
T.G. Andersen, T. Bollerslev, F.X. Diebold, and P. Labys. The distribution of exchange rate volatility.
\emph{Journal of the American Statistical Association}, 96:42-55, 2001.

}
\author{Scott Payseur <scott.payseur@gmail.com> }

\seealso{\code{\link{rv.naive}}, \code{\link{rSignature}}, \code{\link{rRealizedVariance}}}
\examples{
data(sbux.xts)
data(lltc.xts)

# Traditional RC aligned to 60 seconds
#
rc.naive(x = sbux.xts, y=lltc.xts, period = 60, align.by ="seconds", align.period=1)

}
\keyword{methods}
