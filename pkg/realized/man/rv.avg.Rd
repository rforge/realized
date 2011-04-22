\name{rv.avg}
\alias{rv.avg}

\title{Realized Variance: Average Subsample}
\description{
  Realized Variance using average subsample.
}
\usage{
rv.avg(x,period, align.by="seconds", align.period = 1, cts = TRUE, makeReturns = FALSE, ...)
}

\arguments{
   \item{x}{Tick data in xts object.}
  \item{period}{ Sampling period }
  \item{align.by}{Align the tick data to seconds|minutes|hours}
  \item{align.period}{Align the tick data to this many [seconds|minutes|hours]}
  \item{cts}{ Create calendar time sampling if a non realizedObject is passed}
  \item{makeReturns}{ Prices are passed make them into log returns}
    \item{...}{ ...}
}
\value{
 Realized variance using average subsample.
}
\references{
L. Zhang, P.A Mykland, and Y. Ait-Sahalia. A tale of two time scales: Determining integrated volatility
with noisy high-frequency data. \emph{Journal of the American Statistical Association}, 2005.

}
\author{ Scott Payseur <spayseur@u.washington.edu> }
\seealso{\code{\link{rv.avg}}, \code{\link{rRealizedVariance}} }
\examples{

data(sbux.xts)


#
# Average subsampled realized variance for CTS aligned at one second returns at 
# 600 subgrids (10 minutes).
#
rv.avg(x = sbux.xts, period = 600, align.by ="seconds", align.period=1)

#
# Average subsampled realized variance for CTS aligned at one minute returns at 
# 5 subgrids (5 minutes).
#
rv.avg(x = sbux.xts, period = 5, align.by ="minutes", align.period=5)


}
\keyword{methods}