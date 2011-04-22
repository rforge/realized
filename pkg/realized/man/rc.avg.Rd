\name{rc.avg}
\alias{rc.avg}

\title{Realized Covariance: Average Subsample}
\description{
Realized Covariance using average subsample.
}
\usage{
rc.avg(x, y, period, align.by="seconds",align.period = 1, cts = TRUE, makeReturns = FALSE, ...)
}
\arguments{  
  \item{x}{Tick data in xts object.}
  \item{y}{Tick data in xts object.}
  \item{period}{ Sampling period }
  \item{align.by}{Align the tick data to seconds|minutes|hours}
  \item{align.period}{Align the tick data to this many [seconds|minutes|hours]}
  \item{cts}{ Create calendar time sampling if a non realizedObject is passed}
  \item{makeReturns}{ Prices are passed make them into log returns}
  \item{...}{ ...}
}
\value{
Realized covariance using average subsample.
}
\references{ 
L. Zhang, P.A Mykland, and Y. Ait-Sahalia. A tale of two time scales: Determining integrated volatility
with noisy high-frequency data. \emph{Journal of the American Statistical Association}, 2005.

Michiel de Pooter, Martin Martens, and Dick van Dijk. Predicting the daily covariance matrix for sp100
stocks using intraday data - but which frequency to use? \emph{Working Paper}, October 2005.
 }
\author{ Scott Payseur <scott.payseur@gmail.com> }
\seealso{\code{\link{rv.avg}}, \code{\link{rRealizedVariance}} }
\examples{


data(sbux.xts)
data(lltc.xts)


#
# Average subsampled realized variance for CTS aligned at one second returns at 
# 600 subgrids (10 minutes).
#
rc.avg(x = sbux.xts, y=lltc.xts, period = 600, align.by ="seconds", align.period=1)

#
# Average subsampled realized variance for CTS aligned at one minute returns at 
# 5 subgrids (5 minutes).
#
rc.avg(x = sbux.xts,y=lltc.xts, period = 5, align.by ="minutes", align.period=5)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
