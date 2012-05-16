\name{rv.naive}
\alias{rv.naive}
\title{ Realized Variance }
\description{
  Traditional realized variance estimator.
}
\usage{
rv.naive(x, period, align.by = "seconds", align.period = 1, cts = TRUE, makeReturns = FALSE, ...)
}

\arguments{
  \item{x}{Tick data in xts object.}
  \item{period}{ Sampling period }
  \item{align.by}{Align the tick data to seconds|minutes|hours}
  \item{align.period}{Align the tick data to this many [seconds|minutes|hours]}
  \item{cts}{Calendar Time Sampling is used}
  \item{makeReturns}{Convert to Returns}
  \item{...}{...}
  }

\value{
Sum of products of high frequency returns.
}
\references{

O. E. Barndorff-Nielson and N. Sheppard. Econometric analysis of realised volatility and its use in estimating
stochastic volatility models. \emph{Journal of the Royal Statistical Society}, Series B 64:253-280, 2002

T.G. Andersen, T. Bollerslev, F.X. Diebold, and P. Labys. The distribution of exchange rate volatility.
\emph{Journal of the American Statistical Association}, 96:42-55, 2001.

}
\author{Scott Payseur <scott.payseur@gmail.com> }

\seealso{\code{\link{rc.naive}}, \code{\link{rSignature}}, \code{\link{rRealizedVariance}}}
\examples{


data(msft.real.cts)
data(ge.real.cts)

# Traditional RV aligned to 60 seconds
#
rv.naive(x = msft.real.cts[[1]], period = 60)

}
\keyword{methods}