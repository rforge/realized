\name{rv.kernel}
\alias{rv.kernel}

\title{ Realized Variance: Kernel}

\description{
Realized variance calculation using a kernel estimator. 
}


\usage{
rv.kernel(x, kernel.type = "rectangular", kernel.param = 1, kernel.dofadj = TRUE, align.by="seconds",align.period = 1, cts = TRUE, makeReturns = FALSE)
}

\arguments{
  \item{x}{Tick data in xts object.}
  \item{kernel.type}{Kernel name (or number)}
  \item{kernel.param}{Kernel parameter (usually lags)}
  \item{kernel.dofadj}{Kernel Degree of freedom adjustment}
  \item{align.by}{Align the tick data to seconds|minutes|hours}
  \item{align.period}{Align the tick data to this many [seconds|minutes|hours]}
  \item{cts}{Calendar Time Sampling is used}
  \item{makeReturns}{Convert to Returns}
}

\details{
The different types of kernels can be found using rKernel.available().
}
\value{
Kernel estimate of realized variance.
}
\references{
Ole E. Barndorff-Nielsen, Peter Reinhard Hansen, Asger Lunde, and Neil Shephard. Regular and modified kernel-based estimators of integrated variance: The case with independent noise. \emph{Working Paper}, 2004.

B. Zhou. High-frequency data and volatility in foreign-exchange rates. \emph{Journal of Buiness & Economic
Statistics}, 14:45-52, 1996.

P. Hansen and A. Lunde. Realized variance and market microstructure noise. \emph{Journal of Business and
Economic Statistics}, 24:127-218, 2006.
}
\author{Scott Payseur <spayseur@u.washington.edu>}

\seealso{\code{\link{rRealizedVariance}},\code{\link{rc.kernel}}, \code{\link{rKernel.available}},\code{\link{rKernel}}}
\examples{

data(msft.real.cts)
data(ge.real.cts)

# kernel realized variance for CTS aligned at one minute returns
#
rv.kernel(x = msft.real.cts[[1]], kernel.param=1, kernel.type="bartlett", align.period=60)

}
\keyword{methods}
