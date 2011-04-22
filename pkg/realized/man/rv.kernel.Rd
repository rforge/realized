\name{rv.kernel}
\alias{rv.kernel}

\title{ Realized Variance: Kernel}

\description{
Realized variance calculation using a kernel estimator. 
}


\usage{
rv.kernel(x, kernel.type = "rectangular", kernel.param = 1, kernel.dofadj = TRUE, align.by="seconds",align.period = 1, cts = TRUE, makeReturns = FALSE, type = NULL, adj = NULL, q = NULL,...)
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
   \item{type}{Deprecated, use kernel.type}
    \item{adj}{Deprecated, use kernel.dofadj}
  \item{q}{Deprecated, use kernel.param}
  \item{...}{...}
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
\author{Scott Payseur <scott.payseur@gmail.com>}

\seealso{\code{\link{rRealizedVariance}},\code{\link{rc.kernel}}, \code{\link{rKernel.available}},\code{\link{rKernel}}}
\examples{
# kernel realized variance for CTS aligned at one minute returns
data(sbux.xts)
rv.kernel(x = sbux.xts, kernel.param=1, kernel.type="bartlett", align.by ="seconds", align.period=60)
}
\keyword{methods}
