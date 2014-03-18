\name{rRealizedVariance}
\alias{rRealizedVariance}
\alias{rRealisedVariance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Calculate realized variance, covariance, or correlation.}
\description{
Calculate realized variance, covariance, correlation, covariance matrix, or correlation matrix.
}
\usage{
rRealizedVariance(x, y = NULL, type = "naive", period = 1,align.by="seconds", align.period = 1, cor = FALSE, rvargs = list(), cts = TRUE, makeReturns = FALSE, lags=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{x}{Tick data in xts object.}
   \item{y}{Tick data in xts object.}
  \item{period}{ Sampling period }
  \item{align.by}{Align the tick data to seconds|minutes|hours}
    \item{align.period}{Align the tick data to this many [seconds|minutes|hours]}
  \item{type}{ Type of realized estimator to use, a rv. or rc. is appended to this value and that function is called }
  \item{cor}{ T for correlation }
  \item{rvargs}{ List of extra parameters to pass into rv.* or rc.* }
  \item{cts}{ Create calendar time sampling if a non realizedObject is passed}
  \item{makeReturns}{ Prices are passed make them into log returns}
  \item{lags}{Deprecated}
}
\details{
Calculate realized variance, covariance, correlation, covariance matrix, or correlation matrix.
}
\value{
A single numeric value or a matrix if x is multicolumn matrix.
}
\author{ Scott Payseur <scott.payseur@gmail.com>}


\seealso{ \code{\link{rc.avg}}, \code{\link{rc.kernel}}, \code{\link{rc.naive}}, \code{\link{rc.timescale}}, \code{\link{rv.avg}}, \code{\link{rv.kernel}}, \code{\link{rv.naive}}, \code{\link{rv.timescale}}  }
\examples{

data(sbux.xts)
data(lltc.xts)

#
#  VARIANCE
#

# Traditional Estimate at highest frequency
rRealizedVariance(x=sbux.xts, type="naive", period=1, align.by="seconds", align.period=1)

# Traditional Estimate at one minute frequency 
rRealizedVariance(x=sbux.xts, type="naive", period=1, align.by="seconds", align.period=60) 
rRealizedVariance(x=sbux.xts, type="naive", period=1, align.by="mins", align.period=1) 

# Traditional Estimate at 10 minute frequency 
rRealizedVariance(x=sbux.xts, type="naive", period=1, align.by="mins", align.period=10) 

# Bartlett Kernel Estimate with minute aligned data at 20 lags 
rRealizedVariance(x=sbux.xts, type="kernel", align.by="mins", align.period=1, rvargs=list(kernel.param=20,kernel.type="Bartlett"))

# Cubic Kernel Estimate with second aligned data at 400 lags 
rRealizedVariance(x=sbux.xts, type="kernel", align.by="seconds", align.period=1, rvargs=list(kernel.param=400, kernel.type="Cubic"))
 
# Subsample Average Estimate with second aligned data at 600 subgrids 
rRealizedVariance(x=sbux.xts, type="avg", period=600) 

# Timescale Average Estimate with second aligned data at 600 subgrids 
rRealizedVariance(x=sbux.xts, type="timescale", period=600, rvargs=list(adj.type="aa")) 


#
# COVARIANCE
# 

# Traditional Estimate at highest frequency
rRealizedVariance(x=sbux.xts, y=lltc.xts, type="naive", period=1, align.by="seconds", align.period=1)

# Traditional Estimate at one minute frequency 
rRealizedVariance(x=sbux.xts, y=lltc.xts, type="naive", period=1, align.by="seconds", align.period=60) 
rRealizedVariance(x=sbux.xts, y=lltc.xts, type="naive", period=1, align.by="mins", align.period=1) 

# Traditional Estimate at 10 minute frequency 
rRealizedVariance(x=sbux.xts, y=lltc.xts, type="naive", period=1, align.by="mins", align.period=10) 

# Bartlett Kernel Estimate with minute aligned data at 20 lags 
rRealizedVariance(x=sbux.xts, y=lltc.xts, type="kernel", align.by="mins", align.period=1, rvargs=list(kernel.param=20,kernel.type="Bartlett"))

# Cubic Kernel Estimate with second aligned data at 400 lags 
rRealizedVariance(x=sbux.xts, y=lltc.xts, type="kernel", align.by="seconds", align.period=1, rvargs=list(kernel.param=400, kernel.type="Cubic"))
 
# Subsample Average Estimate with second aligned data at 600 subgrids 
rRealizedVariance(x=sbux.xts, y=lltc.xts, type="avg", period=600) 

# Timescale Estimate with second aligned data at 600 subgrids 
rRealizedVariance(x=sbux.xts, y=lltc.xts, type="timescale", period=600, rvargs=list(adj.type="aa")) 


#
# CORRELATION
#
 
# Traditional Estimate at highest frequency
rRealizedVariance(x=sbux.xts, y=lltc.xts, type="naive", period=1, align.by="seconds", align.period=1, cor=TRUE)

# Traditional Estimate at one minute frequency 
rRealizedVariance(x=sbux.xts, y=lltc.xts, type="naive", period=1, align.by="seconds", align.period=60, cor=TRUE) 
rRealizedVariance(x=sbux.xts, y=lltc.xts, type="naive", period=1, align.by="mins", align.period=1, cor=TRUE) 

# Traditional Estimate at 10 minute frequency 
rRealizedVariance(x=sbux.xts, y=lltc.xts, type="naive", period=1, align.by="mins", align.period=10, cor=TRUE) 

# Bartlett Kernel Estimate with minute aligned data at 20 lags 
rRealizedVariance(x=sbux.xts, y=lltc.xts, type="kernel", align.by="mins", align.period=1, rvargs=list(kernel.param=20,kernel.type="Bartlett"), cor=TRUE)

# Cubic Kernel Estimate with second aligned data at 400 lags 
rRealizedVariance(x=sbux.xts, y=lltc.xts, type="kernel", align.by="seconds", align.period=1, rvargs=list(kernel.param=400, kernel.type="Cubic"), cor=TRUE)
 
# Subsample Average Estimate with second aligned data at 600 subgrids 
rRealizedVariance(x=sbux.xts, y=lltc.xts, type="avg", period=600, cor=TRUE) 

# Timescale Estimate with second aligned data at 600 subgrids 
rRealizedVariance(x=sbux.xts, y=lltc.xts,type="timescale", period=600, rvargs=list(adj.type="aa"), cor=TRUE) 
  

#
# COVARIANCE MATRIX
#
rRealizedVariance(x=merge(sbux.xts,lltc.xts, fill=0),type="timescale", period=600, rvargs=list(adj.type="aa"), cor=FALSE) 

#
# CORRELATION MATRIX
#
rRealizedVariance(x=merge(sbux.xts,lltc.xts, lltc.xts, fill=0), type="naive", period=1, align.by="mins", align.period=10, cor=TRUE) 
 

}

\keyword{methods}