#.First.lib <- function(libname, pkgname)
#{
#    library.dynam("realized")
#    require("xts")
#    require("zoo")
#}

.onLoad <- function(libname, pkgname)
{
#    require("xts")
#    require("zoo")
}

.onAttach <- function(libname, pkgname)
{
 #   require("xts")
  #  require("zoo")
}



     #########################################################################
     #
     # Utility Functions
     #
     #########################################################################
     .alignedAccum <- function(x,y, period, cum=TRUE, makeReturns...)
     {
          x<-.accum.naive(x,x, period)
          y<-.accum.naive(y,y, period)
          if(cum)
          {
               ans <- cumsum(x*y)
          }
          else
          {
              ans <- x*y     
          }
          ans
     }


     .accum.naive <- function(x,y, period, ...)
     {
          .C("rv", 
                    as.double(x), #a
                    as.double(y), #b
                    as.integer(length(x)), #na
                    as.integer(period), #period 
                    tmpa = as.double(rep(0,as.integer(length(x)/period +1))), #tmp
                    as.double(rep(0,as.integer(length(x)/period +1))), #tmp
                    as.integer(length(x)/period), #tmpn
                    ans = double(1), 
                    COPY=c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,TRUE), 
                    PACKAGE="realized")$tmpa     
     }


     .alignReturns <- function(x, period, ...)
     {
          .C("rv", 
             as.double(x), #a
             as.double(x), #b
             as.integer(length(x)), #na
             as.integer(period), #period 
             tmpa = as.double(rep(0,as.integer(length(x)/period +1))), #tmp
             as.double(rep(0,as.integer(length(x)/period +1))), #tmp
             as.integer(length(x)/period), #tmpn
             ans = double(1), 
             COPY=c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,TRUE), 
             PACKAGE="realized")$tmpa     
     }

     .getAlignPeriod <- function(align.period, align.by)
     {   
         align.by <- gsub("(^ +)|( +$)", "",align.by) # Trim White
         
         if(casefold(align.by)=="min" || casefold(align.by)=="mins" ||casefold(align.by)=="minute"||casefold(align.by)=="minutes"||casefold(align.by)=="m"){
             ans <- align.period * 60
         }
         if(casefold(align.by)=="sec" || casefold(align.by)=="secs" ||casefold(align.by)=="second"||casefold(align.by)=="seconds"||casefold(align.by)=="s"||casefold(align.by)==""){
	     ans <- align.period
	 }
         if(casefold(align.by)=="hour" || casefold(align.by)=="hours" ||casefold(align.by)=="h"){
	     ans <- align.period * 60 * 60
	 }
         return(ans)
     }


     .alignIndices <- function(x, period, ...)
     {
          .C("rvperiod", 
             as.double(x), #a
             as.double(x), #b
             as.integer(length(x)), #na
             as.integer(period), #period 
             tmpa = as.double(rep(max(x),as.integer(length(x)/period +1))), #tmp
             as.double(rep(0,as.integer(length(x)/period +1))), #tmp
             as.integer(length(x)/period), #tmpn
             ans = double(1), 
             COPY=c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,TRUE), 
             PACKAGE="realized")$tmpa     
     }




     #########################################################################
     #
     # Kernel Estimators 
     # See BNHLS (2006), Zhou (1996), HL (1994), HL (1996)
     #
     #########################################################################


     rKernel <- function(x,type=0)
     {
          type <- .kernel.chartoint(type)
         .C("justKernel", x=as.double(x),type= as.integer(type), ans=as.double(0),PACKAGE="realized")$ans
     }

     .kernel.chartoint <- function(type)
     {
        if(is.character(type))
         {
          ans <- switch(casefold(type), 
                 rectangular=0,
                 bartlett=1,
                 second=2,
                 epanechnikov=3,
                 cubic=4,
                 fifth=5,
                 sixth=6,
                 seventh=7,
                 eighth=8,
                 parzen=9,
                 th=10,
                 mth=11,
                 tukeyhanning=10,
                 modifiedtukeyhanning=11,
                 -99)
           if(ans==-99)
           { 
               warning("Invalid Kernel, using Bartlet")
               1
           }
           else
           {
               ans     
           }
         }
         else
         {
          type
         }
     }

     rKernel.available <- function()
     {
          c("Rectangular", 
            "Bartlett",
            "Second",
            "Epanechnikov",
            "Cubic",
            "Fifth",
            "Sixth",
            "Seventh",
            "Eighth",
            "Parzen",
            "TukeyHanning",
            "ModifiedTukeyHanning")
     }

     # 
     # Realized variance calculation using a kernel estimator.
     #
     rv.kernel <- function(x,                             # Tick Data
                           kernel.type = "rectangular",   # Kernel name (or number)
                           kernel.param = 1,              # Kernel parameter (usually lags)
                           kernel.dofadj = TRUE,          # Kernel Degree of freedom adjustment
                           align.by="seconds",            # Align the tick data to [seconds|minutes|hours]
                           align.period = 1,              # Align the tick data to this many [seconds|minutes|hours]
                           cts = TRUE,                    # Calendar Time Sampling is used
                           makeReturns = FALSE,            # Convert to Returns 
                           type = NULL,                   # Deprectated
                           adj = NULL,                    # Deprectated
                           q = NULL, ...){                     # Deprectated
                           
          #
          # Handle deprication
          #
          if(!is.null(type)){
              warning("type is deprecated, use kernel.type")
              kernel.type=type
          }
          if(!is.null(q)){
	                warning("q is deprecated, use kernel.param")
	                kernel.param=q
          }
          if(!is.null(adj)){
	                warning("adj is deprecated, use kernel.dofadj")
	                kernel.dofadj=adj
          }          
          align.period = .getAlignPeriod(align.period, align.by)         
          cdata <- .convertData(x, cts=cts, makeReturns=makeReturns)
          x <- cdata$data
          x <- .alignReturns(x, align.period)
          type <- .kernel.chartoint(kernel.type)
          .C("kernelEstimator", as.double(x), as.double(x), as.integer(length(x)),
                     as.integer(kernel.param), as.integer(ifelse(kernel.dofadj, 1, 0)),
                     as.integer(type), ab=double(kernel.param + 1),
                     ab2=double(kernel.param + 1),
                     ans=double(1),PACKAGE="realized")$ans
     }

     rc.kernel <- function(x,                             # Tick Data for first asset
                           y,                             # Tick Data for second asset
                           kernel.type = "rectangular",   # Kernel name (or number)
                           kernel.param = 1,              # Kernel parameter (usually lags)
                           kernel.dofadj = TRUE,          # Kernel Degree of freedom adjustment
                           align.by="seconds",            # Align the tick data to [seconds|minutes|hours]
                           align.period = 1,              # Align the tick data to this many [seconds|minutes|hours]
                           cts = TRUE,                    # Calendar Time Sampling is used
                           makeReturns = FALSE,           # Convert to Returns 
                           type = NULL,                   # Deprectated
                           adj = NULL,                    # Deprectated
                           q = NULL,...){                     # Deprectated
          #
          # Handle deprication
          #
          if(!is.null(type)){
              warning("type is deprecated, use kernel.type")
              kernel.type=type
          }
          if(!is.null(q)){
	                warning("q is deprecated, use kernel.param")
	                kernel.param=q
          }
          if(!is.null(adj)){
	                warning("adj is deprecated, use kernel.dofadj")
	                kernel.dofadj=adj
          }
          
          align.period = .getAlignPeriod(align.period, align.by)   
          cdata <- .convertData(x, cts=cts, makeReturns=makeReturns)
          
          x <- cdata$data
          x <- .alignReturns(x, align.period)
          cdatay <- .convertData(y, cts=cts, makeReturns=makeReturns)
          y <- cdatay$data
          y <- .alignReturns(y, align.period)
          type <- .kernel.chartoint(kernel.type)
          .C("kernelEstimator", as.double(x), as.double(y), as.integer(length(x)),
                     as.integer(kernel.param), as.integer(ifelse(kernel.dofadj, 1, 0)),
                     as.integer(type), ab=double(kernel.param + 1),
                     ab2=double(kernel.param + 1),
                     ans=double(1),PACKAGE="realized")$ans
     }


     #########################################################################
     #
     # Subsample based estimators
     # See AMZ (),(),(), 
     #
     #########################################################################

     .rv.subsample <- function(x, period, cts=TRUE, makeReturns=FALSE,...)
     {
          cdata <- .convertData(x, cts=cts, makeReturns=makeReturns)
          x <- cdata$data

          .C("subsample", 

                    as.double(x), #a
                    as.double(x), #na
                    as.integer(length(x)), #na
                    as.integer(length(x)/period),       #m
                    as.integer(period), #period 
                    as.double(rep(0,as.integer(length(x)/period +1))), #tmp
                    as.double(rep(0,as.integer(length(x)/period +1))), #tmp
                    as.integer(length(x)/period), #tmpn
                    ans = double(period), 
                    COPY=c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,TRUE), 
                    PACKAGE="realized")$ans
     }


     .rc.subsample <- function(x, y, period, cts=TRUE, makeReturns=FALSE, ... )
     {
          cdata <- .convertData(x, cts=cts, makeReturns=makeReturns)
          x <- cdata$data

          cdatay <- .convertData(y, cts=cts, makeReturns=makeReturns)
          y <- cdatay$data

          .C("subsample", 
                    as.double(x), #a
                    as.double(y), #na
                    as.integer(length(x)), #na
                    as.integer(length(x)/period),       #m
                    as.integer(period), #period 
                    as.double(rep(0,as.integer(length(x)/period +1))), #tmp
                    as.double(rep(0,as.integer(length(x)/period +1))), #tmp
                    as.integer(length(x)/period), #tmpn
                    ans = double(period), 
                    COPY=c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,TRUE), 
                    PACKAGE="realized")$ans               
     }

     rv.timescale <- function(x, period, align.by="seconds", align.period=1,adj.type="classic", cts=TRUE, makeReturns=FALSE, ...)
     {
          align.period = .getAlignPeriod(align.period, align.by)
          x<- .alignReturns(.convertData(x, cts=cts, makeReturns=makeReturns)$data, align.period)

          n <- dim(as.matrix(x))[[1]]
          nbar <- (n-period+1)/(period)
          adj <- switch(adj.type, classic=1, adj=(1-(nbar/n))^-1, aa= n/((period-1)*nbar))
          adj * (mean(rv.avg(x, period)) - ((nbar/n) * rv.naive(x,1)))
     }

     rc.timescale <- function(x,y, period, align.by="seconds", align.period=1, adj.type="classic", cts=TRUE, makeReturns=FALSE,...)
     {
          align.period = .getAlignPeriod(align.period, align.by)
          x<- .alignReturns(.convertData(x, cts=cts, makeReturns=makeReturns)$data, align.period)
          y<- .alignReturns(.convertData(y, cts=cts, makeReturns=makeReturns)$data, align.period)

          n <- dim(as.matrix(x))[[1]]
          nbar <- (n-period+1)/(period)
          adj <- switch(adj.type, classic=1, adj=(1-(nbar/n))^-1, aa= n/((period-1)*nbar))
          adj * (mean(rc.avg(x,y, period)) - ((nbar/n) * rc.naive(x,y,1)))
     }

     rv.avg <- function(x, period, align.by="seconds", align.period=1, cts=TRUE, makeReturns=FALSE, ...)
     {
          align.period = .getAlignPeriod(align.period, align.by)
          x<- .alignReturns(.convertData(x, cts=cts, makeReturns=makeReturns)$data, align.period)
          mean(.rv.subsample(x, period, ...))
     }

     rc.avg <- function(x, y,  period, align.by="seconds", align.period=1, cts=TRUE, makeReturns=FALSE, ...)
     {
          align.period = .getAlignPeriod(align.period, align.by)
          x<- .alignReturns(.convertData(x, cts=cts, makeReturns=makeReturns)$data, align.period)
          y<- .alignReturns(.convertData(y, cts=cts, makeReturns=makeReturns)$data, align.period)
          mean(.rc.subsample(x, y, period))
     }

     #########################################################################
     #
     # Naive estimators
     # See ABDL, BNS, etc 
     #
     #########################################################################
     rv.naive <- function(x, period, align.by = "seconds",align.period=1, cts=TRUE, makeReturns=FALSE, ...)
     {
     
          align.period = .getAlignPeriod(align.period, align.by)
          x<- .alignReturns(.convertData(x, cts=cts, makeReturns=makeReturns)$data, align.period)
          .C("rv", 
                    as.double(x), #a
                    as.double(x), #b
                    as.integer(length(x)), #na
                    as.integer(period), #period 
                    as.double(rep(0,as.integer(length(x)/period +1))), #tmp
                    as.double(rep(0,as.integer(length(x)/period +1))), #tmp
                    as.integer(length(x)/period), #tmpn
                    ans = double(1), 
                    COPY=c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,TRUE), 
                    PACKAGE="realized")$ans     
     }

     rc.naive <- function(x, y,  period,align.by = "seconds", align.period=1, cts=TRUE, makeReturns=FALSE, ...)
     {
          align.period = .getAlignPeriod(align.period, align.by)
          x<- .alignReturns(.convertData(x, cts=cts, makeReturns=makeReturns)$data, align.period)
          y<- .alignReturns(.convertData(y, cts=cts, makeReturns=makeReturns)$data, align.period)

          .C("rv", 
                    as.double(x), #a
                    as.double(y), #b
                    as.integer(length(x)), #na
                    as.integer(period), #period 
                    as.double(rep(0,as.integer(length(x)/period +1))), #tmp
                    as.double(rep(0,as.integer(length(x)/period +1))), #tmp
                    as.integer(length(x)/period), #tmpn
                    ans = double(1), 
                    COPY=c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,TRUE), 
                    PACKAGE="realized")$ans     
     }


     #########################################
     #
     # Global Call 
     #
     #######################################
     rRealized.variance <- function(x, y=NULL, type="naive", period = 1, lags = 1, cor = FALSE, rvargs = list(), cts=TRUE, makeReturns=FALSE)
     {
          warning("Depricated:  call rRealizedVariance")
          rRealizedVariance(x=x, y=y, type=type, period=period, lags=lags, cor=cor, rvargs=rvargs, cts=cts, makeReturns=makeReturns)
     }

     rRealisedVariance <- function(x, y=NULL, type="naive", period = 1, lags = 1, cor = FALSE, rvargs = list(), cts=TRUE, makeReturns=FALSE)
     {
          rRealizedVariance(x=x, y=y, type=type, period=period, lags=lags, cor=cor, rvargs=rvargs, cts=cts, makeReturns=makeReturns)
     }


     rRealizedVariance <- function(x, y=NULL, type="naive", period = 1, align.by="seconds", align.period = 1, cor = FALSE, rvargs = list(), cts=TRUE, makeReturns=FALSE, lags=NULL)
     {
          
         if(!is.null(lags)){
             warning("lags is deprecated, use the appropriate rv.* or rc.* paramter name in rvargs")
             rvargs$kernel.param = lags
         }   
         # Hack 
         if(length(which(names(rvargs)=="align.period")) > 0){
             warning("align.period in rvargs is deprecated, use as parameter")
             align.period=rvargs$align.period
             rvargs[which(names(rvargs)=="align.period") ] = NULL
             if(length(rvargs)==0)
                rvargs = list(bs=1)
             
         }
         
         
         if((n <- .getDataColNum(x)) > 1)
          {

               ans <- matrix(NA, nrow=n, ncol=n)
               for(i in 1:n)
               {
               for(j in 1:n)
               {
                    if(i == j)
                    {
                         if(cor)
                         {
                              ans[i,j] <- 1
                         }
                         else
                         {
                              ans[i,j] <- .realized.variance(x=.getDataCol(x,j), y=NULL, type=type, period=period, align.by=align.by,align.period=align.period, rvargs=rvargs, cts=cts, makeReturns=makeReturns)
                         }
                    }
                    else
                    {
                        if(j > i)
                        {
                         ans[i,j] <- .realized.variance(x=.getDataCol(x,i), y=.getDataCol(x,j), type=type, period=period,align.by=align.by,align.period=align.period, rvargs=rvargs, cor=cor, cts=cts, makeReturns=makeReturns)
                         ans[j,i]<-ans[i,j]
                        }     
                    }
               }
              }
              ans
          }
          else
          {
               ans <- .realized.variance(x=x, y=y, type=type, period = period, align.by=align.by, align.period=align.period, cor = cor, rvargs = rvargs, cts=cts, makeReturns=makeReturns)
          }
          ans
     }



     .realized.variance <- function(x, y=NULL, type="naive", period = 1, align.by="seconds", align.period = 1, cor = FALSE, rvargs = list(), cts=TRUE,makeReturns=FALSE)
     {   
     
          if(cor)
          {
               rvx <- do.call(paste("rv.", type, sep=""), c(rvargs,list(x=x, period=period, align.by=align.by, align.period=align.period, cts=cts, makeReturns=makeReturns)))
               rvy <- do.call(paste("rv.", type, sep=""), c(rvargs,list(x=y, period=period, align.by=align.by, align.period=align.period, cts=cts, makeReturns=makeReturns)))
               rcxy <- do.call(paste("rc.", type, sep=""),c(rvargs,list(x=x, y=y, period=period, align.by=align.by, align.period=align.period, cts=cts, makeReturns=makeReturns)))
               rcxy/(sqrt(rvx)*sqrt(rvy))
          }
          else
          {
               funct <- paste(ifelse(is.null(y), "rv.", "rc."), type, sep="")
               
               do.call(funct, c(rvargs, list(x=x, y=y, period=period, align.by=align.by, align.period=align.period, cts=cts, makeReturns=makeReturns)))     
          }
     }

     rc.zero <- function(x, y, period, align.by="seconds", align.period=1, cts=TRUE, makeReturns=FALSE, ...)
     {
          align.period = .getAlignPeriod(align.period, align.by)
          y<- .alignReturns(.convertData(y, cts=cts, makeReturns=makeReturns)$data, align.period)
          x<- .alignReturns(.convertData(x, cts=cts, makeReturns=makeReturns)$data, align.period)
          acy <- .accum.naive(x=y,y=y,period=period)
          acx <- .accum.naive(x=x,y=x,period=period)
         sum((acx*acy)==0)/length(acy)
     }

     rv.zero <- function(x, period, align.by="seconds",align.period=1, cts=TRUE, makeReturns=FALSE, ...)
     {
          align.period = .getAlignPeriod(align.period, align.by)   
          x<- .alignReturns(.convertData(x, cts=cts, makeReturns=makeReturns)$data, align.period)
          ac <- .accum.naive(x=x,y=x,period=period)
          sum((ac*ac)==0)/length(ac)
     }





     rCumSum <- function(x, period = 1, align.by="seconds", align.period=1, plotit=FALSE, type='l', cts = TRUE, makeReturns=FALSE)
     {
         align.period = .getAlignPeriod(align.period, align.by)   
         ans <- list(x = NULL, y = NULL)
         ans$x <- .alignIndices(1:length(.convertData(x, cts=cts, makeReturns=makeReturns)$data), align.period)
         ans$x <- .alignIndices(ans$x, period)

          x<- .alignReturns(.convertData(x, cts=cts, makeReturns=makeReturns)$data, align.period)
          x<- .alignReturns(.convertData(x, cts=cts, makeReturns=makeReturns)$data, period)

          ans$y <- cumsum(x)
          if(plotit)
          {
               plot(cumsum(x), xlab="Time", ylab="Cummulative Returns", type=type)
               return(NULL)
          }
          ans
     }


     rScatterReturns <- function(x,y, period, align.by="seconds", align.period=1,numbers=FALSE,xlim= NULL, ylim=NULL, plotit=TRUE, pch=NULL, cts=TRUE, makeReturns=FALSE, scale.size=0, col.change=FALSE,...)
     {
          align.period = .getAlignPeriod(align.period, align.by) 
          y<- .alignReturns(.convertData(y, cts=cts, makeReturns=makeReturns)$data, align.period)
          x<- .alignReturns(.convertData(x, cts=cts, makeReturns=makeReturns)$data, align.period)

          x<-.accum.naive(x, x, period)
          y<-.accum.naive(y, y, period)
          if(is.null(pch))
              pch=1

          it <- table(round(x,4),round(y,4))
          xs <- as.numeric(dimnames(it)[[1]])
          ys <- as.numeric(dimnames(it)[[2]])

          if(is.null(ylim))
              ylim=c(min(ys), max(ys))
          if(is.null(xlim))
              xlim=c(min(xs), max(xs))

          mat <- matrix(it, nrow=length(xs), ncol=length(ys))

          if(plotit)
          {
               plot(0,0, xlim=xlim, ylim=ylim , type='n',...)
               lines(c(0,0), c(-1,2), col="grey", lty=3, lwd=2)
               lines(c(-1,2), c(0,0), col="grey", lty=3, lwd=2)

               maxed <- max(mat)

               for(i in 1:length(xs))
               {
               for(j in 1:length(ys))
               {
               if(mat[i,j]!=0)
               {
                    if(col.change)
                       thecol <- round(runif(1)*100,0)
                    else
                       thecol = 1

                         if(numbers)
                         {

                              if(scale.size ==0)
                                        text(xs[i], ys[j],as.character(mat[i,j]), cex=.7, col=thecol)         
                    else
                         text(xs[i], ys[j], as.character(mat[i,j]), cex = (mat[i,j]/maxed) * scale.size, col=thecol)
                         }
                         else
                         {
                              if(scale.size ==0)
                              points(xs[i], ys[j], pch=pch, cex=.7, col=thecol)         
                    else
                              points(xs[i], ys[j], pch=pch, cex = (mat[i,j]/maxed) * scale.size, col=thecol)
                         }
               }

               }
               }
               return(NULL)

          }     
          mat
     }





     rc.hy <- function(x,y, period=1,align.by="seconds", align.period =1, cts = TRUE, makeReturns=FALSE, ...)
     {
          align.period = .getAlignPeriod(align.period, align.by)
          cdata <- .convertData(x, cts=cts, makeReturns=makeReturns)
          x <- cdata$data
          x.t <- cdata$milliseconds

          cdatay <- .convertData(y, cts=cts, makeReturns=makeReturns)
          y <- cdatay$data
          y.t <- cdatay$milliseconds


          errorCheck <- c(is.null(x.t),is.na(x.t), is.null(y.t), is.na(y.t))
          if(any(errorCheck))
              stop("ERROR: Time data is not in x or y.")


                         sum(     .C("pcovcc", 
                    as.double(x), #a
                    as.double(rep(0,length(x)/(period*align.period)+1)),
                    as.double(y), #b
                    as.double(x.t), #a
                    as.double(rep(0,length(x)/(period*align.period)+1)), #a
                    as.double(y.t), #b
                    as.integer(length(x)), #na
                    as.integer(length(x)/(period*align.period)),
                    as.integer(length(y)), #na
                    as.integer(period*align.period),
                    ans = double(length(x)/(period*align.period)+1), 
                    COPY=c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,TRUE), 
                    PACKAGE="realized")$ans)
     }


     #######################################################################
     #
     # Graphs 
     #
     #######################################################################

     #
     # Signature Plot
     #
     rSignature <- function(range, x, y=NULL, type="naive", cor = FALSE, rvargs = list(), align.by="seconds", align.period =1,xscale=1,  plotit=FALSE, cts=TRUE, makeReturns=FALSE, iteration.funct=NULL, iterations=NULL, lags=NULL)
     {

        if(!is.null(iteration.funct) || !is.null(iterations))
            stop("iterations no longer work, please call this function mutliple times and average the result.")

        if(!is.null(lags)){
             warning("lags is deprecated, use the appropriate rv.* or rc.* paramter name in rvargs")
             rvargs$kernel.param = lags
         }   
         # Hack 
         if(length(which(names(rvargs)=="align.period")) > 0){
             warning("align.period in rvargs is deprecated, use as parameter")
             align.period=rvargs$align.period
             rvargs[which(names(rvargs)=="align.period") ] = NULL
             if(length(rvargs)==0)
                rvargs = list(bs=1)
             
         }
               
               if(type=="kernel"){
               kernfun <- function(r, x, y, type, cor, rvargs){
                   rvargs$kernel.param = r
                   .realized.variance(x=x, y=y, type=type, cor = cor, rvargs = rvargs)
               }
                   ans <- list(x=range*xscale, 
                          y=sapply(range, kernfun,x=x,y=y,type=type, cor=cor, rvargs=rvargs),
                          xgrid=range,
                          type = type,
                          cor = cor,
                          cov = is.null(y),
                          cts= cts)
               }
               else{
                   ans <- list(x=range*xscale, 
		                             y=sapply(range, function(r, x, y, type, period, align.by, align.period, cor, rvargs){.realized.variance(x=x, y=y, type=type, period = r, cor = cor, rvargs = rvargs, align.by=align.by, align.period=align.period)},x=x,y=y,type=type, cor=cor, rvargs=rvargs, align.period=align.period, align.by=align.by),
		                             xgrid=range,
		                             type = type,
		                             cor = cor,
		                             cov = is.null(y),
		                             cts= cts)
               }
               
      
         if(plotit)
         {
          .rSignature.plot(ans)
         }
         ans
     }

     .rSignature.plot <- function(obj)
     {
         if(obj$cov && obj$cor)
         {
          ylab = "Realized Correlation"
         }
         else{
          if(obj$cov)
              ylab = "Realized Covariance"
          else
              ylab = "Realized Variance"
         }
         xlab = "Sampling Frequency"
         main = paste(ylab, ":", obj$type, sep="")
         plot(obj$x, obj$y, xlab=xlab, ylab=ylab, main=main)          
     }


     #
     # Helper functions
     #
     .getDataColNum <- function(x)
     {

          if("xts" %in% class(x))
         {
          return(dim(x)[[2]])     
         }

          if("realizedObject" %in% class(x))
         {
          return(dim(x)[[2]])     
         }
         if(is.null(version$language)) #splus
         {
               if("timeSeries" %in% class(x))
               {
                return(dim(x)[[2]])
               }
         }

         if("list" %in% class(x))
         {
          if(is.null(x$data))
          {
               return(NA)
          }
          else
          {
               if("matrix" %in% class(x) || "data.frame" %in% class(x))
               {
                    return(dim(x$data)[[2]])
               }
               else
               {
                    return(1)
               }
          }
         }
         else
         {     
          if("matrix" %in% class(x) || "data.frame" %in% class(x))
          {
               return(dim(x)[[2]])
          }
          else
          {
               return(1)
          }
         }
     }


     .getDataCol <- function(x,i)
     {
          if(is.null(x))
          {
               return(NULL)
          }

         
	 if("xts" %in% class(x))
	 {
	           return(x[,i])     
         } 
         if("realizedObject" %in% class(x))
         {
          return(x[,i])     
         }     
          if(is.null(version$language)) #splus
         {
               if("timeSeries" %in% class(x))
               {
                return(x[,i])
               }
         }

         if("list" %in% class(x))
         {
          if(is.null(x$data))
          {
               return(x)
          }
          else
          {
               if(class(x$data)=="matrix" || class(x$data)=="data.frame")
               {
                    return(x$data[,i])
               }
               else
               {
                    return(x$data)
               }
          }
         }
         else
         {     
          if("matrix" %in% class(x) || "data.frame" %in% class(x))
          {
               return(x[,i])
          }
          else
          {
               return(x)
          }
         }
     }

     rMarginal <- function(x, y=NULL, period, align.by="seconds", align.period=1, plotit=FALSE, cts=TRUE, makeReturns=FALSE)
     {
         align.period = .getAlignPeriod(align.period, align.by)   
         ans <- list(x = NULL, y = NULL)
         ans$x <- .alignIndices(1:length(x), align.period)
         ans$x <- .alignIndices(ans$x, period)

          if(is.null(y))
              y <- x

          x<- .alignReturns(.convertData(x, cts=cts, makeReturns=makeReturns)$data, align.period)
          y<- .alignReturns(.convertData(y, cts=cts, makeReturns=makeReturns)$data, align.period)


         ans$y <- .alignedAccum(x=x, y=y, period=period, cum=FALSE)

         if(plotit)
         {
          plot(ans, xlab="", ylab="Realized Marginals")
          return(NULL)
         }
         ans
     }

     rAccumulation <- function(x, period=1, y=NULL, align.by="seconds",align.period=1, plotit=FALSE, cts=TRUE, makeReturns=FALSE)
     {
         align.period = .getAlignPeriod(align.period, align.by)   
         ans <- list(x=NULL, y=NULL)
         ans$y <- cumsum(rMarginal(x=x, y=y, period=period, align.period=align.period, cts=cts, makeReturns=makeReturns)$y)
     #    ans$x <- .alignIndices(1:length(x), align.period)
      #   ans$x <- .alignIndices(ans$x, period)
         ans$x <- rCumSum(x=x, period=period, align.period=align.period, cts=cts, makeReturns=makeReturns)$x
         #ans$x <- ans$x[-length(ans$x)]
         if(plotit)
         {
          plot(ans, xlab="", ylab="Realized Accumulation")
          return(NULL)
         }
         ans
     }




     realizedObject <- function(x, cts = TRUE, makeReturns = FALSE, millisstart=NA, millisend=NA)
     {

         y <- .convertData(x, cts=cts, millisstart=millisstart, millisend=millisend, makeReturns=makeReturns)     
         class(y) <- "realizedObject"
         y$cts = cts
         y
     }

     print.realizedObject <- function(x, ...)
     {
         if(class(x$data) == "matrix")
         {
          n <- dim(x$data)[[1]]
          k <- dim(x$data)[[2]]
          display.n <- ifelse(n < 10, n, 10)
          display.obj <- cbind(x$data[1:display.n,], x$milliseconds[1:display.n])
          dims <- dimnames(x$data)[[2]]
          if(is.null(dims))
              dims <- paste("data", 1:k, sep="")

              dimnames(display.obj) <- list(rep("",display.n), c(dims, "milliseconds"))
         }
         else
         {
          n <- length(x$data)
          k <- 1
          display.n <- ifelse(n < 10, n, 10)
          display.obj <- matrix(c(x$data[1:display.n], x$milliseconds[1:display.n]), ncol=1+k)
          dimnames(display.obj) <- list(rep("",display.n), c("data", "milliseconds"))

         }
         cat("Realized Object: (length = ", n, ", cts=", x$cts,")\n")
          print(display.obj)
          cat("...")


     }


     "[.realizedObject" <- function(x, i=NULL,j=NULL, drop=TRUE)
     {
          ret <- x
          ret$milliseconds <- x$milliseconds[i]
          if(class(x$data) == "matrix")
         {
          if(is.null(i))
              i <- 1:(dim(x$data)[[1]]) 
          if(is.null(j))
              j <- 1:(dim(x$data)[[2]]) 
               ret$data <- x$data[i,j]
         }
         else
          ret$data <- x$data[i]
          ret
     }

     dim.realizedObject <- function(x)
     {
         if(class(x$data) == "matrix")
          return(dim(x$data))
         else
          return(c(length(x$data), 1))
     }
     merge.realizedObject <- function(x, y=NULL,...)
     {
          if(is.null(y))
          {
               return(x)
          }

         k <- nargs()
         inputs <- list(x,y, ...)
         inputs.class <- sapply(1:k, function(x, inputs){class(inputs[[x]])}, inputs)
          if(sum(inputs.class != "realizedObject"))
          {
               stop("merge.realizedObject takes object of type realizedObject only.")
          }
          inputs.len <- sapply(1:k, function(x, inputs){length(inputs[[x]]$milliseconds)}, inputs)
          if(sum(inputs.len != inputs.len[[1]]))
          {
               stop("Cannot merge objects with different timings")
          }
          if(k == 1)
          {
               return(x)
          }
          else
          {
               for(i in 2:k)
               {
                   x$data <- cbind(x$data,inputs[[i]]$data)
               }
          }
          x
     }

     #
     # Data Handling
     #
     .convertData <- function(x, cts = TRUE, millisstart=NA, millisend=NA, makeReturns=FALSE)
     {
          if(is.null(x))
          {
               return(NULL)
          }
          if("realizedObject" %in% class(x))
          {
               return(x)
          }
         if(is.null(version$language)) #splus
         {
               if("timeSeries" %in% class(x))
               {
               x <- x[!is.na(x[,1]),1]
               if(cts)
               {
                   return(ts2realized(x, millisstart=millisstart, millisend=millisend, make.returns=makeReturns)$cts)
               }
               else
               {
                   return(ts2realized(x, millisstart=millisstart, millisend=millisend, make.returns=makeReturns)$tts)
               }
               #list(milliseconds = positions(x)@.Data[[2]], data = matrix(seriesData(x), ncol=1))
               }
         }

          if("xts" %in% class(x))
          {
               xtmp <- x
               x <- list() 
               x$data <- as.numeric(xtmp[,1])

               x$milliseconds <- (as.POSIXlt(time(xtmp))$hour*60*60 + as.POSIXlt(time(xtmp))$min*60 + as.POSIXlt(time(xtmp))$sec )*1000
               if(is.na(millisstart))
               {
                   millisstart = x$milliseconds[[1]]
               }
               if(is.na(millisend))
               {
                   millisend = x$milliseconds[[length(x$milliseconds)]]
               }

               cat(paste("xts -> realizedObject [", as.character(time(xtmp[1])), " :: ", as.character(time(xtmp[length(x$milliseconds)])), "]", sep=""),"\n")
          }

          if(is.na(millisstart))
          {
              millisstart=34200000
          }
          if(is.na(millisend))
          {
              millisend=57600000
          }    
          if("list" %in% class(x))
          {
               if(sum(names(x) == c("tts", "cts")) == 2) #realized obj  
               {
                   if(cts)
                   {
                  return(x$cts)
               }
               else
               {
                   return(x$tts)
              }
           }
           if(sum(names(x) == c("data", "milliseconds")) == 2) 
           {
              if(makeReturns)
                   {                                           # only works on non cts prices
                    errcheck <- try(.getReturns(.sameTime(x$data, x$milliseconds)))
                    if(class(errcheck) != "Error")
                    {
                         x$data <- errcheck
                         x$milliseconds <- intersect(x$milliseconds,x$milliseconds)
                    }
                    else
                    {
                         warning("It appears that these are already returns.  Not creating returns")
                    }
              }          
                   else
                   {
               x$data <- .sameTime(x$data, x$milliseconds)
               x$milliseconds <- intersect(x$milliseconds,x$milliseconds)
                   }          
               if(cts)
               {
                   toret <- list(data=.toCts(x=x$data, millis=intersect(x$milliseconds,x$milliseconds), millisstart=millisstart, millisend=millisend),
                           milliseconds=(((millisstart/1000)+1):(millisend/1000))*1000)
                   return(toret)
               }
               else
               {
                   toret <- list(data=x$data, 
                           milliseconds=intersect(x$milliseconds,x$milliseconds))
                   return(toret)
                }
               }
          }


          if("timeSeries" %in% class(x))
         {
          stop("R timeSeries not implmented yet. Convert to realized object")
         }
         return(list(milliseconds = 1:dim(as.matrix(x))[[1]], data = as.matrix(x)))  # not an object, fake the milliseconds and return
     }

     plot.realizedObject <- function(x,y=NULL,...)
     {
         plot(x$milliseconds, x$data, xlab="Milliseconds", ylab="", main="")
         NULL
     }


     .getReturns <- function(x)
     {
             x <- as.numeric(x)
             n <- length(x)[[1]]
             return(log(x[2:n]) - log(x[1:(n-1)]))
     }




   #  timeDate <- function(x, format)
   #  {
   #     warning("This function is for SPLUS and does not work")
   #     x
   #  }



   #  .ts2millis <- function(x,...)
  #   {
#
  #        millis <- as.numeric(as.character((timeDate(as.character(x@positions), format="%H")))) * 60 * 60 * 1000 +
  #               as.numeric(as.character((timeDate(as.character(x@positions), format="%M")))) * 60 * 1000 +
  #               as.numeric(as.character((timeDate(as.character(x@positions), format="%S")))) * 1000 +
  #              as.numeric(as.character((timeDate(as.character(x@positions), format="%N"))))
  #        millis
  #   }


     ts2realized <- function(x, make.returns=TRUE,millisstart=34200000, millisend=57600000)
     {
           warning("SPLUS is no longer supported.")
     #     thedata <- data.sameTime(as.numeric(as.matrix(x@data)), .ts2millis(x))

     #    if(make.returns)
     #    {

     #          thedata <- .getReturns(thedata)

     #          tts <- list(data=as.numeric(thedata), milliseconds=intersect(.ts2millis(x),.ts2millis(x))[-1])
     #          cts <- list(data=.toCts(x=as.numeric(thedata), millis=intersect(.ts2millis(x),.ts2millis(x)), millisstart=millisstart, millisend=millisend),
     #               milliseconds=(((millisstart/1000)+1):(millisend/1000))*1000)
     #    }
     #    else
     #    {
     #          tts <- list(data=as.numeric(thedata), milliseconds=intersect(.ts2millis(x),.ts2millis(x)))
     #          cts <- list(data=.toCts(x=as.numeric(thedata), millis=intersect(.ts2millis(x),.ts2millis(x)), millisstart=millisstart, millisend=millisend),
     #               milliseconds=(((millisstart/1000)+1):(millisend/1000))*1000)


     #    }
     #     ans <- list(tts=tts, cts=cts)     
     #     ans
     }


     tsGetDay<-function(ts, dateString)
     {
           warning("SPLUS is no longer supported.")
   
   #     if(is(ts, "timeSeries"))
     #          pos = ts@positions
     #     else stop("ts must be a timeSeries object")
     #     pos@format = "%02m/%02d/%Y"
     #     poschar = as.character(pos)
     #     inds <- poschar == dateString
     #     ts[inds]
     }


     tsGetDayObject<-function(x, i, cts=TRUE, ...)
     {
     warning("SPLUS is no longer supported.")
     #     ts = x
    #      dateString=i
    #      if(cts)
    #     ts2realized(tsGetDay(ts, dateString))$cts$data
    #     else
    #     ts2realized(tsGetDay(ts, dateString))$tts$data    
     }


     #tsGetDayObject(msftt.ts[,6], "05/01/1997")
               #modeule(finmetrics)
     #          yoyo<-(tsGetDay(msftt.ts, "05/01/1997")[,6])
        #     theday <- tsGetDay(msftt.ts, "05/01/1997")[,6]
         #    yoyo <- data.sameTime(as.numeric(as.matrix(theday@data)), .ts2millis(theday))
          #   yoyo <- getReturns(yoyo)



     data.sameTime <- function(x, millis)
     {
          .sameTime(x=x,millis=millis)
     }

     .sameTime <- function(x, millis)
     {
            .C("sametime", 
                    as.double(x), #a
                    as.integer(length(x)), #na
                    as.integer(millis), #millis
                    ans = double(length(union(millis,millis))), #tts
                    COPY=c(FALSE,FALSE,FALSE,TRUE), 
                    PACKAGE="realized")$ans
     }




     data.toCts <- function(x, millis, millisstart=34200000, millisend=57600000)
     {
          .toCts(x=x, millis=millis, millisstart=millisstart, millisend=millisend)
     }

     .toCts <- function(x, millis, millisstart=34200000, millisend=57600000)
     {
            .C("tocts", 
                    as.double(x), #a
                    as.integer(length(x)),
                    as.integer(millis), #millis
              as.integer(millisstart),
              as.integer(millisend),
                    ans = double(((millisend-millisstart)/1000)), #cts
                    COPY=c(FALSE,FALSE,FALSE,FALSE,TRUE), 
                    PACKAGE="realized")$ans
     }

     data.toReturns <- function(x)
     {
         x <- as.numeric(x)   
         n <- length(x)
         log(x[2:n]) - log(x[1:(n-1)])
     }












