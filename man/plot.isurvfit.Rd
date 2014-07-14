\name{plot.isurvfit}
\alias{plot.isurvfit}
\title{
  Plot method for \code{isurvfit} objects 
}
\usage{
\method{plot}{isurvfit}(x, se.fit=TRUE, ...)
}
\arguments{
  \item{x}{
    an object of class \code{isurvfit}, usually returned by the
    \code{isurvfit} function. 
  }
  \item{se.fit}{
    determines whether confidence intervals will be plotted.  
  }
  \item{...}{other arguments passed to the standard plot function}
}
\description{
A plot of survival curves is produced, one curve for each strata. 
}
  
\seealso{
\code{\link{isurvfit}}   
}

\examples{
leukemia.surv <- isurvfit(Surv(time, status) ~ x, data = aml, display=FALSE) 
plot(leukemia.surv) 
legend(100, .9, c("Maintenance", "No Maintenance"), 
       lty=c(1,1),lwd=c(2,1),col=c('black','red'),pch=c(1,2)) 
title("IDP Curves\nfor AML Maintenance Study") 
}
\keyword{survival}
\keyword{hplot}
