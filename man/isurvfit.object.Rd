\name{isurvfit.object}
\alias{isurvfit.object}
\title{
  IDP Survival Curve Object 
}
\description{
This class of objects is returned by the \code{isurvfit} class of functions
to represent a IDP survival curve. 

Objects of this class have methods for the functions \code{print},
\code{plot}.
}
\arguments{
\item{n}{
total number of subjects in each curve.
}
\item{time}{
the time points at which the curve has a step. 
}
\item{n.risk}{
the number of subjects at risk at t. 
}
\item{n.event}{
the number of events that occur at time t. 
}
\item{n.censor}{
the number of subjects who exit the risk set,
without an event,  at time t. 
(This number can be computed from the successive
values of the number at risk).
}
\item{survUP}{
the estimate of upper expectation of the survival 
probability at time t+0. 
}
\item{survLOW}{
the estimate of lower expectation of the survival 
probability at time t+0. 
}
\item{survLOW0}{
the estimate of lower expectation of the survival 
probability at time t=0. The upper expectation is always 1.
}
\item{std.up}{
the standard deviation of the upper distribution
of the survival probability. 
}
\item{std.low}{
the standard deviation of the lower distribution
of the survival probability. 
}
\item{upper}{
upper confidence limit for the survival curve. 
}
\item{lower}{
lower confidence limit for the survival curve. 
}
\item{lower0}{
lower confidence limit for the survival curve 
at t=0. The upper is always 1. 
}
\item{conf.type}{
the approximation used to compute the confidence limits. 
}
\item{conf.int}{
the level of the confidence limits, e.g. 90 or 95\%. 
}
\item{strata}{
number of elements of the \code{time} vector 
corresponding to each curve.  
The names of the elements are labels for the curves. 
}
\item{call}{
an image of the call that produced the object. 
}
}
\seealso{
\code{\link{isurvfit}}   
\code{\link{plot.isurvfit}}   
}
\keyword{survival}
\keyword{IDP}
