% File IDPSurvival/man/isurvdiff.smax.Rd
% Part of the IDPSurvival package for R (http://www.R-project.org)
% Copyright (C) 2014 Dalle Molle Institute for Artificial Intelligence
% Distributed under FDL 1.3 or later
% http://www.gnu.org/licenses/fdl-1.3.html

\name{isurvdiff.smax}
\alias{isurvdiff.smax}
\title{
Maximum values of s for which the IDP test returns a determinate decision
}
\description{
Search for the maximum values of parameter s for which the IDP test \code{isurvdiff(formula,...)}
issues a determinate decision. 
The function test values of s up to the parameter smax. If for smax the IDP test is still determinate,  
\code{isurvdiff.smax} returns list(smax,testout). If for s=0 the test is already indeterminate,  
\code{isurvdiff.smax} returns list(-1,testout), where testout is the last executed test.
}
\usage{
isurvdiff.smax(formula, ..., verbose=FALSE, accuracy=0.05, smax=12)
}
\arguments{
\item{formula}{a formula expression of the form \code{Surv(time, status) ~ predictor}. 
A single predictor is admitted.}
\item{verbose}{whether to display each value of s that is tried}
\item{accuracy}{to which precision s should be computed}
\item{smax}{to which maximum value s should be tried}
\item{...}{All arguments of \code{isurvdiff.smax} are passed to \code{isurvdiff} to perform the test.
Refer to the help of \code{isurvdiff} for more details about the arguments.}
}

\value{
A list with components:
\item{s}{
The maximum value of s for which the test returns a determinate decision (H=0 or H=1).
}
\item{test0}{
The value returned by \code{isurvdiff(formula,...)} for the last test performed.
Refer to the help of \code{isurvdiff} for more details.
}
}

\section{METHOD}{
This function implements the IDP sum-rank test describe in Mangili and others (2014).
}

\section{References}{
Benavoli, A., Mangili, F., Zaffalon, M. and Ruggeri, F. (2014). Imprecise Dirichlet process with application to the hypothesis test on the probability that X < Y.  \emph{ArXiv e-prints}, \url{http://adsabs.harvard.edu/abs/2014arXiv1402.2755B}. 

Mangili, F., Benavoli, A., Zaffalon, M. and de Campos, C. (2014). Imprecise Dirichlet Process for the estimate and comparison of survival functions with censored data. 
}

\seealso{
\code{\link{Surv}},
\code{\link{isurvdiff}}.  
}

\examples{
data(lung,package='survival')
lung <- lung[1:40,]	# reduced data set just to ensure that the
 					# example is very fast to run for building the package
test <-isurvdiff.smax(Surv(time,status)~sex,lung,groups=c(1,2), 
	 	alternative = 'two.sided', nsamples=1000) 
                    # better to use larger value of nsamples
					# this small value is to run it quickly
print(test$test0)
cat("Maximum s giving the same decision: ",test$s)

}
\keyword{survival}
\keyword{IDP}
