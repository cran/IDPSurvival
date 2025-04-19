% File IDPSurvival/man/isurvdiff.Rd
% Part of the IDPSurvival package for R (http://www.R-project.org)
% Copyright (C) 2014 Dalle Molle Institute for Artificial Intelligence
% Distributed under FDL 1.3 or later
% http://www.gnu.org/licenses/fdl-1.3.html

\name{isurvdiff}
\alias{isurvdiff}
\alias{print.isurvdiff}
\title{
Test Survival Curves Differences for two right censored data
}
\description{
Tests if there is a difference between two survival curves based on two samples (X and Y)
with right censored data. More precisely it test whether the probabiliy P(X<Y) is greater
than, lower than or equal to 1/2.
The prior near-ignorance Dirichlet Process (IDP) rank sum test is used.
It returns the result of the deicison.
H=1 indicates that the alternative hypothesis is true with posterior
probability greater than \code{level}. 
H=0 indicates the hypothesis is not true with posterior greater than \code{level},
H=2 indicates an indeterminate instance. This means that the decision depends on 
the choice of the prior. 
}
\usage{
isurvdiff(formula, data, groups=c(1,2), s=0.25,
          alternative = c("two.sided", "less", "greater"),
          exact=NULL, level = 0.95, display=TRUE, 
          nsamples=10000, rope=0, tmax=NULL)
}
\arguments{
\item{formula}{
a formula expression of the form \code{Surv(time, status) ~ predictor}. 
A single predictor is admitted.
}
\item{data}{
an optional data frame in which to interpret the variables occurring in the
formula.
}
\item{groups}{
a vector of two element indicating which value of the predictor 
represents groups 1 and 2. 
}
\item{s}{
sets the value of the prior strength s of the Dirichlet Process.
}
\item{alternative}{
define the direction of the test:
 "greater" --evaluates the hypothesis  P(X < Y)>1/2, i.e., returns H=1 
             if the lower probability of the hypothesis is larger than 
             \code{level}, H=0 if the upper probability is smaller than 
             \code{level} and H=2 if the lower and upper probabilities 
             encompass \code{level};
 "less" -- evaluates the hypothesis  P(Y < X)>1/2;
 "two.sided" -- performs a two-sided Bayesian test, i.e., returns H=1 if 
                 1/2 is not included between the left bound of the 
                 lower and the right bound of the upper \code{level} 
                 HPD credible intervals, H=0 if 1/2 is included in both 
                 the upper and lower credible intervals, H=2 otherwise.
}
\item{exact}{
computes the posterior probability if value is TRUE, 
or uses a normal approximation if value is FALSE.  
If you omit this argument, \code{isurvdiff} uses the exact method 
if at least one group has less than 100 samples
and the approximate one otherwise.
}
\item{level}{
sets the significance level alpha = 1-\code{level} of the test. 
}
\item{display}{
determines whether the posterior distributions of P(X<Y) have to 
be plotted (TRUE) or not (FALSE).
}
\item{nsamples}{
if \code{exact=TRUE}, sets the number of samples used in the Monte Carlo 
computation of the posterior distributions. For faster but less accurate
results, one can tune down this parameter. For more accurate, one might
increase it.
}
\item{rope}{
introduces a (symmetric) Region of Practical
Equivalence (ROPE) around 1/2, i.e., [1/2-value,1/2+value].
}
\item{tmax}{
whether to consider the difference in survival up to time tmax. 
NULL is the default and means without limit.
}
}

\value{
a list with components:

\item{h}{
The decision of the test:
H=0 -> accept the null hypothesis; 
H=1 -> rejects the null hypothesis; 
H=2 -> indeterminate (a robust decision cannot be made).
}
\item{prob}{
the probability of the alternatice hypotesis P(X<Y)>1/2 
if \code{alternative="greater"} or
P(Y<X)>1/2 if \code{alternative="less"}.
}
\item{Lower.Cred.Int}{
lower HPD credible interval. 
Confidence level defined by \code{level}.
}
\item{Upper.Cred.Int}{
upper HPD credible interval. Confidence level defined by \code{level}.
}
\item{alternative}{
the direction of the test \code{"greater","less"} or \code{"two.sided"}.
}
\item{strata}{
the number of subjects contained in each group.
}
\item{exact}{
logical variable saying if the exact posterior distributions
have been computed (\code{TRUE}) or the Gaussian approximation
has been used (\code{FALSE}).
}
}

\section{METHOD}{
This function implements the IDP sum-rank test describe in Mangili and others (2014).
}

\section{References}{
Benavoli, A., Mangili, F., Zaffalon, M. and Ruggeri, F. (2014). Imprecise Dirichlet process with application to the hypothesis test on the probability that X < Y.  \emph{ArXiv e-prints}, \url{https://ui.adsabs.harvard.edu/abs/2014arXiv1402.2755B/abstract}. 

Mangili, F., Benavoli, A., Zaffalon, M. and de Campos, C. (2014). Imprecise Dirichlet Process for the estimate and comparison of survival functions with censored data. 
}

\seealso{
\code{\link[survival]{Surv}}.  
}

\examples{
test <-isurvdiff(Surv(time,status)~sex,lung,groups=c(1,2), 
	 			 alternative = 'two.sided',s=0.5, nsamples=1000)
print(test)

data(Aids2)
fdata <- Surv(time, status) ~ T.categ
dataset <- Aids2
groups=c("blood","haem")
dataset["time"]<-dataset[4]-dataset[3]
dataset[5]<-as.numeric(unlist(dataset[5]))
test <-isurvdiff(fdata,dataset,groups=groups,
                 alternative = 'greater',s=0.5, nsamples=1000)
print(test)
}
\keyword{survival}
\keyword{IDP}
