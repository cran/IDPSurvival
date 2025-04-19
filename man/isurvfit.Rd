% File IDPSurvival/man/isurvfit.Rd
% Part of the IDPSurvival package for R (http://www.R-project.org)
% Copyright (C) 2014 Dalle Molle Institute for Artificial Intelligence
% Distributed under FDL 1.3 or later
% http://www.gnu.org/licenses/fdl-1.3.html

\name{isurvfit}
\alias{isurvfit}
\alias{print.isurvfit}
\title{Create survival curves based on the IDP model}
\description{
This function creates survival curves from right censored data 
using the prior near-ignorance Dirichlet Process (IDP).
}
\usage{
isurvfit(formula, data, s=0.5, weights, subset, display=TRUE,
         conf.type=c('exact',  'approx', 'none'), nsamples=10000,
         conf.int= .95)
}
\arguments{
  \item{formula}{
    a formula object, which must have a 
    \code{Surv} object as the  
    response on the left of the \code{~} operator and, if desired, terms  
    separated by + operators on the right. 
    For a single survival curve the right hand side should be \code{~ 1}.
  }
  \item{data}{
    a data frame in which to interpret the variables named in the \code{formula}, 
    \code{subset} and \code{weights} arguments. 
  }
  \item{s}{
    sets the value of the prior strength s of the Dirichlet Process.
  }
  \item{weights}{
    the weights must be finite and nonnegative; it is strongly recommended that  
    they be strictly positive, since zero weights are ambiguous, compared 
    to use of the \code{subset} argument.
  }
  \item{subset}{
    expression saying that only a subset of the rows of the data 
    should be used in the fit. 
  }
  \item{display}{
    determines whether the survival curves have to 
    be plotted (TRUE) or not (FALSE).
  }
  \item{conf.type}{
    a variable saying how the credible interval shold be computed: 
    'exact': Monte-Carlo smapling from the exact distribution,
    'approx': Gaussian approximation,
    'none': no credible interval is computed.
  }
  \item{nsamples}{
    number pf samples used to approximate the credible intervals 
    if \code{conf.type='exact'}.
    }
  \item{conf.int}{
    confidence level of the credible interval.
  }
}
\value{
an object of class \code{"isurvfit"}.  

See \code{isurvfit.object} for 
details. Methods defined for survfit objects are  
\code{print} and \code{plot}. 
}
\details{
The estimates are obtained using the IDP estimator by Mangili and others (2014)
based on the prior near-ignorance Dirichlet Process model 
by Benavoli and others (2014).
}

\section{References}{
Benavoli, A., Mangili, F., Zaffalon, M. and Ruggeri, F. (2014). Imprecise Dirichlet process with application to the hypothesis test on the probability that X < Y.  \emph{ArXiv e-prints}, \url{https://ui.adsabs.harvard.edu/abs/2014arXiv1402.2755B/abstract}.

Mangili, F., Benavoli, A., Zaffalon, M. and de Campos, C. (2014). Imprecise Dirichlet Process for the estimate and comparison of survival functions with censored data. 

}

\seealso{
\code{\link{isurvfit.object}},
\code{\link{plot.isurvfit}}, 
\code{\link[survival]{Surv}}.  
}

\examples{
data(aml)
fit <- isurvfit(Surv(time, cens) ~ 1, data=aml, display=TRUE, nsamples=1000) 
legend('topright', c("Lower expectation", 
          "Upper expectation","confidence intervals"), lty=c(1,1,2),lwd=c(1,2,1)) 
title("IDP survival curve (s=0.5) \nAcute Myelogenous Leukemia dataset")

data(Aids2)
dataset <- Aids2
dataset["time"]<-dataset[4]-dataset[3]
dataset[5]<-as.numeric(unlist(dataset[5]))
fit <- isurvfit(Surv(time, status) ~ T.categ, dataset,s=1,
	            subset=(!is.na(match(T.categ, c('blood','haem','het')))),
                nsamples=1000,conf.type='none')
legend('topright',c("Heterosexual contact","Hemophilia","Blood"),
            title="Transmission category:",lty=c(1,1,1),col=c(1,2,3),pch=c(1,2,3))
title("IDP survival curve (s=1) \nAids dataset")
print(fit)

leukemia.surv <- isurvfit(Surv(time, cens) ~ group, data = aml, display=FALSE) 
plot(leukemia.surv) 
legend(100, .9, c("Maintenance", "No Maintenance"), lty=c(1,1),lwd=c(2,1),
       col=c('black','red'),pch=c(1,2)) 
title("IDP Curves\nfor AML Maintenance Study") 
}

\keyword{survival}
\keyword{IDP}
