## Copyright (C) 2014 Dalle Molle Institute for Artificial Intelligence.
##
##    This program is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.
##
##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with this program.  If not, see <http://www.gnu.org/licenses/>.



#######################################################################################
##                              Functions                                            ##
#######################################################################################

repmat = function(X,m,n){
  
  ## R equivalent of repmat (matlab)
  
  mx = dim(X)[1]
  nx = dim(X)[2]
  matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T)
}
#####################


exp.survsamp <- function(n, lambda) {
  
  ## sample survival and censoring times both from the exponential distribution 
  ## with parameter lambda
  
  deaths <- matrix(rexp(n, rate = lambda),n,1)
  censor <- matrix(rexp(n, rate = lambda),n,1)
  status <- matrix((deaths<censor)*1,n,1)
  time <- deaths*status+censor*(1-status)
  temp <- data.frame(time,status)
}
#####################


expunif.survsamp <- function(n, lambda, UB) {
  
  ## sample survival times from the exponential distribution with parameter lambda
  ## and censoring time from the uniform distribution in the range in UB
  
  deaths <- matrix(rexp(n, rate = lambda),n,1)
  censor <- matrix(runif(n,UB[1],UB[2]),n,1)
  status <- matrix((deaths<censor)*1,n,1)
  time <- deaths*status+censor*(1-status)
  temp <- data.frame(time,status)
}
#####################


step.hazard.survsamp <- function(n, lvec, tvec, UB) {
  
  ## sample survival times with step hazard function and censoring time from the 
  ## uniform distribution in the range in UB. 
  ## lvec: vector of hazard values 
  ## tvec: vector of jump times 
  
  tvec <- c(0,tvec)
  dt <- diff(tvec)
  dhaz <- lvec[-length(lvec)]*dt
  chaz <- c(0,cumsum(dhaz))
  chaz.samp <- -log(matrix(runif(n,0,1),n,1))
  tbin <- as.numeric(cut(chaz.samp,chaz))
  tbin[is.na(tbin)]<-length(lvec)
  deaths <- (chaz.samp-chaz[tbin])/lvec[tbin]+tvec[tbin]
  censor <- matrix(runif(n,UB[1],UB[2]),n,1)
  status <- matrix((deaths<censor)*1,n,1)
  time <- deaths*status+censor*(1-status)
  temp <- data.frame(time,status)
}
#####################


maketest <- function(dataset,groups,group_label,s=0.25,display=TRUE,
                     alpha=0.05,nsamples=10000,smax=FALSE){
  
  ## Test the difference between two groups defined bay the values 'groups' of the 
  ## variable 'group_label' using the log-rank, the Peto-Peto and the IDP sumrank tests.
  
  colnames(dataset)[which(colnames(dataset)==group_label)]<-"lab"
  fdata <- Surv(time, status) ~ lab
  testLR <- survdiff(fdata, dataset, lab==groups[1]|lab==groups[2], rho = 0)
  testPP <- survdiff(fdata, dataset, lab==groups[1]|lab==groups[2], rho = 1)
  Xlab <- paste('lab',groups[1],sep='=')
  ix <- which(as.data.frame(testLR$n)[,1]==Xlab)
  zLR <- sign(testLR$obs[ix]-testLR$exp[ix])*sqrt(testLR$chisq)
  pLR <- 1-pnorm(zLR)
  #   pLR2 <- 2*pnorm(-sqrt(testLR$chisq))
  hLR <- 1*(pLR<alpha)
  zPP <- sign(testPP$obs[ix]-testPP$exp[ix])*sqrt(testPP$chisq)
  pPP <- 1-pnorm(zPP)
  #   pPP2 <- 2*pnorm(-sqrt(testPP$chisq))
  hPP <- 1*(pPP<alpha)
  out <- isurvdiff(fdata,dataset,groups=groups, alternative = 'greater',
                   display=FALSE, s=s, nsamples=nsamples)
  hIDP <- out$h
  
  if (display==TRUE) {
    out2 <- isurvfit(fdata,dataset,subset=(lab==groups[1]|lab==groups[2]),s=s)
  }
  
  temp <- list("h"=c(hLR,hPP,hIDP),"p"=c(pLR,pPP,out$prob[1],out$prob[2]), 
               "n"=testLR$n)
  
  if (smax==TRUE) {
    s.out <- isurvdiff.smax(fdata,dataset,groups, alternative = 'greater', 
                            nsamples=nsamples,smax=7)
    temp$s=s.out$s
  } 
  temp
}
#####################


maketestSTD <- function(dataset,groups,group_label,alpha=0.05){
  
  ## Test the difference between two groups defined bay the values 'groups' of the 
  ## variable 'group_label' using the log-rank and the Peto-Peto tests.
  
  colnames(dataset)[which(colnames(dataset)==group_label)]<-"lab"
  fdata <- Surv(time, status) ~ lab
  testLR <- survdiff(fdata, dataset, lab==groups[1]|lab==groups[2], rho = 0)
  testPP <- survdiff(fdata, dataset, lab==groups[1]|lab==groups[2], rho = 1)
  Xlab <- paste('lab',groups[1],sep='=')
  ix <- which(as.data.frame(testLR$n)[,1]==Xlab)
  zLR <- sign(testLR$obs[ix]-testLR$exp[ix])*sqrt(testLR$chisq)
  pLR <- 1-pnorm(zLR)
  pLR2 <- 2*pnorm(-sqrt(testLR$chisq))
  hLR <- 1*(pLR<alpha)
  zPP <- sign(testPP$obs[ix]-testPP$exp[ix])*sqrt(testPP$chisq)
  pPP <- 1-pnorm(zPP)
  pPP2 <- 2*pnorm(-sqrt(testPP$chisq))
  hPP <- 1*(pPP<alpha)
  temp <- list("h"=c(hLR,hPP))
  temp
}
####################