# IDPSurvival package for R (http://www.R-project.org)
# Copyright (C) 2014 Francesca Mangili, Alessio Benavoli, Marco Zaffalon, 
#                    Cassio de Campos.
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

#######
####### Code to generate the results in the paper:
####### Reliable survival analysis based on the Dirichlet Process
#######


library(IDPSurvival)

## IF YOU DONT HAVE THE LIBRARY, PLEASE DOWNLOAD IT FROM
## http://cran.r-project.org/web/packages/IDPSurvival/index.html
## AND INSTALL IT USING:
## R CMD install IDPSurvival.tar.gz
## OR INSIDE R:
## install.packages("IDPSurvival.tar.gz",repos=NULL,type="source")

library(survival)

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
source('functions.r')

# Setting initial seed to have the same results that we obtained.
set.seed(4321)

## This file should run in about 3 to 4 hours in a modern computer (as of April-2015).
initial.time <- proc.time()

#######################################################################################
##                              Section 3 - Simulations                              ##
#######################################################################################

## Figure 1 
## IDP (s=0.25) and Kaplan-Meier survival curves for a sample of both survival and 
## censoring times drawn from the exponential distribution with parameter lambda=5.
################################################

cat('Generating results for Figure 1... \n')

lambda <- 5
s <- 0.25
alpha <- 0.05

## Figure 1 (left): n=25
n <- 25
samp20 <- exp.survsamp(n,lambda)
iout <- isurvfit(Surv(time,status)~1, samp20, conf.int= 1-alpha, s=s, display=TRUE)
out <- survfit(Surv(time,status)~1, samp20, conf.int= 1-alpha,conf.type="log-log")
lines(out, col='red',mark=3)
t<-seq(from = 0., to = 0.6, by = 0.01)
lines(t, exp(-lambda*t), col='green')
legend('bottomleft',c("IDP upper","IDP lower", "IDP cred. int.",
                      "KM","KM conf. int.",expression(exp(-lambda*t))),
       lty=c(1,1,2,1,2,1),lwd=c(2,1,1,1,1,1), 
       col=c('black','black','black','red','red','green'),
       pch=c('o','o','o','+','.','.'), text.width=0.04,inset=0,bty="n")

## Figure 1 (right): n=100
n <- 100
samp20 <- exp.survsamp(n,lambda)
iout <- isurvfit(Surv(time,status)~1, samp20, conf.int= 1-alpha, s=s, display=TRUE)
out <- survfit(Surv(time,status)~1, samp20, conf.int= 1-alpha,conf.type="log-log")
lines(out, col='red',mark=3)
t<-seq(from = 0., to = 0.6, by = 0.01)
lines(t, exp(-lambda*t), col='green')
legend('bottomleft',c("IDP upper","IDP lower", "IDP cred. int.",
                      "KM","KM conf. int.",expression(exp(-lambda*t))),
       lty=c(1,1,2,1,2,1),lwd=c(2,1,1,1,1,1), 
       col=c('black','black','black','red','red','green'),
       pch=c('o','o','o','+','.','.'), text.width=0.04,inset=0,bty="n")
##############################


## Figure 2 and Table 1
## Coverage of the Kaplan-Meier estimator with plain and log-log confidence types 
## and of the IDP estimator with s = 0.25. Evaluated for n=25.
## Figure: s = 0.25
## Table: s= 0.01, 0.25, 0.5, 1 and 2. 
################################################

cat('Generating results for Figure 2 and Table 1... \n')

runs <- 1000
n.eval <- 10
step <- 0.05
sv <- c(0.01,0.25,0.5)
nv <- c(25, 100)
t.eval <-matrix(seq(from = 0.025, to = step*n.eval, by = step),n.eval,1)
ref <- t(exp(-lambda*t.eval))
cov <- vector("list", 2)
sumall <- vector("list", 2)
nall <- vector("list", 2)
sc<-0
for (ni in 1:2) {
  n<-nv[ni]
  cov[[ni]] <- matrix(numeric(0),5,0) 
  sumall[[ni]] <- rep(0,5)
  nall[[ni]] <-0
  t.mat <- repmat(t.eval,1,n)
  error <- matrix(numeric(0),n.eval,0)
  llerror <- matrix(numeric(0),n.eval,0)
  ierror <- list(llerror,llerror,llerror)
  for ( k in 1:runs) {
    samp <- exp.survsamp(n,lambda)
    
    # Estimate coverage for the Kaplan-Meier "plain" confidence interval
    out <- survfit(Surv(time,status)~1, samp, conf.int= 1-alpha,conf.type="log")
    time <- matrix(out$time,1,n)
    mat.time <- repmat(time,n.eval,1)
    time.ref <- apply((mat.time<t.mat),1,sum)+1
    # Set estimates before the first observed death to NA (not considered in the 
    # coverage estimation)
    n.death <- cumsum(out$n.event)
    out$lower[which(n.death==0)] <- NA
    out$upper[which(n.death==0)] <- NA
    # Add estimate in 0
    lower <- c(NA,out$lower)
    upper <-c(NA,out$upper)
    # 
    kerr <- matrix((lower[time.ref]<ref)*(upper[time.ref]>ref),n.eval,1)
    error <- cbind(error, kerr)
    
    # Estimate coverage for the the Kaplan-Meier "log-log" confidence interval
    out <- survfit(Surv(time,status)~1, samp, conf.int= 1-alpha,conf.type="log-log")
    time<-matrix(out$time,1,n)
    mat.time <- repmat(time,n.eval,1)
    time.ref <- apply((mat.time<t.mat),1,sum)+1
    n.death <- cumsum(out$n.event)
    # Set estimates before the first observed death to NA (not considered in the 
    # coverage estimation)
    out$lower[which(n.death==0)] <- NA
    out$upper[which(n.death==0)] <- NA
    # Add estimate in 0
    lower <- c(NA,out$lower)
    upper <-c(NA,out$upper)
    kerr <- matrix((lower[time.ref]<ref)*(upper[time.ref]>ref),n.eval,1)
    llerror <- cbind(llerror, kerr)
    
    # Estimate coverage for the IDP credible interval
    
    for (kk in 1:3){
      iout <- isurvfit(Surv(time,status)~1, samp, conf.int= 1-alpha, s=sv[kk], display=FALSE,nsamples=1000)
      ilower <- c(iout$lower0,iout$lower)
      iupper <-c(1, iout$upper)
      ikerr <- matrix((ilower[time.ref]<ref)*(iupper[time.ref]>ref),n.eval,1)
      ierror[[kk]] <- cbind(ierror[[kk]], ikerr)
    }
  }
 
   
  for (i in 1:n.eval) {
    # Remove NA (i.e., samples for which t.eval(i) is too small or too large)
    if (sum(is.na(error[i,]))>0) {
      i.err <- cbind(error[i,-which(is.na(error[i,]))],
                     llerror[i,-which(is.na(error[i,]))],
                     ierror[[1]][i,-which(is.na(error[i,]))],
                     ierror[[2]][i,-which(is.na(error[i,]))],
                     ierror[[3]][i,-which(is.na(error[i,]))]
      )
    } else {
      i.err <- cbind(error[i,], llerror[i,], ierror[[1]][i,],
                     ierror[[2]][i,],ierror[[3]][i,])
    }
    # Store results
    cov[[ni]] <- cbind(cov[[ni]],apply(i.err,2,mean))
    sumall[[ni]] <- sumall[[ni]]+apply(i.err,2,sum)
    nall[[ni]] <- nall[[ni]]+dim(i.err)[1]
  }
}

# Plot Figure 2 (left) : n = 25
plot(t.eval, cov[[1]][4,], xlab="Time", ylab="Coverage", type = "p",pch='x',ylim=c(0.5,1))
points(t.eval, cov[[1]][1,], type = "p")
points(t.eval, cov[[1]][2,], type = "p",pch='*')
lines(c(0,0.5),c(1-alpha,1-alpha))
legend('bottomleft',c("IDP","Kaplan-Meier plain","Kaplan Meier log-log"),pch=c('x','o','*'), bty="n")

# Plot Figure 2 (right) : n = 100
plot(t.eval, cov[[2]][4,], xlab="Time", ylab="Coverage", type = "p",pch='x',ylim=c(0.5,1))
points(t.eval, cov[[2]][1,], type = "p")
points(t.eval, cov[[2]][2,], type = "p",pch='*')
lines(c(0,0.5),c(1-alpha,1-alpha))
legend('bottomleft',c("IDP","Kaplan-Meier plain","Kaplan Meier log-log"),pch=c('x','o','*'), bty="n")

# Build Table 1
label.row <- c( "KM (plain)", "KM (log-log)", "IDP(s=0.01)",
        "IDP(s=0.25)", "IDP(s=0.5)") 
n25<-sumall[[1]]/nall[[1]]
n100<-sumall[[2]]/nall[[2]]
table1 = data.frame(n25,n100,row.names=label.row)

cat("---------------------------------------------------- \n")
cat("Table 1: Coverage of pointwise KM's confidence intervals and IDP's credible intervals. \n")
print(table1)
cat("---------------------------------------------------- \n")
##############################


## Table 2
## Power (or type I error) of the log-rank, Peto-Peto and IDP tests (with s = 0.25) 
## evaluated for the scenarios SH (same hazard) and PH (proportional hazard). 
################################################

cat('Generating results for Table 2... \n')

runs <- 1000
s <- 0.25
alpha <- 0.05
n <- 30
UB <- c(0,2)

# Same Hazard (SH) scenario
cat('Same Hazard (SH) scenario... \n')
lambda1 <- 1
lambda2 <- 1
Htest <- matrix(numeric(0),runs,3)
for (k in 1:runs) {
  Z1 <- expunif.survsamp(n, lambda1, UB)
  Z1$lab <- rep(0,n)
  Z2 <- expunif.survsamp(n, lambda2, UB)
  Z2$lab <- rep(1,n)
  dataset <- rbind(Z1,Z2)
  temp <-  maketest(dataset,groups=c(0,1),group_label="lab",s=s,display=FALSE,alpha=alpha)
  Htest[k,] <- temp$h
}
detdec <- matrix(Htest[which(Htest[,3]!=2),],sum(Htest[,3]!=2),3)
indetdec <- matrix(Htest[which(Htest[,3]==2),1:2],sum(Htest[,3]==2),2)
errDet <- apply(detdec,2,mean)
errIndet<- c(apply(indetdec,2,mean),mean(1*(Htest[,3]==2)))

# Proportional Hazard (PH) scenario
cat('Proportional Hazard (PH) scenario... \n')
lambda1 <- 2
lambda2 <- 1
Htest <- matrix(numeric(0),runs,3)
for (k in 1:runs) {
  Z1 <- expunif.survsamp(n, lambda1, UB)
  Z1$lab <- rep(0,n)
  Z2 <- expunif.survsamp(n, lambda2, UB)
  Z2$lab <- rep(1,n)
  dataset <- rbind(Z1,Z2)
  temp <-  maketest(dataset,groups=c(0,1),group_label="lab",s=s,display=FALSE,alpha=alpha)
  Htest[k,] <- temp$h
}

detdec <- matrix(Htest[which(Htest[,3]!=2),],sum(Htest[,3]!=2),3)
indetdec <- matrix(Htest[which(Htest[,3]==2),1:2],sum(Htest[,3]==2),2)
errDet <- cbind(errDet,apply(detdec,2,mean))
errIndet<- cbind(errIndet, c(apply(indetdec,2,mean),mean(1*(Htest[,3]==2))))


## Early Hazard Difference (EHD) scenario
cat('Early Hazard Difference (EHD) scenario... \n')
tvec <- c(0.4, 0.6)
lvec1 <- c(3, 0.75, 1)
lvec2 <- c(0.75, 3, 1)
Htest <- matrix(numeric(0),runs,3)
for (k in 1:runs) {
  Z1 <- step.hazard.survsamp(n, lvec1, tvec, UB) 
  Z1$lab <- rep(0,n)
  Z2 <- step.hazard.survsamp(n, lvec2, tvec, UB) 
  Z2$lab <- rep(1,n)
  dataset <- rbind(Z1,Z2)
  temp <-  maketest(dataset,groups=c(0,1),group_label="lab",s=s,display=FALSE,alpha=alpha)
  Htest[k,] <- temp$h
}

detdec <- matrix(Htest[which(Htest[,3]!=2),],sum(Htest[,3]!=2),3)
indetdec <- matrix(Htest[which(Htest[,3]==2),1:2],sum(Htest[,3]==2),2)
errDet <- cbind(errDet,apply(detdec,2,mean))
errIndet<- cbind(errIndet, c(apply(indetdec,2,mean),mean(1*(Htest[,3]==2))))


## Late Hazard Difference (LHD) scenario
cat('Late Hazard Difference (LHD) scenario... \n')
tvec <- c(0.4)
lvec1 <- c(2, 4)
lvec2 <- c(2, 0.4)
Htest <- matrix(numeric(0),runs,3)
for (k in 1:runs) {
  Z1 <- step.hazard.survsamp(n, lvec1, tvec, UB) 
  Z1$lab <- rep(0,n)
  Z2 <- step.hazard.survsamp(n, lvec2, tvec, UB) 
  Z2$lab <- rep(1,n)
  dataset <- rbind(Z1,Z2)
  temp <-  maketest(dataset,groups=c(0,1),group_label="lab",s=s,display=FALSE,alpha=alpha)
  Htest[k,] <- temp$h
}

detdec <- matrix(Htest[which(Htest[,3]!=2),],sum(Htest[,3]!=2),3)
indetdec <- matrix(Htest[which(Htest[,3]==2),1:2],sum(Htest[,3]==2),2)
errDet <- cbind(errDet,apply(detdec,2,mean))
errIndet<- cbind(errIndet, c(apply(indetdec,2,mean),mean(1*(Htest[,3]==2))))

# Build table 2
Log.rank.det <- errDet[1,]
Log.rank.indet <- errIndet[1,]
Peto.Peto.det <- errDet[2,]
Peto.Peto.indet <- errIndet[2,]
IDP.det <- errDet[3,]
IDP.indet <- errIndet[3,]
label.row <- c( "SH", "PH", "EHD","LHD")  
table2 = data.frame(Log.rank.det, Log.rank.indet, Peto.Peto.det, Peto.Peto.indet,
                    IDP.det,IDP.indet, row.names=label.row)

cat("---------------------------------------------------- \n")
cat("Table 2: Fraction of H1 decisions for the different tests in the cases where the IDP 
is determinate or indeterminate. Last column reports the fraction of indeterminate IDP outcomes. \n")
print(table2)
cat("---------------------------------------------------- \n")
##############################


## Table 3
## Type I error of of the log-rank, Peto-Peto (alpha=0.01 and 0.05) and IDP 
## (s=0.25 and alpha=0.05) tests for scenarios A and B.
################################################

cat('Generating results for Table 3... \n')

runs <- 1000
s <- 0.25
UB <- c(0,2)

## Case A scenario
cat('Case A scenario... \n')
# ## Case A scenario
n1 <- 50
n2 <- 25
tvec <- c(0.4)
lvec1 <- c(0.75, 5)
lvec2 <- c(2, 0.25)
Htest <- matrix(numeric(0),runs,5)
for (k in 1:runs) {
  Z1 <- step.hazard.survsamp(n1, lvec1, tvec, UB) 
  Z1$lab <- rep(0,n1)
  Z2 <- step.hazard.survsamp(n2, lvec2, tvec, UB) 
  Z2$lab <- rep(1,n2)
  dataset <- rbind(Z1,Z2)
  temp05 <- maketest(dataset,groups=c(0,1),group_label="lab",s=s,display=FALSE,alpha=0.05)
  temp01 <- maketestSTD(dataset,groups=c(0,1),group_label="lab",alpha=0.01)
  Htest[k,] <-c(temp01$h,temp05$h)
}

detdec <- matrix(Htest[which(Htest[,5]!=2),],sum(Htest[,5]!=2),5)
errDet <- apply(detdec,2,mean)
if (sum(Htest[,5]==2)>0) {
  indetdec <- matrix(Htest[which(Htest[,5]==2),1:4],sum(Htest[,5]==2),4)
  errIndet <-  c(apply(indetdec,2,mean),mean(1*(Htest[,5]==2)))
} else {errIndet <- c(rep(NA,4),0)}


## Case B scenario
cat('Case B scenario... \n')
n1 <- 25
n2 <- 50
tvec <- c(0.4)
lvec1 <- c(2, 0.2)
lvec2 <- c(0.75, 5)
Htest <- matrix(numeric(0),runs,5)
for (k in 1:runs) {
  Z1 <- step.hazard.survsamp(n1, lvec1, tvec, UB) 
  Z1$lab <- rep(0,n1)
  Z2 <- step.hazard.survsamp(n2, lvec2, tvec, UB) 
  Z2$lab <- rep(1,n2)
  dataset <- rbind(Z1,Z2)
  temp05 <- maketest(dataset,groups=c(0,1),group_label="lab",s=s,display=FALSE,alpha=0.05)
  temp01 <- maketestSTD(dataset,groups=c(0,1),group_label="lab",alpha=0.01)
  Htest[k,] <- c(temp01$h,temp05$h)
}

detdec <- matrix(Htest[which(Htest[,5]!=2),],sum(Htest[,5]!=2),5)
errDet <- cbind(errDet,apply(detdec,2,mean))
if (sum(Htest[,5]==2)>0) {
  indetdec <- matrix(Htest[which(Htest[,5]==2),1:4],sum(Htest[,5]==2),4)
  errIndet <-  cbind(errIndet,c(apply(indetdec,2,mean),mean(1*(Htest[,5]==2))))
} else {errIndet <-  cbind(errIndet,c(rep(NA,4),0))}

# Build table 3
LR.det.05 <- errDet[3,]
LR.indet.05 <- errIndet[3,]
LR.det.01 <- errDet[1,]
LR.indet.01 <- errIndet[1,]
PP.det.05 <- errDet[4,]
PP.indet.05 <- errIndet[4,]
PP.det.01 <- errDet[2,]
PP.indet.01 <- errIndet[2,]
IDP.det <- errDet[5,]
IDP.indet <- errIndet[5,]
label.row <- c( "case A","case B")  
table3 = data.frame(LR.det.05,LR.indet.05,LR.det.01,LR.indet.01,
                    PP.det.05,PP.indet.05,PP.det.01,PP.indet.01,
                    IDP.det,IDP.indet, row.names=label.row)

cat("----------------------------------------------------\n ")
cat("Table 3: Type-I error of different tests over the determinate and indeterminate cases.
Last column fraction of IDP indeterminate outcomes. \n")
print(table3)
cat("----------------------------------------------------\n ")
##############################

#######################################################################################
##                Section 4 - Australian AIDS survival data                          ##
#######################################################################################
 
## Figure 4 
##Survival curves estimated with IDP (s = 0:25) for the groups in Table 4.
################################################


## Table 4
## p-values (Log-rank and Peto-Peto tests), posterior probabilities (IDP test with s=0.25) 
## and maximum s giving a determinate decision for the AIDS dataset. 
################################################

cat('Generating results for Figure and Table 4... \n')

s<-0.25

data(Aids2,package='MASS')
dataset <- Aids2
dataset["time"]<-dataset[4]-dataset[3]
dataset["agg"] <- 1*(dataset["T.categ"]=="hs")+1*(dataset["T.categ"]=="het")+
  2*(dataset["T.categ"]=="hsid")+2*(dataset["T.categ"]=="id")  
dataset[5]<-as.numeric(unlist(dataset[5]))

cat("Male < Female\n")
t1 <- maketest(dataset,c("M","F"),"sex",s=s,smax=TRUE)
legend('bottomleft',c("Female (X\')","Male (X)"),lty=c(1,1),
       lwd=c(2,2),col=c('black','red'),pch=c(1,3),bty="n")
table4 <-c(t1$n[2],t1$n[1],t1$p,t1$s) 


cat("NSW < VIC \n")
t1 <- maketest(dataset,c("NSW","VIC"),"state",s=s,smax=TRUE)
legend('bottomleft',c("New South Wales (X)","Victoria (X\')"),
       lty=c(1,1),lwd=c(2,2),col=c('black','red'),pch=c(1,3),bty="n")
table4 <-rbind(table4, c( t1$n,t1$p,t1$s))

cat("QLD < NSW \n")
t1 <- maketest(dataset,c("QLD","NSW"),"state",s=s,smax=TRUE)
legend('bottomleft',c("New South Wales (X\')","Queensland (X)"),
       lty=c(1,1),lwd=c(2,2),col=c('black','red'),pch=c(1,3),bty="n")
table4 <-rbind(table4, c(t1$n[2],t1$n[1],t1$p,t1$s))

cat("no drug user < drug user \n")
t1 <- maketest(dataset,c(1,2),"agg",s=s,smax=TRUE)
legend('bottomleft',c("No drug user (X)","Drug user (X\')"),lty=c(1,1),
       lwd=c(2,2),col=c('black','red'),pch=c(1,3),bty="n")
table4 <-rbind(table4, c(t1$n,t1$p,t1$s))

cat("Blood < Haemophilia \n")
t1 <- maketest(dataset,c("blood","haem"),"T.categ",s=s,smax=TRUE)
legend('bottomleft',c("Hemophilia (X\')","Blood (X)"),lty=c(1,1),
       lwd=c(2,2),col=c('black','red'),pch=c(1,3),bty="n")
table4 <-rbind(table4, c(t1$n[2],t1$n[1],t1$p,t1$s))

label.row <- c("M vs F","NSW vs VIC","QLD vs NSW","no drus vs drug","blood vs haem")
table4<-data.frame(table4,row.names=label.row)
colnames(table4)<-c("n1","n2","LogRank","PetoPeto","IDPlower","IDPupper","smax")

cat("---------------------------------------------------- \n ")
cat("Table 4: p-values (Log-rank and Peto-Peto tests) and posterior probabilities (IDP test) 
for the AIDS dataset.\n") 
print(table4)
cat("----------------------------------------------------\n ")
# ##############################


## finally, just print the total time to run it.


run.time <- (proc.time()[3]-initial.time[3])/60
cat(run.time,"min")



