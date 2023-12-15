### clear all
rm(list=ls()) 

### Import Packages 
library(qvalue)
library(locfdr)
library(mixfdr)

### Functions 
source('ash_functions.R')

### Effects of UA

## Simulation Settings
set.seed(123)   # random seed 
J=10000         # dimension
beta=rnorm(J)   # true effects
est_s=rep(1,J)  # estimated standard error
est_beta=rnorm(J,mean=beta,sd=est_s) # estimated effects
lfdr=rep(0,J)   # local false sign rate

sigmas=get_sigmas(est_beta,est_s)  # sigmas from grid
K=length(sigmas)                   # number of components
EM_result=EM(J,K,est_beta,est_s)   # parameter estimation by EM algorithm
est_pi=EM_result$pi
est_w=EM_result$w
for(j in 1:J)
  lfdr[j]=get_lfdr(j,est_beta,est_s,est_w)
samples=get_samples(est_beta,est_s)

## plots (Figure 1 and Figure 2)
lfdr_l=locfdr(samples$z,nulltype=0)     # locfdr
lfdr_m=mixFdr(samples$z,theonull=TRUE)  # mixfdr
dev.off()
# ash results
hist(samples$p,breaks=seq(0,1,length=100),col='lightblue',main='ash',xlab='p-value')
hist(samples$pp,breaks=seq(0,1,length=100),add=T,col='blue')
hist(samples$z,breaks=seq(-6,6,length=100),col='blue',main='ash',xlab='z-score')
hist(samples$zz,breaks=seq(-6,6,length=100),add=T,col='lightblue')


### Different Measurement Precision
set.seed(123)   # random seed 
J=10000         # dimension
beta=(runif(J)>0.5)*rnorm(J)      # mixture of null and normal
est_s=c(rep(1,J/2),rep(10,J/2))   # different measurement precision
est_beta=rnorm(J,mean=beta,sd=est_s)  # estimated effects

z=est_beta/est_s                # z-scores
p=2*pmin(pnorm(z),1-pnorm(z))   # p-values

## plots (Figure 3)
par(mfrow=c(1,3))  
hist(p[1:(J/2+1)],breaks=seq(0,1,length=100),freq=F,col='lightblue',ylim=c(0,4),
     main='Good-precision',xlab='p',cex.lab=1.3)
hist(p[(J/2+1):J],breaks=seq(0,1,length=100),freq=F,col='lightblue',ylim=c(0,4),
     main='Poor-precision',xlab='p',cex.lab=1.3)
hist(p,breaks=seq(0,1,length=100),freq=F,col='lightblue',ylim=c(0,4),
     main='Combined',xlab='p',cex.lab=1.3)

dev.off()
# Analysing the good-precision data only
x1=qvalue(p[1:J/2+1],0.05)$qvalues
x2=locfdr(z[1:J/2+1],nulltype=0)$fdr
x3=ash_results(J,est_beta[1:J/2+1],est_s[1:J/2+1])
# Analysing the combined data 
y1=qvalue(p,0.05)$qvalues
y2=locfdr(z,nulltype=0)$fdr
y3=ash_results(J,est_beta[1:J/2+1],est_s[1:J/2+1])
x=seq(0,1,by=0.01)    # ideal case for comparison

## plots (Figure 4)
par(mfrow=c(1,3))  
plot(sort(x1),sort(y1),xlab='Good-precision data only',ylab='Combined data',
     main='qvalue',cex.lab=1.3)
lines(x,x,col='red')   # ideal case
plot(sort(x2),sort(y2),xlab='Good-precision data only',yaxt="n",ylab="",
     main='locfdr',cex.lab=1.3)
lines(x,x,col='red')   # ideal case
plot(sort(x3),sort(y3),xlab='Good-precision data only',yaxt="n",ylab="",
     main='ash',cex.lab=1.3)
lines(x,x,col='red')   # ideal case

dev.off()
