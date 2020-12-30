# rm(list=ls())
library('mvtnorm')
set.seed(35)

#basic settings
nrep=5

#get functions
setwd('U:\\GIT_models\\occup_std')
source('aux occup montalvo.R')
source('gibbs_occup_montalvo.R')

#get design matrix for occupancy
setwd('U:\\GIT_models\\occup_std\\fake data')
xmat.occ=data.matrix(read.csv('fake data xmat occ.csv',as.is=T))
xtx.occ=t(xmat.occ)%*%xmat.occ
nparam.occ=ncol(xmat.occ)
nloc=nrow(xmat.occ)

#get data
tmp=read.csv('fake data y.csv',as.is=T)
nspp=nrow(tmp)/(nrep*nloc); nspp
y=array(tmp$V1,dim=c(nloc,nspp,nrep))

#design matrix for detection
tmp=read.csv('fake data xmat det.csv',as.is=T)
nparam.det=nrow(tmp)/(nloc*nrep); nparam.det
xmat.det=array(tmp$V1,dim=c(nloc,nparam.det,nrep))

#basic settings
ngibbs=1000
nburn=ngibbs/2

#priors
tau2.a=0.1; tau2.b=0.1
gamma1=0.01

#run gibbs
mod1=gibbs_occup(y=y,xmat.occ=xmat.occ,xmat.det=xmat.det,
                 tau2.a=tau2.a,tau2.b=tau2.b,gamma1=gamma1,
                 ngibbs=ngibbs,nburn=nburn)