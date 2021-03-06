rm(list=ls())
set.seed(32)

#basic settings
nloc=1000
nrep=5
nspp=150
nparam.det=3

#OCCUPANCY

#parameters
nparam.occ1=4

#get betas1
m.betas1.true=m.betas1=runif(nparam.occ1,min=-1,max=1)
tau2.betas1.true=tau2.betas1=runif(nparam.occ1,min=0,max=1)
betas1=matrix(NA,nparam.occ1,nspp)
for (i in 1:nparam.occ1){
  betas1[i,]=rnorm(nspp,mean=m.betas1.true[i],sd=sqrt(tau2.betas1[i]))
}
betas1.true=betas1

#get covariates for occupancy
xmat.occ=matrix(rnorm(nloc*(nparam.occ1-1)),nloc,nparam.occ1-1)
xmat.occ=cbind(1,xmat.occ)

#generate occupancy status
betas.true=betas=betas1
media=xmat.occ%*%betas
zstar.true=zstar=matrix(rnorm(nloc*nspp,mean=media,sd=1),nloc,nspp)
z.true=z=matrix(ifelse(zstar>0,1,0),nloc,nspp)

#---------------------------------
#---------------------------------
#DETECTION
tau2.gamma.true=tau2.gamma=runif(nparam.det,min=0,max=0.3)
m.gamma.true=m.gamma=runif(nparam.det,min=-1,max=1)
gammas=matrix(NA,nparam.det,nspp)
for (i in 1:nparam.det){
  gammas[i,]=rnorm(nspp,mean=m.gamma[i],sd=sqrt(tau2.gamma[i]))
}
gammas.true=gammas

#get covariates for detection
xmat.det=array(NA,dim=c(nloc,nparam.det,nrep))
for (i in 1:nrep){
  tmp=rnorm(nloc*(nparam.det-1))
  tmp1=cbind(1,matrix(tmp,nloc,nparam.det-1))
  xmat.det[,,i]=tmp1
}

#generate observations
media=ystar=y=array(NA,dim=c(nloc,nspp,nrep))
for (i in 1:nrep){
  media[,,i]=xmat.det[,,i]%*%gammas
  ystar[,,i]=rnorm(nloc*nspp,mean=media[,,i],sd=1)
  #site has to be occupied and species has to have been detected
  y[,,i]=ifelse(ystar[,,i]>0 & z.true==1,1,0) 
}
ystar.true=ystar

#export data
setwd('U:\\GIT_models\\occup_std\\fake data')
y1=matrix(y,nloc*nspp*nrep,1)
write.csv(y1,'fake data y.csv',row.names=F)
write.csv(xmat.occ,'fake data xmat occ.csv',row.names=F)
xmat.det1=matrix(xmat.det,nloc*nparam.det*nrep,1)
write.csv(xmat.det1,'fake data xmat det.csv',row.names=F)