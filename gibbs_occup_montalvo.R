gibbs_occup=function(y,xmat.occ,xmat.det,tau2.a,tau2.b,
                     gamma1,ngibbs,nburn){

  #initial values
  nspp=ncol(y[,,1])
  nloc=nrow(xmat.occ)
  nparam.occ=ncol(xmat.occ)
  nparam.det=ncol(xmat.det[,,1])
  ystar=y
  cond=!is.na(y) & y==1; ystar[cond]=1
  cond=!is.na(y) & y==0; ystar[cond]=-1
  m.betas=rep(0,nparam.occ)
  tau2.betas=rep(1,nparam.occ)
  betas=matrix(0,nparam.occ,nspp)
  z=apply(y,c(1,2),max,na.rm=T)
  zstar=matrix(ifelse(z==1,1,-1),nloc,nspp)
  m.gamma=rep(0,nparam.det)
  tau2.gamma=rep(1,nparam.det)
  gammas=matrix(0,nparam.det,nspp)

  #MCMC settings
  store.m.betas=matrix(NA,ngibbs,nparam.occ)
  store.tau2.betas=matrix(NA,ngibbs,nparam.occ)
  store.gammas=matrix(NA,ngibbs,nparam.det*nspp)  
  store.m.gamma=matrix(NA,ngibbs,nparam.det)
  store.tau2.gamma=matrix(NA,ngibbs,nparam.det)
  store.betas=matrix(NA,ngibbs,nparam.occ*nspp)
  store.llk=matrix(NA,ngibbs,1)
  
  options(warn=2)

  #start gibbs sampler
  for (i in 1:ngibbs){
    print(i)

    #update latent variables
    z=sample.z(xmat.occ=xmat.occ,betas=betas,nloc=nloc,
               nspp=nspp,nrep=nrep,y=y,xmat.det=xmat.det,gammas=gammas)
    # z=z.true
    zstar=sample.zstar(z=z,xmat.occ=xmat.occ,betas=betas,nloc=nloc,nspp=nspp)
    # zstar=zstar.true
    
    ystar=sample.ystar(nrep=nrep,xmat.det=xmat.det,gammas=gammas,
                       y=y,nspp=nspp)
    # ystar=ystar.true

    #sample betas and associated prior parameters
    betas=sample.betas(m.betas=m.betas,tau2.betas=tau2.betas,
                       xmat.occ=xmat.occ,
                       zstar=zstar,nspp=nspp,nparam.occ=nparam.occ,
                       xtx.occ=xtx.occ)

    m.betas=sample.m.betas(nspp=nspp,betas=betas,tau2.betas=tau2.betas,
                           nparam.occ=nparam.occ)
    # m.betas1=m.betas1.true
    
    tau2.betas=sample.tau2.betas(nspp=nspp,betas=betas,
                                 nparam.occ=nparam.occ,
                                 m.betas=m.betas,
                                 tau2.a=tau2.a,tau2.b=tau2.b)

    #update gammas and associated prior parameters
    gammas=sample.gammas(ystar=ystar,xmat.det=xmat.det,z=z,m.gamma=m.gamma,tau2.gamma=tau2.gamma,
                         nparam.det=nparam.det,nspp=nspp,nrep=nrep)
    m.gamma=sample.m.gamma(gammas=gammas,tau2.gamma=tau2.gamma,nspp=nspp,nparam.det=nparam.det)
    tau2.gamma=sample.tau2.gamma(gammas=gammas,m=m.gamma,nspp=nspp,tau2.a=tau2.a,
                                 tau2.b=tau2.b,nparam.det=nparam.det)
    
    llk=get.llk(nloc=nloc,nspp=nspp,betas=betas,
                xmat.occ=xmat.occ,xmat.det=xmat.det,y=y,gammas=gammas)  
      
    #store results
    store.betas[i,]=betas
    store.m.betas[i,]=m.betas
    store.tau2.betas[i,]=tau2.betas
    store.gammas[i,]=gammas    
    store.m.gamma[i,]=m.gamma
    store.tau2.gamma[i,]=tau2.gamma
    store.llk[i]=llk
  }
  
  list(betas=store.betas,m.betas=store.m.betas,
       tau2.betas=store.tau2.betas,
       gammas=store.gammas,m.gamma=store.m.gamma,tau2.gamma=store.tau2.gamma,
       llk=store.llk)
}