compare1=function(estim,true){
  rango=range(c(true,estim))
  plot(true,estim,ylim=rango,xlim=rango)
  lines(rango,rango,col='red',lwd=2)  
}

plot(mod1$llk,type='l')
seq1=500:ngibbs
plot(mod1$llk[seq1],type='l')

#compare betas
betas.estim=matrix(mod1$betas[ngibbs,],nparam.occ,nspp)
compare1(estim=betas.estim,true=betas.true)

#compare m.betas
m.betas.estim=mod1$m.betas[ngibbs,]
compare1(estim=m.betas.estim,true=m.betas1.true)

#compare tau2.betas1
tau2.betas.estim=mod1$tau2.betas[ngibbs,]
compare1(estim=tau2.betas.estim,true=tau2.betas1.true)

#compare gammas
gammas.estim=matrix(mod1$gammas[ngibbs,],nparam.det,nspp)
compare1(estim=gammas.estim,true=gammas.true)

#compare m.gammas
m.gammas.estim=mod1$m.gamma[ngibbs,]
compare1(estim=m.gammas.estim,true=m.gamma.true)

#compare tau2.gammas
tau2.gammas.estim=mod1$tau2.gamma[ngibbs,]
compare1(estim=tau2.gammas.estim,true=tau2.gamma.true)