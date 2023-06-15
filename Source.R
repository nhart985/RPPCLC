
####################################
#MLE Assuming No Outliers (Normal)
#Obs: Observed Cluster-Unit Outcome
#Exp: Expected Cluster-Unit Outcome
#W: Confounding Variable
####################################
mle_no_out=function(Obs,Exp,W) {
  W=as.matrix(W)
  W=W-matrix(apply(W,2,mean),nrow=dim(W)[1],ncol=dim(W)[2],byrow=T)
  get_pdr=function(Wp,arg) {
    return((as.matrix(Wp)%*%arg[1:P])[,1])
  }
  get_mean_var=function(arg,ntilde,W) {
    temp=exp(get_pdr(W,arg)+arg[P+1]/2)
    m=sqrt(ntilde)*(temp-1)
    v=temp*(1+temp*(exp(arg[P+1])-1)*ntilde)
    return(list(m,v))
  }
  P=dim(W)[2]
  Z_FE=(Obs-Exp)/sqrt(Exp)
  W_pre=W*sqrt(Exp)
  sigma_est=lm(Z_FE~0+W_pre)
  cf=sigma_est$coefficients
  pdrW=get_pdr(W,cf)
  pdrWpre=get_pdr(W_pre,cf)
  varphi_init=(summary(sigma_est)[[6]]^2-1-mean(pdrW))/mean(Exp)
  initial=c(cf,varphi_init)
  
  negloglik=function(arg) {
    mean_var=get_mean_var(arg,Exp,W)
    m=mean_var[[1]]
    v=pmax(mean_var[[2]],0.0001)
    loglik=sum(log(dnorm(Z_FE,mean=m,sd=sqrt(v))))
    return(-loglik)
  }
  
  result.optim=optim(par=initial,fn=negloglik)
  
  return(list(nu=result.optim$par[1:P],sig2a=result.optim$par[P+1]))
}

######################################################################
#MLE With Outliers for Normal Outcome Distribution
#Obs: Observed Cluster-Unit Outcome
#Exp: Expected Cluster-Unit Outcome
#Ess: Effective Sample Size (number of observations per unit)
#W: Confounding Variable
#s2e: Error Variance Parameter
#p.grid: null proportion grid
#cutoff: percentile of standard normal distribution for null intervals
######################################################################
nu_est_normal=function(Obs,Exp,Ess,W,s2e,p.grid=seq(0.8,0.999,0.001),cutoff=0.975) {
  W=as.matrix(W)
  W=W-matrix(apply(W,2,mean),nrow=dim(W)[1],ncol=dim(W)[2],byrow=T)
  get_pdr=function(Wp,arg) {
    return((as.matrix(Wp)%*%arg[1:P])[,1])
  }
  get_mean_var=function(arg,Ess,W) {
    m=sqrt(Ess)*get_pdr(W,arg)
    v=1+arg[P+1]*Ess
    return(list(m,v))
  }
  P=dim(W)[2]
  Z_FE=(Obs-Exp)/sqrt(s2e*Ess)
  cval=qnorm(cutoff)
  niter=length(p.grid)
  N=length(Z_FE)
  eval.p.grid=rep(0,niter)
  W_pre=W*sqrt(Ess)
  sigma_est=rlm(Z_FE~0+W_pre,method="M",scale.est="MAD",psi=psi.huber)
  nu_init=sigma_est$coefficients
  pdrW=get_pdr(W,nu_init)
  varphi_init=(sigma_est$s^2-1)/mean(Ess)
  initial=c(nu_init,varphi_init)
  m_init=get_pdr(W_pre,initial)
  v_init=1+varphi_init*Ess
  aorig=m_init-cval*sqrt(v_init)
  borig=m_init+cval*sqrt(v_init)
  
  negloglik=function(arg) {
    Z_FE0=Z_FE[which(Z_FE>=aorig & Z_FE<=borig)]
    Ess0=Ess[which(Z_FE>=aorig & Z_FE<=borig)]
    W0=W[which(Z_FE>=aorig & Z_FE<=borig),]
    null_mean_var=get_mean_var(arg,Ess0,W0)
    m0=null_mean_var[[1]]
    v0=null_mean_var[[2]]
    N0=length(Ess0)
    
    Ess1=Ess[-which(Z_FE>=aorig & Z_FE<=borig)]
    W1=W[-which(Z_FE>=aorig & Z_FE<=borig),]
    out_mean_var=get_mean_var(arg,Ess1,W1)
    m1=out_mean_var[[1]]
    v1=out_mean_var[[2]]
    aorig1=aorig[-which(Z_FE>=aorig & Z_FE<=borig)]
    borig1=borig[-which(Z_FE>=aorig & Z_FE<=borig)]
    
    v0=pmax(v0,0.0001)
    v1=pmax(v1,0.0001)
    Q=pnorm(borig1,mean=m1,sd=sqrt(v1))-pnorm(aorig1,mean=m1,sd=sqrt(v1))
    loglik=N0*log(p0)+sum(log(dnorm(Z_FE0,mean=m0,sd=sqrt(v0))))+sum(log(1-p0*Q))
    return(-loglik)
  }
  for(i in 1:niter) {
    p0=p.grid[i]
    result.optim=optim(par=initial,fn=negloglik)
    eval.p.grid[i]=result.optim$value
  }
  p0=p.grid[which.min(eval.p.grid)]
  result.optim=optim(par=initial,fn=negloglik)
  
  for(i in 1:niter) {
    p0=p.grid[i]
    result.optim=optim(par=initial,fn=negloglik)
    eval.p.grid[i]=result.optim$value
  }
  p0=p.grid[which.min(eval.p.grid)]
  result.optim=optim(par=initial,fn=negloglik)
  
  Ess0=Ess[which(Z_FE>=aorig & Z_FE<=borig)]
  W0=W[which(Z_FE>=aorig & Z_FE<=borig),]
  W_pre0=W_pre[which(Z_FE>=aorig & Z_FE<=borig),]
  pdrW0=get_pdr(W0,result.optim$par)
  Omega=diag(1+pdrW0+result.optim$par[P+1]*Ess0)
  nu_var_est=solve(t(W_pre0)%*%W_pre0)%*%t(W_pre0)%*%Omega%*%W_pre0%*%solve(t(W_pre0)%*%W_pre0)
  return(list(nu_hat=result.optim$par[1:P]*sqrt(s2e),
              sig2_hat=result.optim$par[P+1]*s2e,
              p0=p0,
              nu_var=nu_var_est,
              W=W,
              which(Z_FE < aorig | Z_FE > borig)))
}

#######################################################################################
#MLE With Outliers for Poisson Outcome Distribution
#Obs: Observed Cluster-Unit Outcome
#Exp: Expected Cluster-Unit Outcome
#W: Confounding Variable
#p.grid: null proportion grid
#cutoff: percentile of standard normal distribution for null intervals
#family: Poisson Model or Approximate Poisson Model (based on Taylor Series Expansion)
########################################################################################
nu_est=function(Obs,Exp,W,p.grid=seq(0.4,0.999,0.001),cutoff=0.975,family="poisson") {
  W=as.matrix(W)
  W=W-matrix(apply(W,2,mean),nrow=dim(W)[1],ncol=dim(W)[2],byrow=T)
  get_pdr=function(Wp,arg) {
    return((as.matrix(Wp)%*%arg[1:P])[,1])
  }
  get_mean_var=function(arg,ntilde,W,family) {
    if(family=="poisson") {
      temp=exp(get_pdr(W,arg)+arg[P+1]/2)
      m=sqrt(ntilde)*(temp-1)
      v=temp*(1+temp*(exp(arg[P+1])-1)*ntilde)
    } else if(family=="approx poisson"){
      pdrW=get_pdr(W,arg)
      pdrWpre=get_pdr(W_pre,arg)
      m=pdrWpre
      v=1+pdrW+arg[P+1]*ntilde
    }
    return(list(m,v))
  } 
  P=dim(W)[2]
  Z_FE=(Obs-Exp)/sqrt(Exp)
  cval=qnorm(cutoff)
  niter=length(p.grid)
  N=length(Z_FE)
  eval.p.grid=rep(0,niter)
  W_pre=W*sqrt(Exp)
  sigma_est=rlm(Z_FE~0+W_pre,method="M",scale.est="MAD",psi=psi.huber)
  nu_init=sigma_est$coefficients
  pdrW=get_pdr(W,nu_init)
  varphi_init=(sigma_est$s^2-1-mean(pdrW))/mean(Exp)
  initial=c(nu_init,varphi_init)
  temp=exp(get_pdr(W,initial)+initial[P+1]/2)
  m_init=sqrt(Exp)*(temp-1)
  v_init=temp*(1+temp*(exp(initial[P+1])-1)*Exp)
  aorig=m_init-cval*sqrt(v_init)
  borig=m_init+cval*sqrt(v_init)
  
  negloglik=function(arg) {
    Z_FE0=Z_FE[which(Z_FE>=aorig & Z_FE<=borig)]
    Exp0=Exp[which(Z_FE>=aorig & Z_FE<=borig)]
    W0=W[which(Z_FE>=aorig & Z_FE<=borig),]
    null_mean_var=get_mean_var(arg,Exp0,W0,family)
    m0=null_mean_var[[1]]
    v0=null_mean_var[[2]]
    N0=length(Exp0)
    
    Exp1=Exp[-which(Z_FE>=aorig & Z_FE<=borig)]
    W1=W[-which(Z_FE>=aorig & Z_FE<=borig),]
    out_mean_var=get_mean_var(arg,Exp1,W1,family)
    m1=out_mean_var[[1]]
    v1=out_mean_var[[2]]
    aorig1=aorig[-which(Z_FE>=aorig & Z_FE<=borig)]
    borig1=borig[-which(Z_FE>=aorig & Z_FE<=borig)]
    
    v0=pmax(v0,0.0001)
    v1=pmax(v1,0.0001)
    Q=pnorm(borig1,mean=m1,sd=sqrt(v1))-pnorm(aorig1,mean=m1,sd=sqrt(v1))
    loglik=N0*log(p0)+sum(log(dnorm(Z_FE0,mean=m0,sd=sqrt(v0))))+sum(log(1-p0*Q))
    return(-loglik)
  }
  for(i in 1:niter) {
    p0=p.grid[i]
    result.optim=optim(par=initial,fn=negloglik)
    eval.p.grid[i]=result.optim$value
  }
  p0=p.grid[which.min(eval.p.grid)]
  result.optim=optim(par=initial,fn=negloglik)
  
  Exp0=Exp[which(Z_FE>=aorig & Z_FE<=borig)]
  W0=W[which(Z_FE>=aorig & Z_FE<=borig),]
  W_pre0=W_pre[which(Z_FE>=aorig & Z_FE<=borig),]
  pdrW0=get_pdr(W0,result.optim$par)
  Omega=diag(1+pdrW0+result.optim$par[P+1]*Exp0)
  nu_var_est=solve(t(W_pre0)%*%W_pre0)%*%t(W_pre0)%*%Omega%*%W_pre0%*%solve(t(W_pre0)%*%W_pre0)
  return(list(nu_hat=result.optim$par[1:P],
              sig2_hat=result.optim$par[P+1],
              p0=p0,
              nu_var=nu_var_est,
              W=W,
              which(Z_FE < aorig | Z_FE > borig)))
}

####################################
#Posterior Density Function
#TRRs: Transplant Rate Ratio values
#Obs: Observed Cluster-Unit Outcome
#Exp: Expected Cluster-Unit Outcome
#W: Confounding Variable
####################################
posterior=function(TRRs,Obs,Exp,W) {
  sol=nu_est(Obs,Exp,W)
  hat_cf=sol$nu_hat
  var_hat_cf=sol$nu_var
  W=as.matrix(sol$W)
  var_pri_cf=diag(1,dim(W)[2],dim(W)[2])
  s2_post=solve(solve(var_pri_cf)+solve(var_hat_cf))
  m_post=(var_pri_cf%*%solve(var_pri_cf+var_hat_cf)%*%hat_cf)[,1]
  m=(W%*%m_post)[,1]
  s2=diag(W%*%s2_post%*%t(W))
  f=function(expcf,TRR,Obs_i,Exp_i) {
    return(dgamma(TRR,Obs_i+2,expcf*Exp_i+2)*dlnorm(expcf,m_i,sqrt(s2_i+sol$sig2_hat)))
  }
  out=matrix(NA,nrow=length(Obs),ncol=length(TRRs))
  for(i in 1:length(Obs)) {
    show(i)
    m_i=m[i]
    s2_i=s2[i]
    for(j in 1:length(TRRs)) {
      out_ij=integrate(f,0,Inf,TRR=TRRs[j],Obs_i=Obs[i],Exp_i=Exp[i],stop.on.error=F)$value
      out[i,j]=out_ij
    }
  }
  return(out)
}

######################################################
#Posterior Density Function for One Specific Unit
#TRRs: Transplant Rate Ratio values
#Obs: Observed Cluster-Unit Outcome
#Exp: Expected Cluster-Unit Outcome
#W: Confounding Variable
#i: Index of the unit 
#######################################################
posterior_one_center=function(TRRs,Obs,Exp,W,i) {
  sol=nu_est(Obs,Exp,W)
  hat_cf=sol$nu_hat
  var_hat_cf=sol$nu_var
  W=as.matrix(sol$W)
  var_pri_cf=diag(1,dim(W)[2],dim(W)[2])
  s2_post=solve(solve(var_pri_cf)+solve(var_hat_cf))
  m_post=(var_pri_cf%*%solve(var_pri_cf+var_hat_cf)%*%hat_cf)[,1]
  m=(W%*%m_post)[,1]
  s2=diag(W%*%s2_post%*%t(W))
  f=function(expcf,TRR,Obs_i,Exp_i) {
    return(dgamma(TRR,Obs_i+2,expcf*Exp_i+2)*dlnorm(expcf,m_i,sqrt(s2_i+sol$sig2_hat)))
  }
  out=vector("numeric")
  m_i=m[i]
  s2_i=s2[i]
  for(j in 1:length(TRRs)) {
    out_ij=integrate(f,0,Inf,TRR=TRRs[j],Obs_i=Obs[i],Exp_i=Exp[i],stop.on.error=F)$value
    out[j]=out_ij
  }
  return(out)
}

###################################
#Posterior Mean
#Obs: Observed Cluster-Unit Outcome
#Exp: Expected Cluster-Unit Outcome
#W: Confounding Variable
###################################
posterior_mean=function(Obs,Exp,W) {
  post_fun=function(TRRs,Obs,Exp,W,i) {
    post=posterior_one_center(TRRs,Obs,Exp,W,i)
    return(TRRs*post)
  }
  out=vector("numeric")
  for(i in 1:length(Obs)) {
    show(i)
    out[i]=integrate(post_fun,0,Inf,Obs,Exp,W,i,stop.on.error=F)
  }
  return(out)
}

###################################
#Posterior CDF
#x: TRR value
#Obs: Observed Cluster-Unit Outcome
#Exp: Expected Cluster-Unit Outcome
#W: Confounding Variable
###################################
posterior_CDF=function(x,Obs,Exp,W,var_pri_cf=1) {
  sol=nu_est(Obs,Exp,W)
  hat_cf=sol$nu_hat
  var_hat_cf=sol$nu_var
  W=as.matrix(sol$W)
  Var_pre=diag(W%*%var_hat_cf%*%t(W))
  s2=1/((1/var_pri_cf)+(1/Var_pre))
  m=(var_pri_cf/(Var_pre+var_pri_cf))*((as.matrix(W)%*%hat_cf)[,1])
  f=function(expcf,TRR,Obs_i,Exp_i) {
    return(dgamma(TRR,Obs_i+2,expcf*Exp_i+2)*dlnorm(expcf,m_i,sqrt(s2_i+sol$sig2_hat)))
  }
  post_fun=function(TRRs,Obs_i,Exp_i) {
    out=rep(NA,length(TRRs))
    for(j in 1:length(TRRs)) {
      out[j]=integrate(f,0,Inf,TRR=TRRs[j],Obs_i=Obs_i,Exp_i=Exp_i,stop.on.error=F)$value
    }
    return(out)
  }
  out=matrix(NA,nrow=length(Obs),ncol=length(x))
  for(i in 1:length(Obs)) {
    m_i=m[i]
    s2_i=s2[i]
    for(j in 1:length(x)) {
      out_ij=integrate(post_fun,0,x[j],Obs_i=Obs[i],Exp_i=Exp[i],stop.on.error=F)$value
      out[i,j]=out_ij
    }
  }
  return(out)
}

#######################################
#Posterior Quantile
#P: Percentile Value
#Obs: Observed Cluster-Unit Outcome
#Exp: Expected Cluster-Unit Outcome
#W: Confounding Variable
#var_pri_cf: Prior Coefficient Variance
#######################################
posterior_quantile=function(P,Obs,Exp,W,var_pri_cf=1) {
  sol=nu_est(Obs,Exp,W)
  hat_cf=sol$nu_hat
  var_hat_cf=sol$nu_var
  W=as.matrix(sol$W)
  Var_pre=diag(W%*%var_hat_cf%*%t(W))
  s2=1/((1/var_pri_cf)+(1/Var_pre))
  m=(var_pri_cf/(Var_pre+var_pri_cf))*((as.matrix(W)%*%hat_cf)[,1])
  f=function(expcf,TRR,Obs_i,Exp_i) {
    return(dgamma(TRR,Obs_i+2,expcf*Exp_i+2)*dlnorm(expcf,m_i,sqrt(s2_i+sol$sig2_hat)))
  }
  post_fun=function(TRRs,Obs_i,Exp_i) {
    out=rep(NA,length(TRRs))
    for(j in 1:length(TRRs)) {
      out[j]=integrate(f,0,Inf,TRR=TRRs[j],Obs_i=Obs_i,Exp_i=Exp_i,stop.on.error=F)$value
    }
    return(out)
  }
  posterior_CDF_fun=function(x,Obs_i,Exp_i) {
    out=rep(NA,length(x))
    for(k in 1:length(x)) {
      out[k]=integrate(post_fun,0,x,Obs_i=Obs_i,Exp_i=Exp_i,stop.on.error=F)$value
    }
    return(out)
  }
  out=rep(NA,length(Obs))
  for(i in 1:length(Obs)) {
    show(i)
    m_i=m[i]
    s2_i=s2[i]
    f_root=function(x) {
      return(posterior_CDF_fun(x,Obs_i=Obs[i],Exp_i=Exp[i])-P)
    }
    ratio=(Obs[i]+2)/(Exp[i]+2)
    out[i]=uniroot(f_root,c(0,100))$root
  }
  return(out)
}

#########################################
#Posterior Quantile for One Specific Unit 
#P: Percentile Value
#Obs: Observed Cluster-Unit Outcome
#Exp: Expected Cluster-Unit Outcome
#W: Confounding Variable
#var_pri_cf: Prior Coefficient Variance
#########################################
posterior_quantile_sim=function(P,Obs,Exp,W,var_pri_cf=1) {
  sol=nu_est(Obs,Exp,W)
  hat_cf=sol$nu_hat
  var_hat_cf=sol$nu_var
  W=as.matrix(sol$W)
  Var_pre=diag(W%*%var_hat_cf%*%t(W))
  s2=1/((1/var_pri_cf)+(1/Var_pre))
  m=(var_pri_cf/(Var_pre+var_pri_cf))*((as.matrix(W)%*%hat_cf)[,1])
  f=function(expcf,TRR,Obs_i,Exp_i) {
    return(dgamma(TRR,Obs_i+2,expcf*Exp_i+2)*dlnorm(expcf,m_i,sqrt(s2_i+sol$sig2_hat)))
  }
  post_fun=function(TRRs,Obs_i,Exp_i) {
    out=rep(NA,length(TRRs))
    for(j in 1:length(TRRs)) {
      out[j]=integrate(f,0,Inf,TRR=TRRs[j],Obs_i=Obs_i,Exp_i=Exp_i,stop.on.error=F)$value
    }
    return(out)
  }
  posterior_CDF_fun=function(x,Obs_i,Exp_i) {
    out=rep(NA,length(x))
    for(k in 1:length(x)) {
      out[k]=integrate(post_fun,0,x,Obs_i=Obs_i,Exp_i=Exp_i,stop.on.error=F)$value
    }
    return(out)
  }
  i=1
  m_i=m[i]
  s2_i=s2[i]
  f_root=function(x) {
    return(posterior_CDF_fun(x,Obs_i=Obs[i],Exp_i=Exp[i])-P)
  }
  ratio=(Obs[i]+2)/(Exp[i]+2)
  out=uniroot(f_root,c(0,20))$root
  return(out)
}

