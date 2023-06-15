library(ggplot2)
library(grid)
library(gridExtra)
library(ggpubr)

#####################
#Simulation Function
#####################
#outliers: boolean for whether outlying cluster unit effects are simulated
#outlier_prop: if outliers=T, proportion of units that are outliers
#outlier_effect: if outliers=T, treatment effect coefficient for outlying units
sim=function(outliers=T,outlier_prop=0.1,outlier_effect=2) {
  nit=1000 #iterations
  B=100 
  cl=parallel::makeCluster(tail(numbers::divisors(B)[numbers::divisors(B)<=parallel::detectCores()],1))
  doParallel::registerDoParallel(cl)
  `%dopar%`=foreach::`%dopar%`
  result=foreach::foreach(vals=iterators::icount(nit),.combine=rbind,.errorhandling="remove") %dopar% {
    library(mvtnorm)
    library(MASS)
    source("Source.R")
    m=200 #number of centers
    mu=-6 #population norm
    ct=rep(1:m,each=100) #center ID
    nu=0.25
    W=rep(rnorm(m,0,1),each=100)
    W200=W[seq(1,100*m,100)]
    gamma=rep(0,100*m)
    if(outliers) {
      indices1=(100*m-100*outlier_prop*m+1):(100*m-50*outlier_prop*m)
      indices2=(100*m-50*outlier_prop*m+1):(100*m)
      gamma[indices1]=outlier_effect+0.5*W[indices1]
      gamma[indices2]=-outlier_effect+0.5*W[indices2]
    }
    P=rep(rnorm(m,-0.4,sqrt(0.5)),each=100)
    t=log(rep(rexp(m,0.001),each=100))
    alpha=rep(rnorm(m,0,sqrt(0.1)),each=100)
    Y=rpois(100*m,exp(mu+gamma+t+P+alpha+W*nu))
    Obs=sapply(split(Y,ct),sum)
    Exp=sapply(split(exp(mu+t+P),ct),sum)
    Z_FE=(Obs-Exp)/sqrt(Exp)
    #Estimators
    est1=mle_no_out(Obs,Exp,W200)
    nu1=est1[[1]]
    s2a1=est1[[2]]
    est2=nu_est(Obs,Exp,W200)
    nu2=est2[[1]]
    s2a2=est2[[2]]
    nu_vec=c(nu1,nu2,s2a1,s2a2)
    names(nu_vec)=c("Normal Nu","Truncated Normal Nu","Normal S2a","Truncated Normal S2a")
    nu_vec
  }
  parallel::stopCluster(cl)
  result
}

####################
#Outlier Simulations
####################
nu=0.25
s2a=0.1
outlier_props=seq(0,0.5,0.05)
bias=matrix(0,nrow=length(outlier_props),ncol=4)
for(i in 1:length(outlier_props)) {
  show(i)
  sim_result=sim(outliers=T,outlier_prop=outlier_props[i])
  bias[i,]=abs(apply(sim_result,2,mean)-c(nu,nu,s2a,s2a))
}
write.csv(bias,"Bias_Simulation_Outlier_Prop_2.csv")

nu=0.25
s2a=0.1
outlier_effects=seq(0,3.5,0.25)
bias=matrix(0,nrow=length(outlier_effects),ncol=4)
for(i in 1:length(outlier_effects)) {
  show(i)
  sim_result=sim(outliers=T,outlier_effect=outlier_effects[i])
  bias[i,]=abs(apply(sim_result,2,mean)-c(nu,nu,s2a,s2a))
}
write.csv(bias,"Bias_Simulation_Outlier_Effect_2.csv")
 
bias=read.csv("Bias_Simulation_Outlier_Prop_2.csv")
outlier_props=seq(0,0.5,0.05)
plot_data=data.frame(outlier_props=rep(outlier_props,2),bias=c(bias$V1,bias$V2),
                     method=c(rep("Normal MLE",length(outlier_props)),rep("RPP-CLC",length(outlier_props))))
g=ggplot(plot_data)+geom_line(aes(x=outlier_props,y=bias,linetype=method,colour=method),size=1)
g=g+scale_colour_manual(values=c("Normal MLE"="black","RPP-CLC"="dark green"))
g=g+scale_linetype_manual(values=c("Normal MLE"="solid","RPP-CLC"="dashed"))
g=g+theme_classic()+xlab("Outlier Proportion")+ylab("Absolute Bias")
g=g+ggtitle(expression(paste("(a) ",nu)))
g=g+ylim(c(0,2))
g=g+theme(text=element_text(size=20),legend.position="top",legend.title=element_blank(),legend.key.width=unit(3.5,"line"),legend.text=element_text(size=20))

bias=read.csv("Bias_Simulation_Outlier_Effect_2.csv")
outlier_effects=seq(0,3.5,0.25)
plot_data=data.frame(outlier_effects=rep(outlier_effects,2),bias=c(bias$V1,bias$V2),
                     method=c(rep("Normal MLE",length(outlier_effects)),rep("RPP-CLC",length(outlier_effects))))
g2=ggplot(plot_data)+geom_line(aes(x=outlier_effects,y=bias,linetype=method,colour=method),size=1)
g2=g2+scale_colour_manual(values=c("Normal MLE"="black","RPP-CLC"="dark green"))
g2=g2+scale_linetype_manual(values=c("Normal MLE"="solid","RPP-CLC"="dashed"))
g2=g2+theme_classic()+xlab("Outlier Quality of Care Effect Size")+ylab("Absolute Bias")
g2=g2+ggtitle(expression(paste("(b) ",nu)))
g2=g2+ylim(c(0,2))
g2=g2+theme(text=element_text(size=20),legend.position="top",legend.title=element_blank(),legend.key.width=unit(3.5,"line"),legend.text=element_text(size=20))

bias=read.csv("Bias_Simulation_Outlier_Prop_2.csv")
outlier_props=seq(0,0.5,0.05)
plot_data=data.frame(outlier_props=rep(outlier_props,2),bias=c(bias$V3,bias$V4),
                     method=c(rep("Normal MLE",length(outlier_props)),rep("RPP-CLC",length(outlier_props))))
g3=ggplot(plot_data)+geom_line(aes(x=outlier_props,y=bias,linetype=method,colour=method),size=1)
g3=g3+scale_colour_manual(values=c("Normal MLE"="black","RPP-CLC"="dark green"))
g3=g3+scale_linetype_manual(values=c("Normal MLE"="solid","RPP-CLC"="dashed"))
g3=g3+theme_classic()+xlab("Outlier Proportion")+ylab("Absolute Bias")
g3=g3+ggtitle(expression(paste("(c) ",sigma[alpha]^2)))
g3=g3+ylim(c(0,2))
g3=g3+theme(text=element_text(size=20),legend.position="top",legend.title=element_blank(),legend.key.width=unit(3.5,"line"),legend.text=element_text(size=20))

bias=read.csv("Bias_Simulation_Outlier_Effect_2.csv")
outlier_effects=seq(0,3.5,0.25)
plot_data=data.frame(outlier_effects=rep(outlier_effects,2),bias=c(bias$V3,bias$V4),
                     method=c(rep("Normal MLE",length(outlier_effects)),rep("RPP-CLC",length(outlier_effects))))
g4=ggplot(plot_data)+geom_line(aes(x=outlier_effects,y=bias,linetype=method,colour=method),size=1)
g4=g4+scale_colour_manual(values=c("Normal MLE"="black","RPP-CLC"="dark green"))
g4=g4+scale_linetype_manual(values=c("Normal MLE"="solid","RPP-CLC"="dashed"))
g4=g4+theme_classic()+xlab("Outlier Quality of Care Effect Size")+ylab("Absolute Bias")
g4=g4+ggtitle(expression(paste("(d) ",sigma[alpha]^2)))
g4=g4+ylim(c(0,2))
g4=g4+theme(text=element_text(size=20),legend.position="top",legend.title=element_blank(),legend.key.width=unit(3.5,"line"),legend.text=element_text(size=20))

ggarrange(g,g2,g3,g4,nrow=2,ncol=2,common.legend=T,legend="top")



