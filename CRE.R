library(ggplot2)
library(grid)
library(gridExtra)
library(ggpubr)

#####################
#Simulation Function
#####################
#outlier_prop: proportion of units that are outliers
#outlier_effect: treatment effect coefficient for outlying units
sim=function(outlier_prop=0.1,outlier_effect=2) {
  nit=5000 #iterations
  B=100 
  cl=parallel::makeCluster(tail(numbers::divisors(B)[numbers::divisors(B)<=parallel::detectCores()],1))
  doParallel::registerDoParallel(cl)
  `%dopar%`=foreach::`%dopar%`
  result=foreach::foreach(vals=iterators::icount(nit),.combine=rbind,.errorhandling="remove") %dopar% {
    library(mvtnorm)
    library(MASS)
    library(lme4)
    library(robustlmm)
    library(rlme)
    m=200 #number of centers
    mu=-6 #population norm
    ct=rep(1:m,each=100) #center ID
    source("Source.R")
    means=rep(rnorm(m,-0.4,sqrt(0.25)),each=100)
    X=rnorm(m*100,means,sqrt(0.25))
    Xbar=rep(tapply(X,ct,mean),each=100)
    Xdiff=X-Xbar
    alpha=rep(rnorm(m,0,sqrt(0.1)),each=100)
    tau=rep(rnorm(m,0,sqrt(0.1)),each=100)
    indices=(100*m-50*outlier_prop*m+1):(100*m)
    tau[indices]=rnorm(length(indices),-outlier_effect,sqrt(0.1))
    gamma=0.25*Xbar+tau
    Y=rnorm(100*m,mu+gamma+alpha+Xdiff)
    Obs=sapply(split(Y,ct),sum)
    temp=mu-0.25*0.4+Xdiff #0.25*0.4 to account for non-zero mean of Xbar term
    Exp=sapply(split(temp,ct),sum)
    cre=lmer(Y~Xdiff+Xbar+(1|ct))
    Xbar200=Xbar[seq(1,100*m,100)]
    Ess=rep(100,200)
    nu=nu_est_normal(Obs,Exp,Ess,Xbar200,1)
    result_vec=c(summary(cre)[[10]][,1][3],
                 nu$nu_hat)
    result_vec
  }
  parallel::stopCluster(cl)
  result
}

nu=0.25
outlier_props=seq(0,0.5,0.05)/2
bias=matrix(0,nrow=length(outlier_props),ncol=2)
sds=matrix(0,nrow=length(outlier_props),ncol=2)
see=matrix(0,nrow=length(outlier_props),ncol=2)
for(i in 1:length(outlier_props)) {
  show(i)
  sim_result=sim(outlier_prop=outlier_props[i],outlier_effect=2)
  bias[i,]=abs(apply(sim_result[,c(1,2)],2,mean)-nu)
  sds[i,]=apply(sim_result[,c(1,2)],2,sd)
}
mse=100*(bias^2+sds^2)
results=data.frame(bias=c(bias[,1],bias[,2]),mse=c(mse[,1],mse[,2]))
results$method=c(rep("CRE Model",length(outlier_props)),
                 rep("Proposed Two-Stage Model",length(outlier_props)))
write.csv(results,"CRE_Simulation_Outlier_Prop.csv")

plot_data=read.csv("CRE_Simulation_Outlier_Prop.csv")
plot_data$method=as.character(plot_data$method)
plot_data$method[plot_data$method=="CRE Model"]="CRE"
plot_data$method[plot_data$method=="Proposed Two-Stage Model"]="RPP-CLC"
plot_data$outlier_props=outlier_props
g=ggplot(plot_data)+geom_line(aes(x=outlier_props,y=bias,linetype=method,colour=method),size=1.25)
g=g+scale_colour_manual(values=c("CRE"="black",
                                 "RPP-CLC"="dark green"))
g=g+scale_linetype_manual(values=c("CRE"="solid",
                                   "RPP-CLC"="dashed"))
g=g+theme_classic()+xlab("Outlier Proportion")+ylab("Absolute Bias")
g=g+ylim(0,0.05)
g=g+theme(text=element_text(size=20),legend.position="top",legend.title=element_blank(),legend.key.width=unit(3.5,"line"),legend.text=element_text(size=20))
g=g+ggtitle("(a)")


plot_data=read.csv("CRE_Simulation_Outlier_Prop.csv")
plot_data$method=as.character(plot_data$method)
plot_data$method[plot_data$method=="CRE Model"]="CRE"
plot_data$method[plot_data$method=="Proposed Two-Stage Model"]="RPP-CLC"
plot_data$outlier_props=outlier_props
g2=ggplot(plot_data)+geom_line(aes(x=outlier_props,y=mse,linetype=method,colour=method),size=1.25)
g2=g2+scale_colour_manual(values=c("CRE"="black",
                                 "RPP-CLC"="dark green"))
g2=g2+scale_linetype_manual(values=c("CRE"="solid",
                                   "RPP-CLC"="dashed"))
g2=g2+theme_classic()+xlab("Outlier Proportion")+ylab("100MSE")
g2=g2+theme(text=element_text(size=20),legend.position="top",legend.title=element_blank(),legend.key.width=unit(3.5,"line"),legend.text=element_text(size=20))
g2=g2+ggtitle("(b)")

ggarrange(g,g2,nrow=1,common.legend=T)
























