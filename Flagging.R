library(ggplot2)
library(grid)
library(gridExtra)
library(ggpubr)

#####################
#Simulation Function
#####################
#gammaVal: cluster unit treatment effect for unit of interest
#WVal: value for cluster unit confounding variable
#s: 1=flag in either direction, 2=flag in the true effect size direction
sim=function(gammaVal=0,WVal=0,s=1) {
  nit=1000 #iterations
  m=500 #number of centers
  nu=0.25
  mu=-6 #population norm
  ct=rep(1:m,each=100) #center ID
  outlier_prop=0.1 #outliers
  B=100 
  cl=parallel::makeCluster(tail(numbers::divisors(B)[numbers::divisors(B)<=parallel::detectCores()],1))
  doParallel::registerDoParallel(cl)
  `%dopar%`=foreach::`%dopar%`
  result=foreach::foreach(vals=iterators::icount(nit),.combine=rbind,.errorhandling="remove") %dopar% {
    library(mvtnorm)
    library(MASS)
    source("Source.R")
    gamma=rep(0,100*m)
    gamma[(100*m-100*outlier_prop*m+1):(100*m-50*outlier_prop*m)]=-2
    gamma[(100*m-50*outlier_prop*m+1):(100*m)]=2
    gamma[1:100]=gammaVal
    P=rep(rnorm(m,-0.4,sqrt(0.5)),each=100)
    t=log(rep(rexp(m,0.001),each=100))
    W=rep(rnorm(m,0,1),each=100)
    W[1:100]=WVal
    alpha=rep(rnorm(m,0,sqrt(0.1)),each=100)
    Y=rpois(100*m,exp(mu+gamma+t+P+alpha+W*nu))
    Obs=sapply(split(Y,ct),sum)
    Exp=sapply(split(exp(mu+t+P),ct),sum)
    W200=W[seq(1,100*m,100)]
    Z_FE=(Obs-Exp)/sqrt(Exp)
    W_pre=sqrt(Exp)*W200
    #Frequentist
    est=nu_est(Obs[2:m],Exp[2:m],W200[2:m])
    nu_hat=est[[1]]
    s2_hat=est[[2]]
    temp=exp(W200*nu_hat+s2_hat/2)
    mnull=sqrt(Exp)*(temp-1)
    vnull=temp*(1+temp*(exp(s2_hat)-1)*Exp)
    Z_EN=(Z_FE-mnull)/sqrt(vnull)
    if(s==1) {
      f1=(abs(Z_FE[1]) > 1.96)
      f2=(abs(Z_EN[1]) > 1.96)
    } else {
      f1=(Z_FE[1] < -1.96)
      f2=(Z_EN[1] < -1.96)
    }
    #Bayes
    lcl_srtr=qgamma(0.025,Obs[1]+2,Exp[1]+2)
    ucl_srtr=qgamma(0.975,Obs[1]+2,Exp[1]+2)
    lcl_prop=posterior_quantile_sim(0.025,Obs,Exp,W200)
    if(s==1) {
      f3=(ucl_srtr < 1 | lcl_srtr > 1)
      if(lcl_prop > 1) {
        f4=TRUE
      } else {
        ucl_prop=posterior_quantile_sim(0.975,Obs,Exp,W200)
        f4=(ucl_prop < 1)
      }
    } else {
      ucl_prop=posterior_quantile_sim(0.975,Obs,Exp,W200)
      f3=(ucl_srtr < 1)
      f4=(ucl_prop < 1)
    }
    flag=c(f1,f2,f3,f4)
    names(flag)=c("Naive Frequentist","Adjusted Frequentist","Naive Bayes","Adjusted Bayes")
    flag
  }
  parallel::stopCluster(cl)
  result
}

gamma_vals=-seq(0,3,0.25)
flag_result=matrix(0,nrow=length(gamma_vals),ncol=4)
for(i in 1:length(gamma_vals)) {
  show(i)
  sim_result=sim(gammaVal=gamma_vals[i],WVal=4,s=1)
  flag_result[i,]=apply(sim_result,2,mean)
}
write.csv(flag_result,"Flagging_Simulation_m=500.csv")
  
W_vals=-seq(0,12,1)
flag_result=matrix(0,nrow=length(W_vals),ncol=4)
for(i in 1:length(W_vals)) {
  show(i)
  sim_result=sim(gammaVal=0,WVal=W_vals[i],s=2)
  flag_result[i,]=apply(sim_result,2,mean)
}
write.csv(flag_result,"Flagging_Simulation_2_m=500.csv")

                                                  
get_plot1=function(flag_result,label1,label2) {
  plot_data=data.frame(flag=c(flag_result[,1],flag_result[,2],flag_result[,3],flag_result[,4]),
                       W_vals=rep(0.25*seq(0,12,1),4),
                       method=c(rep(c("Naive Frequentist","Adjusted Frequentist","Naive Pseudo-Bayesian","Adjusted Pseudo-Bayesian"),each=length(W_vals))))
  g=ggplot(plot_data)+geom_line(aes(x=W_vals,y=flag,linetype=method,colour=method),size=1.25)
  g=g+scale_colour_manual(values=c("Naive Frequentist"="black","Adjusted Frequentist"="blue",
                                   "Naive Pseudo-Bayesian"="dark green","Adjusted Pseudo-Bayesian"="dark orange"))
  g=g+scale_linetype_manual(values=c("Naive Frequentist"="solid","Adjusted Frequentist"="dotted",
                                     "Naive Pseudo-Bayesian"="dashed","Adjusted Pseudo-Bayesian"="dotdash"))
  g=g+theme_classic()+xlab("Observed Confounding Magnitude")+ylab("FFP")
  g=g+ggtitle(label1)
  g=g+theme(text=element_text(size=20),legend.position="top",legend.title=element_blank(),legend.key.width=unit(3.5,"line"),legend.text=element_text(size=20))
  g=g+guides(colour=guide_legend(nrow=2))
  return(g)
}

get_plot2=function(flag_result,label1,label2) {
  plot_data=data.frame(flag=c(flag_result[,1],flag_result[,2],flag_result[,3],flag_result[,4]),
                       gamma_vals=rep(seq(0,3,0.25),4),
                       method=c(rep(c("Naive Frequentist","Adjusted Frequentist","Naive Pseudo-Bayesian","Adjusted Pseudo-Bayesian"),each=length(gamma_vals))))
  g=ggplot(plot_data)+geom_line(aes(x=gamma_vals,y=flag,linetype=method,colour=method),size=1.25)
  g=g+scale_colour_manual(values=c("Naive Frequentist"="black","Adjusted Frequentist"="blue",
                                   "Naive Pseudo-Bayesian"="dark green","Adjusted Pseudo-Bayesian"="dark orange"))
  g=g+scale_linetype_manual(values=c("Naive Frequentist"="solid","Adjusted Frequentist"="dotted",
                                     "Naive Pseudo-Bayesian"="dashed","Adjusted Pseudo-Bayesian"="dotdash"))
  g=g+theme_classic()+xlab("Quality of Care Effect Size")+ylab("TFP")
  g=g+ggtitle(label2)
  g=g+theme(text=element_text(size=20),legend.position="top",legend.title=element_blank(),legend.key.width=unit(3.5,"line"),legend.text=element_text(size=20))
  g=g+guides(colour=guide_legend(nrow=2))
  return(g)
}

get_plots=function(name1,name2,label1,label2) {
  flag_result=read.csv(name1,head=T)
  flag_result=flag_result[,-1]
  gg1=get_plot1(flag_result,label1,label2)
  flag_result=read.csv(name2,head=T)
  flag_result=flag_result[,-1]
  gg2=get_plot2(flag_result,label1,label2)
  return(list(gg1,gg2))
}

space=ggplot()+theme_void()
g1=get_plots("Flagging_Simulation_2_m=50.csv","Flagging_Simulation_m=50.csv",
          label1="(a) I=50",label2="(b) I=50")
g2=get_plots("Flagging_Simulation_2.csv","Flagging_Simulation.csv",
          label1="(c) I=200",label2="(d) I=200")
g3=get_plots("Flagging_Simulation_2_m=500.csv","Flagging_Simulation_m=500.csv",
          label1="(e) I=500",label2="(f) I=500")
ggarrange(g1[[1]],g1[[2]],
          space,space,
          g2[[1]],g2[[2]],
          space,space,
          g3[[1]],g3[[2]],nrow=5,ncol=2,
          heights=c(1,0.1,1,0.1,1),widths=c(1,1),common.legend=T)


g1=get_plots("Flagging_Simulation_2_m=50.csv","Flagging_Simulation_m=50.csv",
             label1="(a) I=50",label2="(b) I=50")
g2=get_plots("Flagging_Simulation_2.csv","Flagging_Simulation.csv",
             label1="(a) I=200",label2="(b) I=200")
g3=get_plots("Flagging_Simulation_2_m=500.csv","Flagging_Simulation_m=500.csv",
             label1="(a) I=500",label2="(b) I=500")
ggarrange(g1[[1]],g1[[2]],widths=c(1,1),common.legend=T)
ggarrange(g2[[1]],g2[[2]],widths=c(1,1),common.legend=T)
ggarrange(g3[[1]],g3[[2]],widths=c(1,1),common.legend=T)


