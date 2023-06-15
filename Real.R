library(MASS)
source("Source.R")

#########
#CTR Data
#########
library(readxl)
CTR=as.data.frame(read_excel("csrs_final_tables_2111_KI.xls",sheet=6))
CTR=CTR[-1,]
Obs=as.numeric(CTR$TMR_Ad_Txed_c)
Exp=as.numeric(CTR$TMR_Ad_Txexp_c)
CT=substr(CTR$center,1,4)
LST=as.numeric(CTR$TMR_Ad_TxN_o)
DF=data.frame(CT,Obs,Exp,LST)

#########
#OPO Data
#########
OPO=as.data.frame(read_excel("OSR_final_tables2111.xlsx",sheet=1))
OPO2=as.data.frame(read_excel("OSR_final_tables2111.xlsx",sheet=4))
OPO3=as.data.frame(read_excel("OSR_final_tables2111.xlsx",sheet=3))
OPO4=as.data.frame(read_excel("OSR_final_tables2111.xlsx",sheet=2))
OPO2=OPO2[OPO2$`Begin date`=="07/01/2020",]
OPO=merge(merge(merge(OPO,OPO2,"OPO code"),OPO3,"OPO code"),OPO4,"OPO code")
OPO=data.frame("OPO"=OPO$`OPO code`,"CT"=OPO$`Center code`, 
               "DON"=OPO$`OPO observed KI donation rate per 100 eligible deaths`,
               "D"=OPO$`OPO deceased donors meeting eligibility criteria`,
               "DD"=OPO$`2019 Deaths`)
OPO$CT=substr(OPO$CT,1,4)

#######################
#Merge CTR and OPO Data
#######################
DF=merge(DF,OPO,by="CT")
DF=DF[order(DF$OPO),]
N_OPO=sapply(split(DF$OPO,DF$OPO),length)
Obs=DF$Obs
Exp=DF$Exp
W=DF$D/DF$LST

##########
#Estimates
##########
nu=nu_est(Obs,Exp,W,family="poisson")
nu$nu_hat
nu$sig2_hat
nu$p0
nu$nu_hat-1.96*sqrt(diag(nu$nu_var))
nu$nu_hat+1.96*sqrt(diag(nu$nu_var))

############################
#SRTR Posterior Distribution
############################
post_medians_srtr=qgamma(0.5,Obs+2,Exp+2)
post_q025_srtr=qgamma(0.025,Obs+2,Exp+2)
post_q975_srtr=qgamma(0.975,Obs+2,Exp+2)

#################################
#Proposed Posterior Distribution
#################################
post_medians_mod=posterior_quantile(0.5,Obs,Exp,W)
post_q025_mod=posterior_quantile(0.025,Obs,Exp,W)
post_q975_mod=posterior_quantile(0.975,Obs,Exp,W)

########################
#Save Posterior Results
########################
Dat=data.frame(post_medians_srtr,post_medians_mod,post_q025_srtr,post_q025_mod,post_q975_srtr,post_q975_mod)
write.csv(Dat,"Result.csv")

#######################
#Plot Posterior Results
#######################
TRRs=seq(0,2.5,0.005)
post=posterior(TRRs,Obs,Exp,W)
index=126
par(mai=c(1.25,1.25,0.82,0.42))
post_SRTR=dgamma(TRRs,Obs[index]+2,Exp[index]+2)
plot(post[index,]~TRRs,type="l",lwd=3,xlab="TRR Value",ylab="Posterior Density",lty=2,ylim=c(0,4),cex.axis=2,cex.lab=2,cex=2)
lines(post_SRTR~TRRs,lwd=3,col="blue")
legend("topright",legend=c("Original","Proposed"),lty=c(1,2),lwd=c(3,3),col=c("blue",1),cex=2)

#################
#Flagging Results
#################
flag_srtr=rep(0,length(post_q975_srtr))
flag_srtr[post_q975_srtr < 1]=-1
flag_srtr[post_q025_srtr > 1]=1

flag_mod=rep(0,length(post_q975_mod))
flag_mod[post_q975_mod < 1]=-1
flag_mod[post_q025_mod > 1]=1
table(flag_srtr,flag_mod)

which(flag_srtr==1 & flag_mod==-1)

##############
#OPO Dataset
##############
OPO_CD=DF$OPO
m_srtr=tapply(post_medians_srtr,OPO_CD,mean)
m_mod=tapply(post_medians_mod,OPO_CD,mean)
W1_opo=tapply(W1,OPO_CD,mean)
W2_opo=tapply(W2,OPO_CD,mean)
opo_ctr_cd=unique(OPO_CD)
O_MAP=data.frame(opo_ctr_cd,m_srtr,m_mod,W1_opo,W2_opo)
rownames(O_MAP)=NULL
write.csv(O_MAP,"O_MAP.csv")

dat=data.frame(OPO=DF$OPO,W,post=post_medians_srtr)
dat=dat[order(dat$W),]
dat$group="High"
dat$group[dat$W < quantile(unique(dat$W),1/3)]="Low"
dat$group[dat$W >= quantile(unique(dat$W),1/3) & dat$W < quantile(unique(dat$W),2/3)]="Medium"
dat$group=factor(dat$group,levels=c("Low","Medium","High"))
boxplot(post~group,data=dat,xlab="DSA Donor Organ Availability",ylab="TRR Posterior Median",cex.axis=1.5,cex.lab=1.5,cex=1.5)
abline(h=1,col="dark blue",lwd=3,lty=2)






