#+++++++++++++++++++++++++++++++++++++++++++++++++++++#
#     Spline MDM                                      #
#     DLM for 9188A                                   #
#     Set/2011                                        #
#+++++++++++++++++++++++++++++++++++++++++++++++++++++#
library(chron);library(splines);library(car)
source("/home/terrahome/oaj9/_codes/MDM-VLaw-v12.r")
load("/home/terrahome/oaj9/_data/all2010/t5.9188A.Rdata")
load("/home/terrahome/oaj9/_data/all2010/t5.9193J.Rdata")
load("/home/terrahome/oaj9/_data/all2010/t5.1436M.Rdata")
load("/home/terrahome/oaj9/_data/all2010/stats/betas.9200B.Rdata")
load("/home/terrahome/oaj9/_data/all2010/stats/betas.9206B.Rdata")


#**************************************************************************************************#
#-------------------------------------9188A MODELLING----------------------------------------------#
#**************************************************************************************************#

time.filter<-6:20
test.dates<-c("18/08/2010","25/08/2010")
prior.months<- c("Jul","Aug")
#learning.months<-c("Sep","Oct","Nov")
learning.months<-c("Sep","Oct")
wkday<-"Wed"
#-------------------------#
#     Data Selection      #
#-------------------------#
# prior and learning data:
t5.9188A$Flow.idx<-rep(1:288,dim(table(t5.9188A$Date)))
c.pri.par<-months(t5.9188A$Date1) %in% prior.months & !t5.9188A$Date %in% test.dates & t5.9188A$Hour %in% time.filter
c.par<-(months(t5.9188A$Date1) %in% learning.months | t5.9188A$Date %in% test.dates)  & t5.9188A$Hour %in% time.filter  & t5.9188A$Dweek==wkday
period<-dim(table(t5.9188A[c.pri.par,]$Time5))
t5.9188A[is.na(t5.9188A$Hway_lag1)==T,]$Flow<-NA
t5.9188A[is.na(t5.9188A$Occ_lag1)==T,]$Flow<-NA
n.test<-period*length(test.dates)
with(t5.9188A[c.pri.par,],table(months(Date1)))
with(t5.9188A[c.par,],as.matrix(table(Date1)))


#----------#
#  priors  #
#----------#
prior.9188A<-list()
prior.9188A$sigma2<-with(t5.9188A[c.pri.par,], prior.5min(Flow,t5.9188A[c.pri.par,],dim(table(Date1)),period)$S0)
prior.9188A$mean.prior<-with(t5.9188A[c.pri.par & t5.9188A$Time5==t5.9188A[c.pri.par,]$Time5[1],],
                      mean(Flow,na.rm=T))
# total flow per period (denominator of the vlaw):
prior.9188A$total.flow<-sum(t5.9188A[t5.9188A$Date =="04/10/2010",]$Flow,na.rm=T)
#--------------------------------------#
#          Spline fitting              #
#--------------------------------------#
c.spl.par<-c.pri.par  & is.na(t5.9188A$Flow)==F & is.na(t5.9188A$Occ_lag1)==F  & is.na(t5.9188A$Hway_lag1)==F
c.spl.flow<-c.pri.par &  t5.9188A$Dweek==wkday & is.na(t5.9188A$Flow)==F
with(t5.9188A[c.spl.par,],table(months(Date1)))
with(t5.9188A[c.spl.flow,],table(months(Date1)))
prior.9188A$sp.flow<-with(t5.9188A[c.spl.flow,],ns(Flow.idx,knots=c(seq(40,130,by=5),seq(180,220,by=10))))
prior.9188A$sp.occ<-with(t5.9188A[c.spl.par,],ns(Occ_lag1,df=4))
prior.9188A$sp.hway<-with(t5.9188A[c.spl.par & t5.9188A$Hway_lag1>1,],ns(Hway_lag1,df=4))
prior.9188A$sp.spd<-with(t5.9188A[c.spl.par,],ns(Speed_lag1,df=4))
prior.9188A$sp.occsp<-with(t5.9188A[c.spl.par,],ns(Occ_lag1*Speed_lag1,df=4))


#--------------------------------------#
#          Spline prediction           #
#--------------------------------------#
# removing missing (otherwise predict function won't work - obs won't be used anyway)

x.flow<-t5.9188A[c.par,]$Flow.idx
for (i in 1:length(x.flow)) { if (is.na(x.flow[i])==T) x.flow[i]=x.flow[i-1]}
spl.flow<-predict(prior.9188A$sp.flow,x.flow)

x.occ<-t5.9188A[c.par,]$Occ_lag1
for (i in 1:length(x.occ)) { if (is.na(x.occ[i])==T) x.occ[i]=x.occ[i-1]}
spl.occ<-predict(prior.9188A$sp.occ,x.occ)

x.hway<-t5.9188A[c.par & t5.9188A$Hway_lag1>1,]$Hway_lag1
for (i in 1:length(x.hway)) { if (is.na(x.hway[i])==T) x.hway[i]=x.hway[i-1]}
spl.hway<-predict(prior.9188A$sp.hway,x.hway)

x.spd<-t5.9188A[c.par ,]$Speed_lag1
for (i in 1:length(x.spd)) { if (is.na(x.spd[i])==T) x.spd[i]=x.spd[i-1]}
spl.spd<-predict(prior.9188A$sp.spd,x.spd)

x.occsp<-t5.9188A[c.par,]$Occ_lag1*t5.9188A[c.par,]$Speed_lag1
for (i in 1:length(x.occsp)) { if (is.na(x.occsp[i])==T) x.occsp[i]=x.occsp[i-1]}
spl.occsp<-predict(prior.9188A$sp.occ,x.occsp)



#***************************#
#   9188A MODEL FITTING     #
#***************************#
dlm.9188A<-list()
VDF=rep(0.99,3)  #original results
#VDF=rep(0.96,3)   #optimzed



#   S:  Flow spline trend
dlm.9188A[[1]]<-MDM.VLaw.v12(info=t5.9188A[c.par,c("Flow","Date1","Hour")],nparents=0,beta.var=betas.rnd[[3]],Vdf=VDF,
                              X_t=as.matrix(spl.flow[1:nrow(spl.flow),]),prior.obj=prior.9188A,n.test=n.test)
#---------individual models----------------------------------#
#   SO  (S + occupancy)
dlm.9188A[[2]]<-MDM.VLaw.v12(info=t5.9188A[c.par,c("Flow","Date1","Hour")],,nparents=0,beta.var=betas.rnd[[3]],Vdf=VDF,
                             X_t=as.matrix(cbind(spl.flow[1:nrow(spl.flow),],spl.occ[1:nrow(spl.occ),])),
                             prior.obj=prior.9188A,n.test=n.test)
#   SH (S + headway)
dlm.9188A[[3]]<-MDM.VLaw.v12(info=t5.9188A[c.par,c("Flow","Date1","Hour")],,nparents=0,beta.var=betas.rnd[[3]],Vdf=VDF,
                             X_t=as.matrix(cbind(spl.flow[1:nrow(spl.flow),],spl.hway[1:nrow(spl.hway),])),
                             prior.obj=prior.9188A,n.test=n.test)
#   SHinv  (S + 1/headway)
dlm.9188A[[4]]<-MDM.VLaw.v12(info=t5.9188A[c.par,c("Flow","Date1","Hour")],,nparents=0,beta.var=betas.rnd[[3]],Vdf=VDF,
                             X_t=as.matrix(cbind(spl.flow[1:nrow(spl.flow),],1/t5.9188A[c.par,]$Hway_lag1)),
                             prior.obj=prior.9188A,n.test=n.test)
# - SSp (S + Speed)
dlm.9188A[[5]]<-MDM.VLaw.v12(info=t5.9188A[c.par,c("Flow","Date1","Hour")],nparents=0,beta.var=betas.rnd[[3]],Vdf=VDF,
                             X_t=as.matrix(cbind(spl.flow[1:nrow(spl.flow),],spl.spd[1:nrow(spl.spd),])),
                             prior.obj=prior.9188A,n.test=n.test)
#------Interaction test---------------------------------------#
# - SOSp: (S + occupancy + speed)
dlm.9188A[[6]]<-MDM.VLaw.v12(info=t5.9188A[c.par,c("Flow","Date1","Hour")],nparents=0,beta.var=betas.rnd[[3]],Vdf=VDF,
                             X_t=as.matrix(cbind(spl.flow[1:nrow(spl.flow),],spl.occ[1:nrow(spl.occ),],spl.spd[1:nrow(spl.spd),])),
                             prior.obj=prior.9188A,n.test=n.test)
# - S(OSp) : (S + occupancy + speed + occ x speed)
dlm.9188A[[7]]<-MDM.VLaw.v12(info=t5.9188A[c.par,c("Flow","Date1","Hour")],nparents=0,beta.var=betas.rnd[[3]],Vdf=VDF,
                             X_t=as.matrix(cbind(spl.flow[1:nrow(spl.flow),],spl.occ[1:nrow(spl.occ),],spl.spd[1:nrow(spl.spd),],
                                                  t5.9188A[c.par,]$Occ_lag1*t5.9188A[c.par,]$Speed_lag1)),
                             prior.obj=prior.9188A,n.test=n.test)
#------models with all variables:------------------------------#
# - Spline Flow Level + O + H + Sp
dlm.9188A[[8]]<-MDM.VLaw.v12(info=t5.9188A[c.par,c("Flow","Date1","Hour")],nparents=0,beta.var=betas.rnd[[3]],Vdf=VDF,
                             X_t=as.matrix(cbind(spl.flow[1:nrow(spl.flow),],spl.occ[1:nrow(spl.occ),],
                                                 spl.hway[1:nrow(spl.hway),],spl.spd[1:nrow(spl.spd),])),
                             prior.obj=prior.9188A,n.test=n.test)
# - Spline Flow Level +  Sp + Occ + (OSp) +  Hway
dlm.9188A[[9]]<-MDM.VLaw.v12(info=t5.9188A[c.par,c("Flow","Date1","Hour")],nparents=0,beta.var=betas.rnd[[3]],Vdf=VDF,
                             X_t=as.matrix(cbind(spl.flow[1:nrow(spl.flow),],spl.hway[1:nrow(spl.hway),],spl.occ[1:nrow(spl.occ),],
                                                  spl.spd[1:nrow(spl.spd),],t5.9188A[c.par,]$Occ_lag1*t5.9188A[c.par,]$Speed_lag1)),
                             prior.obj=prior.9188A,n.test=n.test)


save(dlm.9188A,file="/home/terrahome/oaj9/splines/dlm.9188A.Rdata")
save(prior.9188A,file="/home/terrahome/oaj9/splines/prior.9188A.Rdata")


(res2.9188A<-data.frame(Model=c("S","SO","SH","SHinv","SSp","SOSp","S(OSp)","SOSpH","S(OSp)H"),
            LPL=c(dlm.9188A[[1]]$perf[5,2],dlm.9188A[[2]]$perf[5,2], dlm.9188A[[3]]$perf[5,2],
                  dlm.9188A[[4]]$perf[5,2],dlm.9188A[[5]]$perf[5,2], dlm.9188A[[6]]$perf[5,2],
                  dlm.9188A[[7]]$perf[5,2],dlm.9188A[[8]]$perf[5,2], dlm.9188A[[9]]$perf[5,2])))

#**************************************************************************************************#
#-------------------------------------9193J MODELLING----------------------------------------------#
#**************************************************************************************************#
#load("/home/terrahome/oaj9/splines/dlm.9188A.Rdata")
#-------------------------#
#     Data Selection      #
#-------------------------#
t5.9193J$Flow.idx<-rep(1:288,dim(table(t5.9193J$Date)))
# prior and learning data:
period<-dim(table(t5.9193J[c.pri.chld,]$Time5))
c.pri.chld<-months(t5.9193J$Date1) %in% prior.months & !t5.9193J$Date %in% test.dates & t5.9193J$Hour %in% time.filter
c.chld<-(months(t5.9193J$Date1) %in% learning.months | t5.9193J$Date %in% test.dates) & t5.9193J$Hour %in% time.filter  & t5.9193J$Dweek==wkday
n.test<-period*length(test.dates)

with(t5.9193J[c.pri.chld,],table(months(Date1)))
with(t5.9193J[c.chld,],as.matrix(table(Date1)))

#----------#
#  priors  #
#----------#
prior.9193J<-MDM.DLM.find.prior(t5.9193J[c.pri.chld,]$Flow,t5.9188A[c.pri.par,]$Flow,period,0.001)
prior.9193J$sigma2<-sum((t5.9193J[c.pri.chld,]$Flow-t5.9188A[c.pri.par,]$Flow)^2,na.rm=T)/
                         length(t5.9193J[c.pri.chld,]$Flow)
prior.9193J$mean.prior<-with(t5.9193J[c.pri.chld & t5.9193J$Time5==t5.9193J[c.pri.chld,]$Time5[1],],
                      mean(Flow,na.rm=T))

# total flow per period (denominator of the vlaw):
prior.9193J$total.flow<-sum(t5.9193J[t5.9193J$Date =="04/10/2010",]$Flow,na.rm=T)

#--------------------------------------#
#          Spline fitting              #
#--------------------------------------#
c.chld.spl<-c.pri.chld & is.na(t5.9193J$Flow_lag1)==F  &  is.na(t5.9193J$Flow)==F & is.na(t5.9193J$Occ_lag1)==F &  is.na(t5.9193J$Hway_lag1)==F &  is.na(t5.9193J$Speed_lag1)==F
prior.9193J$sp.occ<-with(t5.9193J[c.chld.spl ,],ns(Occ_lag1,df=4))
prior.9193J$sp.hway<-with(t5.9193J[c.chld.spl & t5.9193J$Hway_lag1>1.3 ,],ns(Hway_lag1,df=4))
prior.9193J$sp.spd<-with(t5.9193J[c.chld.spl,],ns(Speed_lag1,df=4))
aux.prop<-data.frame(Flow.idx=t5.9193J[c.pri.chld,]$Flow.idx,Dweek=t5.9193J[c.pri.chld,]$Dweek,
                     Flow.par=t5.9188A[c.pri.par,]$Flow,Flow.chd=t5.9193J[c.pri.chld,]$Flow,
                     Flow.prop=t5.9193J[c.pri.chld,]$Flow/t5.9188A[c.pri.par,]$Flow)
aux.prop<-aux.prop[is.na(aux.prop$Flow.prop)==F & aux.prop$Dweek==wkday ,]
prior.9193J$sp.prop<-with(aux.prop,ns(Flow.idx,knots=c(seq(60,120,by=10),seq(180,220,by=10))))

#--------------------------------------#
#          Spline prediction           #
#--------------------------------------#
spl.prop<-predict(prior.9193J$sp.prop,t5.9193J[c.chld,]$Flow.idx)
# removing missing (otherwise predict function won't work - obs won't be used anyway)
x.occ<-t5.9193J[c.chld,]$Occ_lag1
for (i in 1:length(x.occ)) { if (is.na(x.occ[i])==T) x.occ[i]=x.occ[i-1]}
spl.occ<-predict(prior.9193J$sp.occ,x.occ)
x.hway<-t5.9193J[c.chld & t5.9193J$Hway_lag1>1.3,]$Hway_lag1
x.hway[1:204]<-x.hway[205:408]
for (i in 1:length(x.hway)) { if (is.na(x.hway[i])==T) x.hway[i]=x.hway[i-1]}
spl.hway<-predict(prior.9193J$sp.hway,x.hway)
x.spd<-t5.9193J[c.chld ,]$Speed_lag1
for (i in 1:length(x.spd)) { if (is.na(x.spd[i])==T) x.spd[i]=x.spd[i-1]}
spl.spd<-predict(prior.9193J$sp.spd,x.spd)


#***************************#
#   9193J MODEL FITTING     #
#***************************#
#  Spline proportion (S) and parent with spline flow (S) : S/S Model
VDF=rep(0.85,3)
MDM.9193J<-list()
occsp<-t5.9193J[c.chld,]$Occ_lag1*t5.9193J[c.chld,]$Speed_lag1
for (i in 1:length(occsp)) if (is.na(occsp[i])==T) occsp[i]<-occsp[i-1]
MDM.9193J[[1]]<-MDM.VLaw.v12(info=t5.9193J[c.chld,c("Flow","Date1","Hour")],nparents=1,beta.var=betas.ch.rnd[[3]],Vdf=VDF,
             X_t=t5.9188A[c.par,]$Flow*spl.prop[1:nrow(spl.prop),],F_t.star=spl.prop[1:nrow(spl.prop),],
             parent.f.mean=dlm.9188A[[1]]$f,parent.f.var=dlm.9188A[[1]]$Q,prior.obj=prior.9193J,n.test=n.test)
#  S/F Model
MDM.9193J[[2]]<-MDM.VLaw.v12(info=t5.9193J[c.chld,c("Flow","Date1","Hour")],nparents=1,beta.var=betas.ch.rnd[[3]],Vdf=VDF,
             X_t=t5.9188A[c.par,]$Flow*spl.prop[1:nrow(spl.prop),],F_t.star=spl.prop[1:nrow(spl.prop),],
             parent.f.mean=dlm.9188A[[9]]$f,parent.f.var=dlm.9188A[[9]]$Q,prior.obj=prior.9193J,n.test=n.test)
#   F/S Model
MDM.9193J[[3]]<-MDM.VLaw.v12(info=t5.9193J[c.chld,c("Flow","Date1","Hour")],nparents=1,beta.var=betas.ch.rnd[[3]],Vdf=VDF,
             X_t=as.matrix(cbind(t5.9188A[c.par,]$Flow*spl.prop[1:nrow(spl.prop),],spl.hway[1:nrow(spl.hway),],
             spl.occ[1:nrow(spl.occ),],spl.spd[1:nrow(spl.spd),],occsp)),F_t.star=spl.prop[1:nrow(spl.prop),],
             parent.f.mean=dlm.9188A[[1]]$f,parent.f.var=dlm.9188A[[1]]$Q,prior.obj=prior.9193J,n.test=n.test)
#  F/F Model
MDM.9193J[[4]]<-MDM.VLaw.v12(info=t5.9193J[c.chld,c("Flow","Date1","Hour")],nparents=1,beta.var=betas.ch.rnd[[3]], Vdf=VDF,
               X_t=as.matrix(cbind(t5.9188A[c.par,]$Flow*spl.prop[1:nrow(spl.prop),],spl.hway[1:nrow(spl.hway),],
               spl.occ[1:nrow(spl.occ),],spl.spd[1:nrow(spl.spd),],occsp)),F_t.star=spl.prop[1:nrow(spl.prop),],
               parent.f.mean=dlm.9188A[[9]]$f,parent.f.var=dlm.9188A[[9]]$Q,prior.obj=prior.9193J,n.test=n.test)
# DLM F/F
sp.flow.9193J<-with(aux.prop,ns(Flow.idx,knots=c(seq(40,100,by=10),seq(180,220,by=10))))
spl.flow<-predict(sp.flow.9193J,t5.9193J[c.chld,]$Flow.idx)
MDM.9193J[[5]]<-MDM.VLaw.v12(info=t5.9193J[c.chld,c("Flow","Date1","Hour")],nparents=0,beta.var=betas.ch.rnd[[3]], Vdf=VDF,
                X_t=as.matrix(cbind(spl.flow[1:nrow(spl.flow),],spl.hway[1:nrow(spl.hway),],
                spl.occ[1:nrow(spl.occ),],spl.spd[1:nrow(spl.spd),],occsp)),prior.obj=prior.9193J,n.test=n.test)


r<--(1:n.test)
#-----------------------#
# Joint LPL 9188A/9193J:#
#-----------------------#
aux.l<-length(t5.9193J[c.chld,]$Date[r])
jlpl.9193J<-data.frame(Date=t5.9193J[c.chld,]$Date[r],Hour=t5.9193J[c.chld,]$Hour[r],lpl.DLM=rep(NA,aux.l),
                       lpl.m1=rep(NA,aux.l),lpl.m2=rep(NA,aux.l),lpl.m3=rep(NA,aux.l),lpl.m4=rep(NA,aux.l))

#"S/S"
jlpl.9193J$lpl.m1<-dnorm(t5.9188A[c.par,]$Flow[r],mean=dlm.9188A[[1]]$f[r],sd=sqrt(dlm.9188A[[1]]$Q[r]),log=TRUE)+
                   dnorm(t5.9193J[c.chld,]$Flow[r],mean=MDM.9193J[[1]]$f[r],sd=sqrt(MDM.9193J[[1]]$Q[r]),log=TRUE)
#"S/F"
jlpl.9193J$lpl.m2<-dnorm(t5.9188A[c.par,]$Flow[r],mean=dlm.9188A[[9]]$f[r],sd=sqrt(dlm.9188A[[9]]$Q[r]),log=TRUE)+
                   dnorm(t5.9193J[c.chld,]$Flow[r],mean=MDM.9193J[[2]]$f[r],sd=sqrt(MDM.9193J[[2]]$Q[r]),log=TRUE)
#"F/S"
jlpl.9193J$lpl.m3<-dnorm(t5.9188A[c.par,]$Flow[r],mean=dlm.9188A[[1]]$f[r],sd=sqrt(dlm.9188A[[1]]$Q[r]),log=TRUE)+
                   dnorm(t5.9193J[c.chld,]$Flow[r],mean=MDM.9193J[[3]]$f[r],sd=sqrt(MDM.9193J[[3]]$Q[r]),log=TRUE)
#"F/F"
jlpl.9193J$lpl.m4<-dnorm(t5.9188A[c.par,]$Flow[r],mean=dlm.9188A[[9]]$f[r],sd=sqrt(dlm.9188A[[9]]$Q[r]),log=TRUE)+
                   dnorm(t5.9193J[c.chld,]$Flow[r],mean=MDM.9193J[[4]]$f[r],sd=sqrt(MDM.9193J[[4]]$Q[r]),log=TRUE)

#"F/F DLM"
jlpl.9193J$lpl.DLM<-dnorm(t5.9188A[c.par,]$Flow[r],mean=dlm.9188A[[9]]$f[r],sd=sqrt(dlm.9188A[[9]]$Q[r]),log=TRUE)+
                   dnorm(t5.9193J[c.chld,]$Flow[r],mean=MDM.9193J[[5]]$f[r],sd=sqrt(MDM.9193J[[5]]$Q[r]),log=TRUE)

(jlpl.9193J<-data.frame(model=c("S/S","S/F","F/S","F/F","F/F DLM"),
                           lpl=c(sum(jlpl.9193J$lpl.m1,na.rm="T"),sum(jlpl.9193J$lpl.m2,na.rm="T"),
                                 sum(jlpl.9193J$lpl.m3,na.rm="T"),sum(jlpl.9193J$lpl.m4,na.rm="T"),
                                 sum(jlpl.9193J$lpl.DLM,na.rm="T"))))

save(MDM.9193J,file="/home/terrahome/oaj9/splines/MDM.9193J.Rdata")
save(prior.9193J,file="/home/terrahome/oaj9/splines/prior.9193J.Rdata")

#**************************************************************************************************#
#-------------------------------------1436M MODELLING----------------------------------------------#
#**************************************************************************************************#
#load("/home/terrahome/oaj9/splines/MDM.9193J.Rdata")
#-------------------------#
#     Data Selection      #
#-------------------------#
t5.1436M$Flow.idx<-rep(1:288,dim(table(t5.1436M$Date)))
# prior and learning data:
c.pri.gchld<-months(t5.1436M$Date1) %in% prior.months & !t5.1436M$Date %in% test.dates & t5.1436M$Hour %in% time.filter
c.gchld<-(months(t5.1436M$Date1) %in% learning.months | t5.1436M$Date %in% test.dates) & t5.1436M$Hour %in% time.filter  & t5.1436M$Dweek==wkday
period<-dim(table(t5.1436M[c.pri.gchld,]$Time5))
n.test<-period*length(test.dates)
with(t5.1436M[c.pri.gchld,],table(months(Date1)))
with(t5.1436M[c.gchld,],as.matrix(table(Date1)))
#----------#
#  priors  #
#----------#
#prior.1436M<-MDM.DLM.find.prior(t5.1436M[c.pri.gchld,]$Flow,
#                                     t5.9193J[c.pri.chld,]$Flow, period,0.001)
prior.1436M<-list()
prior.1436M$sigma2<-sum((t5.1436M[c.pri.gchld,]$Flow-t5.9193J[c.pri.chld,]$Flow)^2,na.rm=T)/
                  length(t5.1436M[c.pri.gchld,]$Flow)
prior.1436M$mean.prior<-with(t5.1436M[c.pri.gchld & t5.1436M$Time5==t5.1436M[c.pri.gchld,]$Time5[1],],
                      mean(Flow,na.rm=T))
# total flow per period (denominator of the vlaw):
prior.1436M$total.flow<-sum(t5.1436M[t5.1436M$Date =="04/10/2010",]$Flow,na.rm=T)
#--------------------------------------#
#          Spline fitting              #
#--------------------------------------#
c.gchld.spl<-c.pri.gchld & is.na(t5.1436M$Flow_lag1)==F  &  is.na(t5.1436M$Flow)==F & is.na(t5.1436M$Occ_lag1)==F &  is.na(t5.1436M$Hway_lag1)==F &  is.na(t5.1436M$Speed_lag1)==F
prior.1436M$sp.occ<-with(t5.1436M[c.gchld.spl ,],ns(Occ_lag1,df=4))
prior.1436M$sp.hway<-with(t5.1436M[c.gchld.spl & t5.1436M$Hway_lag1>1.3 ,],ns(Hway_lag1,df=4))
prior.1436M$sp.spd<-with(t5.1436M[c.gchld.spl,],ns(Speed_lag1,df=4))
aux.prop<-data.frame(Flow.idx=t5.1436M[c.pri.gchld,]$Flow.idx,Dweek=t5.1436M[c.pri.gchld,]$Dweek,
                     Flow.par=t5.9193J[c.pri.gchld,]$Flow,Flow.chd=t5.1436M[c.pri.gchld,]$Flow,
                     Flow.prop=t5.1436M[c.pri.gchld,]$Flow/t5.9193J[c.pri.chld,]$Flow)
aux.prop<-aux.prop[is.na(aux.prop$Flow.prop)==F & aux.prop$Dweek==wkday,]
prior.1436M$sp.prop<-with(aux.prop,ns(Flow.idx,knots=c(seq(60,120,by=5),seq(180,220,by=5))))
#--------------------------------------#
#          Spline prediction           #
#--------------------------------------#
spl.prop<-predict(prior.1436M$sp.prop,t5.1436M[c.gchld,]$Flow.idx)
# removing missing (otherwise predict function won't work - obs won't be used anyway)
x.occ<-t5.1436M[c.gchld,]$Occ_lag1
for (i in 1:length(x.occ)) { if (is.na(x.occ[i])==T) x.occ[i]=x.occ[i-1]}
spl.occ<-predict(prior.1436M$sp.occ,x.occ)
x.hway<-t5.1436M[c.gchld,]$Hway_lag1
for (i in 1:length(x.hway)) { if (is.na(x.hway[i])==T) x.hway[i]=x.hway[i-1]}
spl.hway<-predict(prior.1436M$sp.hway,x.hway)
x.spd<-t5.1436M[c.par ,]$Speed_lag1
for (i in 1:length(x.spd)) { if (is.na(x.spd[i])==T) x.spd[i]=x.spd[i-1]}
spl.spd<-predict(prior.1436M$sp.spd,x.spd)


#***************************#
#   1436M MODEL FITTING     #
#***************************#
#S/S/S
MDM.1436M<-list()
occsp<-t5.1436M[c.gchld,]$Occ_lag1*t5.1436M[c.gchld,]$Speed_lag1
for (i in 1:length(occsp)) if (is.na(occsp[i])==T) occsp[i]<-occsp[i-1]
VDF=rep(0.89,3)
#source("/home/terrahome/oaj9/_codes/VLaw_Function-with_DLM_pack.r")
#S/S/S
MDM.1436M<-list()
MDM.1436M[[1]]<-MDM.VLaw.v12(info=t5.1436M[c.gchld,c("Flow","Date1","Hour")],nparents=1,beta.var=betas.ch.rnd[[3]],Vdf=VDF,
            X_t=t5.9193J[c.chld,]$Flow*spl.prop[1:nrow(spl.prop),],prior.obj=prior.1436M,n.test=n.test,
            parent.f.mean=MDM.9193J[[1]]$f.marg,parent.f.var=MDM.9193J[[1]]$marg.var,
            F_t.star=spl.prop[1:nrow(spl.prop),])

#S/S/F
MDM.1436M[[2]]<-MDM.VLaw.v12(info=t5.1436M[c.gchld,c("Flow","Date1","Hour")],nparents=1,beta.var=betas.ch.rnd[[3]],Vdf=VDF,
            X_t=t5.9193J[c.chld,]$Flow*spl.prop[1:nrow(spl.prop),],prior.obj=prior.1436M,n.test=n.test,
            parent.f.mean=MDM.9193J[[2]]$f.marg,parent.f.var=MDM.9193J[[2]]$marg.var,
            F_t.star=spl.prop[1:nrow(spl.prop),])
#S/F/F
MDM.1436M[[3]]<-MDM.VLaw.v12(info=t5.1436M[c.gchld,c("Flow","Date1","Hour")],nparents=1,beta.var=betas.ch.rnd[[3]],Vdf=VDF,
            X_t=t5.9193J[c.chld,]$Flow*spl.prop[1:nrow(spl.prop),],prior.obj=prior.1436M,n.test=n.test,
            parent.f.mean=MDM.9193J[[4]]$f.marg,parent.f.var=MDM.9193J[[4]]$marg.var,
            F_t.star=spl.prop[1:nrow(spl.prop),])
# F/F/F
MDM.1436M[[4]]<-MDM.VLaw.v12(info=t5.1436M[c.gchld,c("Flow","Date1","Hour")],nparents=1,beta.var=betas.ch.rnd[[3]],Vdf=VDF,
             X_t=as.matrix(cbind(t5.9193J[c.chld,]$Flow*spl.prop[1:nrow(spl.prop),],spl.hway[1:nrow(spl.hway),],
             spl.occ[1:nrow(spl.occ),],spl.spd[1:nrow(spl.spd),],occsp)),prior.obj=prior.1436M,n.test=n.test,
             parent.f.mean=MDM.9193J[[4]]$f.marg,parent.f.var=MDM.9193J[[4]]$marg.var,F_t.star=spl.prop[1:nrow(spl.prop),])

#single DLM with F:
sp.flow.1436M<-with(aux.prop,ns(Flow.idx,knots=c(seq(40,100,by=10),seq(180,220,by=10))))
spl.flow<-predict(sp.flow.1436M,t5.1436M[c.gchld,]$Flow.idx)
MDM.1436M[[5]]<-MDM.VLaw.v12(info=t5.1436M[c.gchld,c("Flow","Date1","Hour")],nparents=0,beta.var=betas.ch.rnd[[3]],Vdf=VDF,
   X_t=as.matrix(cbind(spl.flow[1:nrow(spl.flow),],spl.hway[1:nrow(spl.hway),],spl.occ[1:nrow(spl.occ),],
   spl.spd[1:nrow(spl.spd),],occsp)),prior.obj=prior.9193J,n.test=n.test)



r<--(1:n.test)
# Joint LPL 9188A/9193J/1436M:
aux.l<-length(t5.9193J[c.gchld,]$Date[r])
jlpl.1436M<-data.frame(Date=t5.9193J[c.chld,]$Date[r],Hour=t5.9193J[c.chld,]$Hour[r],lpl.DLM=rep(NA,aux.l),
                       lpl.m1=rep(NA,aux.l),lpl.m2=rep(NA,aux.l),lpl.m3=rep(NA,aux.l),lpl.m4=rep(NA,aux.l))
#"S/S/S"
jlpl.1436M$lpl.m1<-dnorm(t5.9188A[c.par,]$Flow[r],mean=dlm.9188A[[1]]$f[r],sd=sqrt(dlm.9188A[[1]]$Q[r]),log=TRUE)+
                   dnorm(t5.9193J[c.chld,]$Flow[r],mean=MDM.9193J[[1]]$f[r],sd=sqrt(MDM.9193J[[1]]$Q[r]),log=TRUE)+
                   dnorm(t5.1436M[c.gchld,]$Flow[r],mean=MDM.1436M[[1]]$f[r],sd=sqrt(MDM.1436M[[1]]$Q[r]),log=TRUE)
# "S/S/F"
jlpl.1436M$lpl.m2<-dnorm(t5.9188A[c.par,]$Flow[r],mean=dlm.9188A[[9]]$f[r],sd=sqrt(dlm.9188A[[9]]$Q[r]),log=TRUE)+
                   dnorm(t5.9193J[c.chld,]$Flow[r],mean=MDM.9193J[[2]]$f[r],sd=sqrt(MDM.9193J[[2]]$Q[r]),log=TRUE)+
                   dnorm(t5.1436M[c.gchld,]$Flow[r],mean=MDM.1436M[[2]]$f[r],sd=sqrt(MDM.1436M[[2]]$Q[r]),log=TRUE)
#"S/F/F"
jlpl.1436M$lpl.m3<-dnorm(t5.9188A[c.par,]$Flow[r],mean=dlm.9188A[[9]]$f[r],sd=sqrt(dlm.9188A[[9]]$Q[r]),log=TRUE)+
                   dnorm(t5.9193J[c.chld,]$Flow[r],mean=MDM.9193J[[4]]$f[r],sd=sqrt(MDM.9193J[[4]]$Q[r]),log=TRUE)+
                   dnorm(t5.1436M[c.gchld,]$Flow[r],mean=MDM.1436M[[3]]$f[r],sd=sqrt(MDM.1436M[[3]]$Q[r]),log=TRUE)
#"F/F/F"
jlpl.1436M$lpl.m4<-dnorm(t5.9188A[c.par,]$Flow[r],mean=dlm.9188A[[9]]$f[r],sd=sqrt(dlm.9188A[[9]]$Q[r]),log=TRUE)+
                   dnorm(t5.9193J[c.chld,]$Flow[r],mean=MDM.9193J[[4]]$f[r],sd=sqrt(MDM.9193J[[4]]$Q[r]),log=TRUE)+
                   dnorm(t5.1436M[c.gchld,]$Flow[r],mean=MDM.1436M[[4]]$f[r],sd=sqrt(MDM.1436M[[4]]$Q[r]),log=TRUE)

#"F/F/F DLM"
jlpl.1436M$lpl.DLM<-dnorm(t5.9188A[c.par,]$Flow[r],mean=dlm.9188A[[9]]$f[r],sd=sqrt(dlm.9188A[[9]]$Q[r]),log=TRUE)+
                   dnorm(t5.9193J[c.chld,]$Flow[r],mean=MDM.9193J[[5]]$f[r],sd=sqrt(MDM.9193J[[5]]$Q[r]),log=TRUE)+
                   dnorm(t5.1436M[c.gchld,]$Flow[r],mean=MDM.1436M[[5]]$f[r],sd=sqrt(MDM.1436M[[5]]$Q[r]),log=TRUE)

(jlpl.1436M<-data.frame(model=c("S/S/S","S/S/F","S/F/F","F/F/F","F/F/F DLM"),
                           lpl=c(sum(jlpl.1436M$lpl.m1,na.rm="T"),sum(jlpl.1436M$lpl.m2,na.rm="T"),
                                 sum(jlpl.1436M$lpl.m3,na.rm="T"),sum(jlpl.1436M$lpl.m4,na.rm="T"),
                                 sum(jlpl.1436M$lpl.DLM,na.rm="T"))))



save(MDM.1436M,file="/home/terrahome/oaj9/splines/MDM.1436M.Rdata")
save(prior.1436M,file="/home/terrahome/oaj9/splines/prior.1436M.Rdata")






