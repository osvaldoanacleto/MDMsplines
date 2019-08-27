#+++++++++++++++++++++++++++++++++++++++++++++++++++++#
#     Spline MDM                                      #
#     DLM for 1431A                                   #
#     Set/2011                                        #
#+++++++++++++++++++++++++++++++++++++++++++++++++++++#
library(chron);library(splines);library(car)
source("/home/terrahome/oaj9/_codes/MDM-VLaw-v12.r")
load("/home/terrahome/oaj9/_data/all2010/t5.1431A.Rdata")
load("/home/terrahome/oaj9/_data/all2010/t5.1437A.Rdata")
load("/home/terrahome/oaj9/_data/all2010/t5.1441A.Rdata")
#Defining periods
time.filter<-6:20
test.dates<-c("18/08/2010","25/08/2010")
prior.months<- c("Jul","Aug")
#learning.months<-c("Sep","Oct","Nov")
learning.months<-c("Sep","Oct")
wkday<-"Wed"

#**************************************************************************************************#
#-------------------------------------1431A MODELLING----------------------------------------------#
#**************************************************************************************************#

#-------------------------#
#     Data Selection      #
#-------------------------#
# prior and learning data:
t5.1431A$Flow.idx<-rep(1:288,dim(table(t5.1431A$Date)))
c.pri.par<-months(t5.1431A$Date1) %in% prior.months & !t5.1431A$Date %in% test.dates & t5.1431A$Hour %in% time.filter
c.par<-(months(t5.1431A$Date1) %in% learning.months | t5.1431A$Date %in% test.dates)  & t5.1431A$Hour %in% time.filter  & t5.1431A$Dweek==wkday
period<-dim(table(t5.1431A[c.pri.par,]$Time5))
t5.1431A[is.na(t5.1431A$Hway_lag1)==T,]$Flow<-NA
t5.1431A[is.na(t5.1431A$Occ_lag1)==T,]$Flow<-NA
n.test<-period*length(test.dates)

with(t5.1431A[c.par,],as.matrix(table(Date1)))
with(t5.1431A[c.pri.par,],table(months(Date1)))


#----------#
#  priors  #
#----------#
prior.1431A<-list()
prior.1431A$sigma2<-with(t5.1431A[c.pri.par,], prior.5min(Flow,t5.1431A[c.pri.par,],dim(table(Date1)),period)$S0)
prior.1431A$mean.prior<-with(t5.1431A[c.pri.par & t5.1431A$Time5==t5.1431A[c.pri.par,]$Time5[1],],
                      mean(Flow,na.rm=T))
# total flow per period (denominator of the vlaw):
prior.1431A$total.flow<-sum(t5.1431A[t5.1431A$Date =="04/10/2010",]$Flow,na.rm=T)
#--------------------------------------#
#          Spline fitting              #
#--------------------------------------#

#c.spl.par<-months(t5.1431A$Date1) %in% c("Jun") & !t5.1431A$Date %in% test.dates & t5.1431A$Hour %in% time.filter  &
# is.na(t5.1431A$Flow)==F & is.na(t5.1431A$Occ_lag1)==F  & is.na(t5.1431A$Hway_lag1)==F
c.spl.par<-c.pri.par  & is.na(t5.1431A$Flow)==F & is.na(t5.1431A$Occ_lag1)==F  & is.na(t5.1431A$Hway_lag1)==F
c.spl.flow<-c.pri.par &  t5.1431A$Dweek==wkday & is.na(t5.1431A$Flow)==F
with(t5.1431A[c.spl.par,],table(months(Date1)))
with(t5.1431A[c.spl.flow,],table(months(Date1)))
prior.1431A$sp.flow<-with(t5.1431A[c.spl.flow,],ns(Flow.idx,knots=c(seq(40,130,by=5),seq(180,220,by=10))))
prior.1431A$sp.occ<-with(t5.1431A[c.spl.par,],ns(Occ_lag1,df=4))
prior.1431A$sp.hway<-with(t5.1431A[c.spl.par & t5.1431A$Hway_lag1>1,],ns(Hway_lag1,df=4))
prior.1431A$sp.spd<-with(t5.1431A[c.spl.par,],ns(Speed_lag1,df=4))
prior.1431A$sp.occsp<-with(t5.1431A[c.spl.par,],ns(Occ_lag1*Speed_lag1,df=4))



#--------------------------------------#
#          Spline prediction           #
#--------------------------------------#
# removing missing (otherwise predict function won't work - obs won't be used anyway)

x.flow<-t5.1431A[c.par,]$Flow.idx
for (i in 1:length(x.flow)) { if (is.na(x.flow[i])==T) x.flow[i]=x.flow[i-1]}
spl.flow<-predict(prior.1431A$sp.flow,x.flow)

x.occ<-t5.1431A[c.par,]$Occ_lag1
for (i in 1:length(x.occ)) { if (is.na(x.occ[i])==T) x.occ[i]=x.occ[i-1]}
spl.occ<-predict(prior.1431A$sp.occ,x.occ)

x.hway<-t5.1431A[c.par & t5.1431A$Hway_lag1>1,]$Hway_lag1
for (i in 1:length(x.hway)) { if (is.na(x.hway[i])==T) x.hway[i]=x.hway[i-1]}
spl.hway<-predict(prior.1431A$sp.hway,x.hway)

x.spd<-t5.1431A[c.par ,]$Speed_lag1
for (i in 1:length(x.spd)) { if (is.na(x.spd[i])==T) x.spd[i]=x.spd[i-1]}
spl.spd<-predict(prior.1431A$sp.spd,x.spd)

x.occsp<-t5.1431A[c.par,]$Occ_lag1*t5.1431A[c.par,]$Speed_lag1
for (i in 1:length(x.occsp)) { if (is.na(x.occsp[i])==T) x.occsp[i]=x.occsp[i-1]}
spl.occsp<-predict(prior.1431A$sp.occ,x.occsp)



#***************************#
#   1431A MODEL FITTING     #
#***************************#



dlm.1431A<-list()
VDF=rep(0.99,3)  #original results

#   S:  Flow spline trend
dlm.1431A[[1]]<-MDM.VLaw.v12(info=t5.1431A[c.par,c("Flow","Date1","Hour")],nparents=0,beta.var=betas.rnd[[4]],Vdf=VDF,
                              X_t=as.matrix(spl.flow[1:nrow(spl.flow),]),prior.obj=prior.1431A,n.test=n.test)
#---------individual models----------------------------------#
#   SO  (S + occupancy)
dlm.1431A[[2]]<-MDM.VLaw.v12(info=t5.1431A[c.par,c("Flow","Date1","Hour")],,nparents=0,beta.var=betas.rnd[[4]],Vdf=VDF,
                             X_t=as.matrix(cbind(spl.flow[1:nrow(spl.flow),],spl.occ[1:nrow(spl.occ),])),
                             prior.obj=prior.1431A,n.test=n.test)
#   SH (S + headway)
dlm.1431A[[3]]<-MDM.VLaw.v12(info=t5.1431A[c.par,c("Flow","Date1","Hour")],,nparents=0,beta.var=betas.rnd[[4]],Vdf=VDF,
                             X_t=as.matrix(cbind(spl.flow[1:nrow(spl.flow),],spl.hway[1:nrow(spl.hway),])),
                             prior.obj=prior.1431A,n.test=n.test)
#   SHinv  (S + 1/headway)
dlm.1431A[[4]]<-MDM.VLaw.v12(info=t5.1431A[c.par,c("Flow","Date1","Hour")],,nparents=0,beta.var=betas.rnd[[4]],Vdf=VDF,
                             X_t=as.matrix(cbind(spl.flow[1:nrow(spl.flow),],1/t5.1431A[c.par,]$Hway_lag1)),
                             prior.obj=prior.1431A,n.test=n.test)
# - SSp (S + Speed)
dlm.1431A[[5]]<-MDM.VLaw.v12(info=t5.1431A[c.par,c("Flow","Date1","Hour")],nparents=0,beta.var=betas.rnd[[4]],Vdf=VDF,
                             X_t=as.matrix(cbind(spl.flow[1:nrow(spl.flow),],spl.spd[1:nrow(spl.spd),])),
                             prior.obj=prior.1431A,n.test=n.test)
#------Interaction test---------------------------------------#
# - SOSp: (S + occupancy + speed)
dlm.1431A[[6]]<-MDM.VLaw.v12(info=t5.1431A[c.par,c("Flow","Date1","Hour")],nparents=0,beta.var=betas.rnd[[4]],Vdf=VDF,
                             X_t=as.matrix(cbind(spl.flow[1:nrow(spl.flow),],spl.occ[1:nrow(spl.occ),],spl.spd[1:nrow(spl.spd),])),
                             prior.obj=prior.1431A,n.test=n.test)
# - S(OSp) : (S + occupancy + speed + occ x speed)
dlm.1431A[[7]]<-MDM.VLaw.v12(info=t5.1431A[c.par,c("Flow","Date1","Hour")],nparents=0,beta.var=betas.rnd[[4]],Vdf=VDF,
                             X_t=as.matrix(cbind(spl.flow[1:nrow(spl.flow),],spl.occ[1:nrow(spl.occ),],spl.spd[1:nrow(spl.spd),],
                                                  t5.1431A[c.par,]$Occ_lag1*t5.1431A[c.par,]$Speed_lag1)),
                             prior.obj=prior.1431A,n.test=n.test)
#------models with all variables:------------------------------#
# - Spline Flow Level + O + H + Sp
dlm.1431A[[8]]<-MDM.VLaw.v12(info=t5.1431A[c.par,c("Flow","Date1","Hour")],nparents=0,beta.var=betas.rnd[[4]],Vdf=VDF,
                             X_t=as.matrix(cbind(spl.flow[1:nrow(spl.flow),],spl.occ[1:nrow(spl.occ),],
                                                 spl.hway[1:nrow(spl.hway),],spl.spd[1:nrow(spl.spd),])),
                             prior.obj=prior.1431A,n.test=n.test)
# - Spline Flow Level +  Sp + Occ + (OSp) +  Hway
dlm.1431A[[9]]<-MDM.VLaw.v12(info=t5.1431A[c.par,c("Flow","Date1","Hour")],nparents=0,beta.var=betas.rnd[[4]],Vdf=VDF,
                             X_t=as.matrix(cbind(spl.flow[1:nrow(spl.flow),],spl.hway[1:nrow(spl.hway),],spl.occ[1:nrow(spl.occ),],
                                                  spl.spd[1:nrow(spl.spd),],t5.1431A[c.par,]$Occ_lag1*t5.1431A[c.par,]$Speed_lag1)),
                             prior.obj=prior.1431A,n.test=n.test)


save(dlm.1431A,file="/home/terrahome/oaj9/splines/dlm.1431A.Rdata")
save(prior.1431A,file="/home/terrahome/oaj9/splines/prior.1431A.Rdata")




(res2.1431A<-data.frame(Model=c("S","SO","SH","SHinv","SSp","SOSp","S(OSp)","SOSpH","S(OSp)H"),
            LPL=c(dlm.1431A[[1]]$perf[5,2],dlm.1431A[[2]]$perf[5,2], dlm.1431A[[3]]$perf[5,2],
                  dlm.1431A[[4]]$perf[5,2],dlm.1431A[[5]]$perf[5,2], dlm.1431A[[6]]$perf[5,2],
                  dlm.1431A[[7]]$perf[5,2],dlm.1431A[[8]]$perf[5,2], dlm.1431A[[9]]$perf[5,2])))
res.part2<-data.frame(model=res2.9206B[,1],lpl.9206B=res2.9206B[,2],lpl.6013B=res2.6013B[,2],
                                          lpl.9188A=res2.9188A[,2],lpl.1431A=res2.1431A[,2])
                                          
                                          
write.csv(res.part2,file="/home/terrahome/oaj9/splines/res.part2.csv")

#**************************************************************************************************#
#-------------------------------------1437A MODELLING----------------------------------------------#
#**************************************************************************************************#
#load("/home/terrahome/oaj9/splines/dlm.1431A.Rdata")
#-------------------------#
#     Data Selection      #
#-------------------------#
t5.1437A$Flow.idx<-rep(1:288,dim(table(t5.1437A$Date)))
# prior and learning data:
period<-dim(table(t5.1437A[c.pri.chld,]$Time5))
c.pri.chld<-months(t5.1437A$Date1) %in% prior.months & !t5.1437A$Date %in% test.dates & t5.1437A$Hour %in% time.filter
c.chld<-(months(t5.1437A$Date1) %in% learning.months | t5.1437A$Date %in% test.dates) & t5.1437A$Hour %in% time.filter  & t5.1437A$Dweek==wkday
n.test<-period*length(test.dates)

with(t5.1437A[c.pri.chld,],table(months(Date1)))
with(t5.1437A[c.chld,],as.matrix(table(Date1)))

#----------#
#  priors  #
#----------#
prior.1437A<-MDM.DLM.find.prior(t5.1437A[c.pri.chld,]$Flow,t5.1431A[c.pri.par,]$Flow,period,0.001)
prior.1437A$sigma2<-sum((t5.1437A[c.pri.chld,]$Flow-t5.1431A[c.pri.par,]$Flow)^2,na.rm=T)/
                         length(t5.1437A[c.pri.chld,]$Flow)
prior.1437A$mean.prior<-with(t5.1437A[c.pri.chld & t5.1437A$Time5==t5.1437A[c.pri.chld,]$Time5[1],],
                      mean(Flow,na.rm=T))

# total flow per period (denominator of the vlaw):
prior.1437A$total.flow<-sum(t5.1437A[t5.1437A$Date =="04/10/2010",]$Flow,na.rm=T)

#--------------------------------------#
#          Spline fitting              #
#--------------------------------------#
c.chld.spl<-c.pri.chld & is.na(t5.1437A$Flow_lag1)==F  &  is.na(t5.1437A$Flow)==F & is.na(t5.1437A$Occ_lag1)==F &  is.na(t5.1437A$Hway_lag1)==F &  is.na(t5.1437A$Speed_lag1)==F
prior.1437A$sp.occ<-with(t5.1437A[c.chld.spl ,],ns(Occ_lag1,df=4))
prior.1437A$sp.hway<-with(t5.1437A[c.chld.spl & t5.1437A$Hway_lag1>1.3 ,],ns(Hway_lag1,df=4))
prior.1437A$sp.spd<-with(t5.1437A[c.chld.spl,],ns(Speed_lag1,df=4))

aux.prop<-data.frame(Flow.idx=t5.1437A[c.pri.chld,]$Flow.idx,Dweek=t5.1437A[c.pri.chld,]$Dweek,
                     Flow.par=t5.1431A[c.pri.par,]$Flow,Flow.chd=t5.1437A[c.pri.chld,]$Flow,
                     Flow.prop=t5.1437A[c.pri.chld,]$Flow/t5.1431A[c.pri.par,]$Flow)
aux.prop<-aux.prop[is.na(aux.prop$Flow.prop)==F & aux.prop$Dweek==wkday ,]
prior.1437A$sp.prop<-with(aux.prop,ns(Flow.idx,knots=c(seq(60,120,by=5),seq(180,220,by=5))))

#--------------------------------------#
#          Spline prediction           #
#--------------------------------------#
spl.prop<-predict(prior.1437A$sp.prop,t5.1437A[c.chld,]$Flow.idx)
# removing missing (otherwise predict function won't work - obs won't be used anyway)
x.occ<-t5.1437A[c.chld,]$Occ_lag1
for (i in 1:length(x.occ)) { if (is.na(x.occ[i])==T) x.occ[i]=x.occ[i-1]}
spl.occ<-predict(prior.1437A$sp.occ,x.occ)
x.hway<-t5.1437A[c.chld & t5.1437A$Hway_lag1>1.3,]$Hway_lag1
for (i in 1:length(x.hway)) { if (is.na(x.hway[i])==T) x.hway[i]=x.hway[i-1]}
spl.hway<-predict(prior.1437A$sp.hway,x.hway)
x.spd<-t5.1437A[c.chld,]$Speed_lag1
for (i in 1:length(x.spd)) { if (is.na(x.spd[i])==T) x.spd[i]=x.spd[i-1]}
spl.spd<-predict(prior.1437A$sp.spd,x.spd)


#***************************#
#   1437A MODEL FITTING     #
#***************************#
#  S/S Model
MDM.1437A<-list()
VDF=rep(0.9,3)
occsp<-t5.1437A[c.chld,]$Occ_lag1*t5.1437A[c.chld,]$Speed_lag1
for (i in 1:length(occsp)) if (is.na(occsp[i])==T) occsp[i]<-occsp[i-1]
MDM.1437A[[1]]<-MDM.VLaw.v12(info=t5.1437A[c.chld,c("Flow","Date1","Hour")],nparents=1,beta.var=betas.ch.rnd[[4]],Vdf=VDF,
             X_t=t5.1431A[c.par,]$Flow*spl.prop[1:nrow(spl.prop),],F_t.star=spl.prop[1:nrow(spl.prop),],
             parent.f.mean=dlm.1431A[[1]]$f,parent.f.var=dlm.1431A[[1]]$Q,prior.obj=prior.1437A,n.test=n.test)
#   S/F Model
MDM.1437A[[2]]<-MDM.VLaw.v12(info=t5.1437A[c.chld,c("Flow","Date1","Hour")],nparents=1,beta.var=betas.ch.rnd[[4]],Vdf=VDF,
             X_t=t5.1431A[c.par,]$Flow*spl.prop[1:nrow(spl.prop),],F_t.star=spl.prop[1:nrow(spl.prop),],
             parent.f.mean=dlm.1431A[[9]]$f,parent.f.var=dlm.1431A[[9]]$Q,prior.obj=prior.1437A,n.test=n.test)
# F/S Model
MDM.1437A[[3]]<-MDM.VLaw.v12(info=t5.1437A[c.chld,c("Flow","Date1","Hour")],nparents=1,beta.var=betas.ch.rnd[[4]],Vdf=VDF,
             X_t=as.matrix(cbind(t5.1431A[c.par,]$Flow*spl.prop[1:nrow(spl.prop),],spl.hway[1:nrow(spl.hway),],
             spl.occ[1:nrow(spl.occ),],spl.spd[1:nrow(spl.spd),],occsp)),F_t.star=spl.prop[1:nrow(spl.prop),],
             parent.f.mean=dlm.1431A[[1]]$f,parent.f.var=dlm.1431A[[1]]$Q,prior.obj=prior.1437A,n.test=n.test)
# F/F Model
MDM.1437A[[4]]<-MDM.VLaw.v12(info=t5.1437A[c.chld,c("Flow","Date1","Hour")],nparents=1,beta.var=betas.ch.rnd[[4]],Vdf=VDF,
               X_t=as.matrix(cbind(t5.1431A[c.par,]$Flow*spl.prop[1:nrow(spl.prop),],spl.hway[1:nrow(spl.hway),], delta=1,
               spl.occ[1:nrow(spl.occ),],spl.spd[1:nrow(spl.spd),],occsp)),F_t.star=spl.prop[1:nrow(spl.prop),],
               parent.f.mean=dlm.1431A[[9]]$f,parent.f.var=dlm.1431A[[9]]$Q,prior.obj=prior.1437A,n.test=n.test)
# DLM F/F
sp.flow.1437A<-with(aux.prop,ns(Flow.idx,knots=c(seq(40,100,by=10),seq(180,220,by=10))))
spl.flow<-predict(sp.flow.1437A,t5.1437A[c.chld,]$Flow.idx)
MDM.1437A[[5]]<-MDM.VLaw.v12(info=t5.1437A[c.chld,c("Flow","Date1","Hour")],nparents=0,beta.var=betas.ch.rnd[[4]],Vdf=VDF,
     X_t=as.matrix(cbind(spl.flow[1:nrow(spl.flow),],spl.hway[1:nrow(spl.hway),],spl.occ[1:nrow(spl.occ),],
     spl.spd[1:nrow(spl.spd),],occsp)),prior.obj=prior.1437A,n.test=n.test)


r<--(1:n.test)
#-----------------------#
# Joint LPL 1431A/1437A:#
#-----------------------#
aux.l<-length(t5.1437A[c.chld,]$Date[r])
jlpl.1437A<-data.frame(Date=t5.1437A[c.chld,]$Date[r],Hour=t5.1437A[c.chld,]$Hour[r],lpl.DLM=rep(NA,aux.l),
                       lpl.m1=rep(NA,aux.l),lpl.m2=rep(NA,aux.l),lpl.m3=rep(NA,aux.l),lpl.m4=rep(NA,aux.l))

#"S/S"
jlpl.1437A$lpl.m1<-dnorm(t5.1431A[c.par,]$Flow[r],mean=dlm.1431A[[1]]$f[r],sd=sqrt(dlm.1431A[[1]]$Q[r]),log=TRUE)+
                   dnorm(t5.1437A[c.chld,]$Flow[r],mean=MDM.1437A[[1]]$f[r],sd=sqrt(MDM.1437A[[1]]$Q[r]),log=TRUE)
#"S/F"
jlpl.1437A$lpl.m2<-dnorm(t5.1431A[c.par,]$Flow[r],mean=dlm.1431A[[9]]$f[r],sd=sqrt(dlm.1431A[[9]]$Q[r]),log=TRUE)+
                   dnorm(t5.1437A[c.chld,]$Flow[r],mean=MDM.1437A[[2]]$f[r],sd=sqrt(MDM.1437A[[2]]$Q[r]),log=TRUE)
#"F/S"
jlpl.1437A$lpl.m3<-dnorm(t5.1431A[c.par,]$Flow[r],mean=dlm.1431A[[1]]$f[r],sd=sqrt(dlm.1431A[[1]]$Q[r]),log=TRUE)+
                   dnorm(t5.1437A[c.chld,]$Flow[r],mean=MDM.1437A[[3]]$f[r],sd=sqrt(MDM.1437A[[3]]$Q[r]),log=TRUE)
#"F/F"
jlpl.1437A$lpl.m4<-dnorm(t5.1431A[c.par,]$Flow[r],mean=dlm.1431A[[9]]$f[r],sd=sqrt(dlm.1431A[[9]]$Q[r]),log=TRUE)+
                   dnorm(t5.1437A[c.chld,]$Flow[r],mean=MDM.1437A[[4]]$f[r],sd=sqrt(MDM.1437A[[4]]$Q[r]),log=TRUE)

#"F/F DLM"
jlpl.1437A$lpl.DLM<-dnorm(t5.1431A[c.par,]$Flow[r],mean=dlm.1431A[[9]]$f[r],sd=sqrt(MDM.1437A[[5]]$Q[r]),log=TRUE)+
                   dnorm(t5.1437A[c.chld,]$Flow[r],mean=MDM.1437A[[5]]$f[r],sd=sqrt(MDM.1437A[[5]]$Q[r]),log=TRUE)

(jlpl.1437A<-data.frame(model=c("S/S","S/F","F/S","F/F","F/F DLM"),
                           lpl=c(sum(jlpl.1437A$lpl.m1,na.rm="T"),sum(jlpl.1437A$lpl.m2,na.rm="T"),
                                 sum(jlpl.1437A$lpl.m3,na.rm="T"),sum(jlpl.1437A$lpl.m4,na.rm="T"),
                                 sum(jlpl.1437A$lpl.DLM,na.rm="T"))))


save(MDM.1437A,file="/home/terrahome/oaj9/splines/MDM.1437A.Rdata")
save(prior.1437A,file="/home/terrahome/oaj9/splines/prior.1437A.Rdata")

#**************************************************************************************************#
#-------------------------------------1441A MODELLING----------------------------------------------#
#**************************************************************************************************#
#load("/home/terrahome/oaj9/splines/MDM.1437A.Rdata")
#-------------------------#
#     Data Selection      #
#-------------------------#
t5.1441A$Flow.idx<-rep(1:288,dim(table(t5.1441A$Date)))
# prior and learning data:
c.pri.gchld<-months(t5.1441A$Date1) %in% prior.months & !t5.1441A$Date %in% test.dates & t5.1441A$Hour %in% time.filter
c.gchld<-(months(t5.1441A$Date1) %in% learning.months | t5.1441A$Date %in% test.dates) & t5.1441A$Hour %in% time.filter  & t5.1441A$Dweek==wkday
period<-dim(table(t5.1441A[c.pri.gchld,]$Time5))
n.test<-period*length(test.dates)
with(t5.1441A[c.pri.gchld,],table(months(Date1)))
with(t5.1441A[c.gchld,],as.matrix(table(Date1)))
#----------#
#  priors  #
#----------#
#prior.1441A<-MDM.DLM.find.prior(t5.1441A[c.pri.gchld,]$Flow,
#                                     t5.1437A[c.pri.chld,]$Flow, period,0.001)
prior.1441A<-list()
prior.1441A$sigma2<-sum((t5.1441A[c.pri.gchld,]$Flow-t5.1437A[c.pri.chld,]$Flow)^2,na.rm=T)/
                  length(t5.1441A[c.pri.gchld,]$Flow)
prior.1441A$mean.prior<-with(t5.1441A[c.pri.gchld & t5.1441A$Time5==t5.1441A[c.pri.gchld,]$Time5[1],],
                      mean(Flow,na.rm=T))
# total flow per period (denominator of the vlaw):
prior.1441A$total.flow<-sum(t5.1441A[t5.1441A$Date =="04/10/2010",]$Flow,na.rm=T)

#--------------------------------------#
#          Spline fitting              #
#--------------------------------------#
c.gchld.spl<-c.pri.gchld & is.na(t5.1441A$Flow_lag1)==F  &  is.na(t5.1441A$Flow)==F & is.na(t5.1441A$Occ_lag1)==F &  is.na(t5.1441A$Hway_lag1)==F &  is.na(t5.1441A$Speed_lag1)==F
prior.1441A$sp.occ<-with(t5.1441A[c.gchld.spl ,],ns(Occ_lag1,df=4))
prior.1441A$sp.hway<-with(t5.1441A[c.gchld.spl & t5.1441A$Hway_lag1>1.3 ,],ns(Hway_lag1,df=4))
prior.1441A$sp.spd<-with(t5.1441A[c.gchld.spl,],ns(Speed_lag1,df=4))
aux.prop<-data.frame(Flow.idx=t5.1441A[c.pri.gchld,]$Flow.idx,Dweek=t5.1441A[c.pri.gchld,]$Dweek,
                     Flow.par=t5.1437A[c.pri.gchld,]$Flow,Flow.chd=t5.1441A[c.pri.gchld,]$Flow,
                     Flow.prop=t5.1441A[c.pri.gchld,]$Flow/t5.1437A[c.pri.chld,]$Flow)
aux.prop<-aux.prop[is.na(aux.prop$Flow.prop)==F & aux.prop$Dweek==wkday ,]
prior.1441A$sp.prop<-with(aux.prop,ns(Flow.idx,knots=c(seq(60,120,by=5),seq(180,220,by=5))))
#--------------------------------------#
#          Spline prediction           #
#--------------------------------------#
spl.prop<-predict(prior.1441A$sp.prop,t5.1441A[c.gchld,]$Flow.idx)
# removing missing (otherwise predict function won't work - obs won't be used anyway)
x.occ<-t5.1441A[c.gchld,]$Occ_lag1
for (i in 1:length(x.occ)) { if (is.na(x.occ[i])==T) x.occ[i]=x.occ[i-1]}
spl.occ<-predict(prior.1441A$sp.occ,x.occ)

x.hway<-t5.1441A[c.gchld,]$Hway_lag1
for (i in 1:length(x.hway)) { if (is.na(x.hway[i])==T) x.hway[i]=x.hway[i-1]}
spl.hway<-predict(prior.1441A$sp.hway,x.hway)

x.spd<-t5.1441A[c.gchld ,]$Speed_lag1
for (i in 1:length(x.spd)) { if (is.na(x.spd[i])==T) x.spd[i]=x.spd[i-1]}
spl.spd<-predict(prior.1441A$sp.spd,x.spd)


#***************************#
#   1437A MODEL FITTING     #
#***************************#
#S/S/S
VDF=rep(0.85,3)
MDM.1441A<-list()
occsp<-t5.1441A[c.gchld,]$Occ_lag1*t5.1441A[c.gchld,]$Speed_lag1
for (i in 1:length(occsp)) if (is.na(occsp[i])==T) occsp[i]<-occsp[i-1]
#source("/home/terrahome/oaj9/_codes/VLaw_Function-with_DLM_pack.r")
#S/S/S
MDM.1441A<-list()
MDM.1441A[[1]]<-MDM.VLaw.v12(info=t5.1441A[c.gchld,c("Flow","Date1","Hour")],nparents=1,beta.var=betas.ch.rnd[[4]],
            X_t=t5.1437A[c.chld,]$Flow*spl.prop[1:nrow(spl.prop),],prior.obj=prior.1441A,n.test=n.test,Vdf=VDF,
            parent.f.mean=MDM.1437A[[1]]$f.marg,parent.f.var=MDM.1437A[[1]]$marg.var,
            F_t.star=spl.prop[1:nrow(spl.prop),])

#S/S/F
MDM.1441A[[2]]<-MDM.VLaw.v12(info=t5.1441A[c.gchld,c("Flow","Date1","Hour")],nparents=1,beta.var=betas.ch.rnd[[4]],
            X_t=t5.1437A[c.chld,]$Flow*spl.prop[1:nrow(spl.prop),],prior.obj=prior.1441A,n.test=n.test,Vdf=VDF,
            parent.f.mean=MDM.1437A[[2]]$f.marg,parent.f.var=MDM.1437A[[2]]$marg.var,
            F_t.star=spl.prop[1:nrow(spl.prop),])
#S/F/F
MDM.1441A[[3]]<-MDM.VLaw.v12(info=t5.1441A[c.gchld,c("Flow","Date1","Hour")],nparents=1,beta.var=betas.ch.rnd[[4]],
            X_t=t5.1437A[c.chld,]$Flow*spl.prop[1:nrow(spl.prop),],prior.obj=prior.1441A,n.test=n.test,Vdf=VDF,
            parent.f.mean=MDM.1437A[[4]]$f.marg,parent.f.var=MDM.1437A[[4]]$marg.var,
            F_t.star=spl.prop[1:nrow(spl.prop),])
# F/F/F
MDM.1441A[[4]]<-MDM.VLaw.v12(info=t5.1441A[c.gchld,c("Flow","Date1","Hour")],nparents=1,beta.var=betas.ch.rnd[[4]],
             X_t=as.matrix(cbind(t5.1437A[c.chld,]$Flow*spl.prop[1:nrow(spl.prop),],spl.hway[1:nrow(spl.hway),],
             spl.occ[1:nrow(spl.occ),],spl.spd[1:nrow(spl.spd),])),prior.obj=prior.1441A,n.test=n.test,Vdf=VDF,
             parent.f.mean=MDM.1437A[[4]]$f.marg,parent.f.var=MDM.1437A[[4]]$marg.var,
             F_t.star=spl.prop[1:nrow(spl.prop),])
#single DLM with F:
sp.flow.1441A<-with(aux.prop,ns(Flow.idx,knots=c(seq(40,100,by=10),seq(180,220,by=10))))
spl.flow<-predict(sp.flow.1441A,t5.1441A[c.gchld,]$Flow.idx)
MDM.1441A[[5]]<-MDM.VLaw.v12(info=t5.1441A[c.gchld,c("Flow","Date1","Hour")],nparents=0,beta.var=betas.ch.rnd[[4]],
   X_t=as.matrix(cbind(spl.flow[1:nrow(spl.flow),],spl.hway[1:nrow(spl.hway),],spl.occ[1:nrow(spl.occ),],
  spl.spd[1:nrow(spl.spd),])),Vdf=VDF,prior.obj=prior.1437A,n.test=n.test)


r<--(1:n.test)
# Joint LPL 1431A/1437A:
aux.l<-length(t5.1437A[c.gchld,]$Date[r])
jlpl.1441A<-data.frame(Date=t5.1437A[c.chld,]$Date[r],Hour=t5.1437A[c.chld,]$Hour[r],lpl.DLM=rep(NA,aux.l),
                       lpl.m1=rep(NA,aux.l),lpl.m2=rep(NA,aux.l),lpl.m3=rep(NA,aux.l),lpl.m4=rep(NA,aux.l))
#"S/S/S"
jlpl.1441A$lpl.m1<-dnorm(t5.1431A[c.par,]$Flow[r],mean=dlm.1431A[[1]]$f[r],sd=sqrt(dlm.1431A[[1]]$Q[r]),log=TRUE)+
                   dnorm(t5.1437A[c.chld,]$Flow[r],mean=MDM.1437A[[1]]$f[r],sd=sqrt(MDM.1437A[[1]]$Q[r]),log=TRUE)+
                   dnorm(t5.1441A[c.gchld,]$Flow[r],mean=MDM.1441A[[1]]$f[r],sd=sqrt(MDM.1441A[[1]]$Q[r]),log=TRUE)
# "S/S/F"
jlpl.1441A$lpl.m2<-dnorm(t5.1431A[c.par,]$Flow[r],mean=dlm.1431A[[9]]$f[r],sd=sqrt(dlm.1431A[[9]]$Q[r]),log=TRUE)+
                   dnorm(t5.1437A[c.chld,]$Flow[r],mean=MDM.1437A[[2]]$f[r],sd=sqrt(MDM.1437A[[2]]$Q[r]),log=TRUE)+
                   dnorm(t5.1441A[c.gchld,]$Flow[r],mean=MDM.1441A[[2]]$f[r],sd=sqrt(MDM.1441A[[2]]$Q[r]),log=TRUE)
#"S/F/F"
jlpl.1441A$lpl.m3<-dnorm(t5.1431A[c.par,]$Flow[r],mean=dlm.1431A[[9]]$f[r],sd=sqrt(dlm.1431A[[9]]$Q[r]),log=TRUE)+
                   dnorm(t5.1437A[c.chld,]$Flow[r],mean=MDM.1437A[[4]]$f[r],sd=sqrt(MDM.1437A[[4]]$Q[r]),log=TRUE)+
                   dnorm(t5.1441A[c.gchld,]$Flow[r],mean=MDM.1441A[[3]]$f[r],sd=sqrt(MDM.1441A[[3]]$Q[r]),log=TRUE)
#"F/F/F"
jlpl.1441A$lpl.m4<-dnorm(t5.1431A[c.par,]$Flow[r],mean=dlm.1431A[[9]]$f[r],sd=sqrt(dlm.1431A[[9]]$Q[r]),log=TRUE)+
                   dnorm(t5.1437A[c.chld,]$Flow[r],mean=MDM.1437A[[4]]$f[r],sd=sqrt(MDM.1437A[[4]]$Q[r]),log=TRUE)+
                   dnorm(t5.1441A[c.gchld,]$Flow[r],mean=MDM.1441A[[4]]$f[r],sd=sqrt(MDM.1441A[[4]]$Q[r]),log=TRUE)

#"F/F/F DLM"
jlpl.1441A$lpl.DLM<-dnorm(t5.1431A[c.par,]$Flow[r],mean=dlm.1431A[[9]]$f[r],sd=sqrt(dlm.1431A[[9]]$Q[r]),log=TRUE)+
                   dnorm(t5.1437A[c.chld,]$Flow[r],mean=MDM.1437A[[5]]$f[r],sd=sqrt(MDM.1437A[[5]]$Q[r]),log=TRUE)+
                   dnorm(t5.1441A[c.gchld,]$Flow[r],mean=MDM.1441A[[5]]$f[r],sd=sqrt(MDM.1441A[[5]]$Q[r]),log=TRUE)

(jlpl.1441A<-data.frame(model=c("S/S/S","S/S/F","S/F/F","F/F/F","F/F/F DLM"),
                           lpl=c(sum(jlpl.1441A$lpl.m1,na.rm="T"),sum(jlpl.1441A$lpl.m2,na.rm="T"),
                                 sum(jlpl.1441A$lpl.m3,na.rm="T"),sum(jlpl.1441A$lpl.m4,na.rm="T"),
                                 sum(jlpl.1441A$lpl.DLM,na.rm="T"))))



save(MDM.1441A,file="/home/terrahome/oaj9/splines/MDM.1441A.Rdata")
save(prior.1441A,file="/home/terrahome/oaj9/splines/prior.1441A.Rdata")







