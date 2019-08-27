#+++++++++++++++++++++++++++++++++++++++++++++++++++++#
#     Spline MDM                                      #
#     DLM for 6013B                                   #
#     Set/2011                                        #
#+++++++++++++++++++++++++++++++++++++++++++++++++++++#
library(chron);library(splines);library(car)
source("/home/terrahome/oaj9/_codes/MDM-VLaw-v12.r")
load("/home/terrahome/oaj9/_data/all2010/t5.6013B.Rdata")
load("/home/terrahome/oaj9/_data/all2010/t5.6007L.Rdata")
load("/home/terrahome/oaj9/_data/all2010/t5.6004L.Rdata")
load("/home/terrahome/oaj9/_data/all2010/stats/betas.9200B.Rdata")
load("/home/terrahome/oaj9/_data/all2010/stats/betas.9206B.Rdata")

time.filter<-6:20
test.dates<-c("18/08/2010","25/08/2010")
prior.months<- c("Jul","Aug")
#learning.months<-c("Sep","Oct","Nov")
learning.months<-c("Sep","Oct")
wkday<-"Wed"
#**************************************************************************************************#
#-------------------------------------6013B MODELLING----------------------------------------------#
#**************************************************************************************************#

#-------------------------#
#     Data Selection      #
#-------------------------#
# prior and learning data:
t5.6013B$Flow.idx<-rep(1:288,dim(table(t5.6013B$Date)))
c.pri.par<-months(t5.6013B$Date1) %in% prior.months & !t5.6013B$Date %in% test.dates & t5.6013B$Hour %in% time.filter
c.par<-(months(t5.6013B$Date1) %in% learning.months | t5.6013B$Date %in% test.dates)  & t5.6013B$Hour %in% time.filter  & t5.6013B$Dweek==wkday
period<-dim(table(t5.6013B[c.pri.par,]$Time5))
t5.6013B[is.na(t5.6013B$Hway_lag1)==T,]$Flow<-NA
t5.6013B[is.na(t5.6013B$Occ_lag1)==T,]$Flow<-NA
n.test<-period*length(test.dates)
with(t5.6013B[c.pri.par,],table(months(Date1)))
with(t5.6013B[c.par,],as.matrix(table(Date1)))



#----------#
#  priors  #
#----------#
prior.6013B<-list()
prior.6013B$sigma2<-with(t5.6013B[c.pri.par,], prior.5min(Flow,t5.6013B[c.pri.par,],dim(table(Date1)),period)$S0)
prior.6013B$mean.prior<-with(t5.6013B[c.pri.par & t5.6013B$Time5==t5.6013B[c.pri.par,]$Time5[1],],
                      mean(Flow,na.rm=T))
# total flow per period (denominator of the vlaw):
prior.6013B$total.flow<-sum(t5.6013B[t5.6013B$Date =="04/10/2010",]$Flow,na.rm=T)
#--------------------------------------#
#          Spline fitting              #
#--------------------------------------#
c.spl.par<-c.pri.par  & is.na(t5.6013B$Flow)==F & is.na(t5.6013B$Occ_lag1)==F  & is.na(t5.6013B$Hway_lag1)==F & is.na(t5.6013B$Speed_lag1)==F
c.spl.flow<-c.pri.par &  t5.6013B$Dweek==wkday & is.na(t5.6013B$Flow)==F
with(t5.6013B[c.spl.par,],table(months(Date1)))
with(t5.6013B[c.spl.flow,],table(months(Date1)))
prior.6013B$sp.flow<-with(t5.6013B[c.spl.flow,],ns(Flow.idx,knots=c(seq(60,130,by=5),seq(180,220,by=10))))
prior.6013B$sp.hway<-with(t5.6013B[c.spl.par & t5.6013B$Hway_lag1>1,],ns(Hway_lag1,df=4))
prior.6013B$sp.spd<-with(t5.6013B[c.spl.par,],ns(Speed_lag1,df=4))
prior.6013B$sp.occ<-with(t5.6013B[c.spl.par,],ns(Occ_lag1,df=4))

#--------------------------------------#
#          Spline prediction           #
#--------------------------------------#
# removing missing (otherwise predict function won't work - obs won't be used anyway)

x.flow<-t5.6013B[c.par,]$Flow.idx
for (i in 1:length(x.flow)) { if (is.na(x.flow[i])==T) x.flow[i]=x.flow[i-1]}
spl.flow<-predict(prior.6013B$sp.flow,x.flow)

x.occ<-t5.6013B[c.par,]$Occ_lag1
for (i in 1:length(x.occ)) { if (is.na(x.occ[i])==T) x.occ[i]=x.occ[i-1]}
spl.occ<-predict(prior.6013B$sp.occ,x.occ)

#x.hway<-t5.6013B[c.par & t5.6013B$Hway_lag1>1,]$Hway_lag1
x.hway<-t5.6013B[c.par,]$Hway_lag1
for (i in 1:length(x.hway)) { if (is.na(x.hway[i])==T) x.hway[i]=x.hway[i-1]}
spl.hway<-predict(prior.6013B$sp.hway,x.hway)

x.spd<-t5.6013B[c.par ,]$Speed_lag1
for (i in 1:length(x.spd)) { if (is.na(x.spd[i])==T) x.spd[i]=x.spd[i-1]}
spl.spd<-predict(prior.6013B$sp.spd,x.spd)


#***************************#
#   6013B MODEL FITTING     #
#***************************#
dlm.6013B<-list()
VDF=rep(0.99,3)

#   S:  Flow spline trend
dlm.6013B[[1]]<-MDM.VLaw.v12(info=t5.6013B[c.par,c("Flow","Date1","Hour")],nparents=0,beta.var=betas.rnd[[2]],Vdf=VDF,
                              X_t=as.matrix(spl.flow[1:nrow(spl.flow),]),prior.obj=prior.6013B,n.test=n.test)
#---------individual models----------------------------------#
#   SO  (S + occupancy)
dlm.6013B[[2]]<-MDM.VLaw.v12(info=t5.6013B[c.par,c("Flow","Date1","Hour")],,nparents=0,beta.var=betas.rnd[[2]],Vdf=VDF,
                             X_t=as.matrix(cbind(spl.flow[1:nrow(spl.flow),],spl.occ[1:nrow(spl.occ),])),
                             prior.obj=prior.6013B,n.test=n.test)
#   SH (S + headway)
dlm.6013B[[3]]<-MDM.VLaw.v12(info=t5.6013B[c.par,c("Flow","Date1","Hour")],,nparents=0,beta.var=betas.rnd[[2]],Vdf=VDF,
                             X_t=as.matrix(cbind(spl.flow[1:nrow(spl.flow),],spl.hway[1:nrow(spl.hway),])),
                             prior.obj=prior.6013B,n.test=n.test)
#   SHinv  (S + 1/headway)
dlm.6013B[[4]]<-MDM.VLaw.v12(info=t5.6013B[c.par,c("Flow","Date1","Hour")],,nparents=0,beta.var=betas.rnd[[2]],Vdf=VDF,
                             X_t=as.matrix(cbind(spl.flow[1:nrow(spl.flow),],1/(t5.6013B[c.par,]$Hway_lag1+0.1))),
                             prior.obj=prior.6013B,n.test=n.test)
# - SSp (S + Speed)
dlm.6013B[[5]]<-MDM.VLaw.v12(info=t5.6013B[c.par,c("Flow","Date1","Hour")],nparents=0,beta.var=betas.rnd[[2]],Vdf=VDF,
                             X_t=as.matrix(cbind(spl.flow[1:nrow(spl.flow),],spl.spd[1:nrow(spl.spd),])),
                             prior.obj=prior.6013B,n.test=n.test)
#------Interaction test---------------------------------------#
# - SOSp: (S + occupancy + speed)
dlm.6013B[[6]]<-MDM.VLaw.v12(info=t5.6013B[c.par,c("Flow","Date1","Hour")],nparents=0,beta.var=betas.rnd[[2]],Vdf=VDF,
                             X_t=as.matrix(cbind(spl.flow[1:nrow(spl.flow),],spl.occ[1:nrow(spl.occ),],spl.spd[1:nrow(spl.spd),])),
                             prior.obj=prior.6013B,n.test=n.test)
# - S(OSp) : (S + occupancy + speed + occ x speed)
dlm.6013B[[7]]<-MDM.VLaw.v12(info=t5.6013B[c.par,c("Flow","Date1","Hour")],nparents=0,beta.var=betas.rnd[[2]],Vdf=VDF,
                             X_t=as.matrix(cbind(spl.flow[1:nrow(spl.flow),],spl.occ[1:nrow(spl.occ),],spl.spd[1:nrow(spl.spd),],
                                                  t5.6013B[c.par,]$Occ_lag1*t5.6013B[c.par,]$Speed_lag1)),
                             prior.obj=prior.6013B,n.test=n.test)
#------models with all variables:------------------------------#
# - Spline Flow Level + O + H + Sp
dlm.6013B[[8]]<-MDM.VLaw.v12(info=t5.6013B[c.par,c("Flow","Date1","Hour")],nparents=0,beta.var=betas.rnd[[2]],Vdf=VDF,
                             X_t=as.matrix(cbind(spl.flow[1:nrow(spl.flow),],spl.occ[1:nrow(spl.occ),],
                                                 spl.hway[1:nrow(spl.hway),],spl.spd[1:nrow(spl.spd),])),
                             prior.obj=prior.6013B,n.test=n.test)
# - Spline Flow Level +  Sp + Occ + (OSp) +  Hway
dlm.6013B[[9]]<-MDM.VLaw.v12(info=t5.6013B[c.par,c("Flow","Date1","Hour")],nparents=0,beta.var=betas.rnd[[2]],Vdf=VDF,
                             X_t=as.matrix(cbind(spl.flow[1:nrow(spl.flow),],spl.hway[1:nrow(spl.hway),],spl.occ[1:nrow(spl.occ),],
                                                  spl.spd[1:nrow(spl.spd),],t5.6013B[c.par,]$Occ_lag1*t5.6013B[c.par,]$Speed_lag1)),
                             prior.obj=prior.6013B,n.test=n.test)



(res2.6013B<-data.frame(Model=c("S","SO","SH","SHinv","SSp","SOSp","S(OSp)","SOSpH","S(OSp)H"),
            LPL=c(dlm.6013B[[1]]$perf[5,2],dlm.6013B[[2]]$perf[5,2], dlm.6013B[[3]]$perf[5,2],
                  dlm.6013B[[4]]$perf[5,2],dlm.6013B[[5]]$perf[5,2], dlm.6013B[[6]]$perf[5,2],
                  dlm.6013B[[7]]$perf[5,2],dlm.6013B[[8]]$perf[5,2], dlm.6013B[[9]]$perf[5,2])))


save(dlm.6013B,file="/home/terrahome/oaj9/splines/dlm.6013B.Rdata")
save(prior.6013B,file="/home/terrahome/oaj9/splines/prior.6013B.Rdata")


#**************************************************************************************************#
#-------------------------------------6007L MODELLING----------------------------------------------#
#**************************************************************************************************#
#load("/home/terrahome/oaj9/splines/dlm.6013B.Rdata")
#-------------------------#
#     Data Selection      #
#-------------------------#
t5.6007L$Flow.idx<-rep(1:288,dim(table(t5.6007L$Date)))
# prior and learning data:
period<-dim(table(t5.6007L[c.pri.chld,]$Time5))
c.pri.chld<-months(t5.6007L$Date1) %in% prior.months & !t5.6007L$Date %in% test.dates & t5.6007L$Hour %in% time.filter
c.chld<-(months(t5.6007L$Date1) %in% learning.months | t5.6007L$Date %in% test.dates) & t5.6007L$Hour %in% time.filter  & t5.6007L$Dweek==wkday
n.test<-period*length(test.dates)

with(t5.6007L[c.pri.chld,],table(months(Date1)))
with(t5.6007L[c.chld,],as.matrix(table(Date1)))

#----------#
#  priors  #
#----------#
prior.6007L<-MDM.DLM.find.prior(t5.6007L[c.pri.chld,]$Flow,t5.6013B[c.pri.par,]$Flow,period,0.001)
prior.6007L$sigma2<-sum((t5.6007L[c.pri.chld,]$Flow-t5.6013B[c.pri.par,]$Flow)^2,na.rm=T)/
                         length(t5.6007L[c.pri.chld,]$Flow)
prior.6007L$mean.prior<-with(t5.6007L[c.pri.chld & t5.6007L$Time5==t5.6007L[c.pri.chld,]$Time5[1],],
                      mean(Flow,na.rm=T))

# total flow per period (denominator of the vlaw):
prior.6007L$total.flow<-sum(t5.6007L[t5.6007L$Date =="04/10/2010",]$Flow,na.rm=T)

#
#--------------------------------------#
#          Spline fitting              #
#--------------------------------------#
c.chld.spl<-c.pri.chld & is.na(t5.6007L$Flow)==F & is.na(t5.6007L$Flow)==F & is.na(t5.6007L$Occ_lag1)==F & is.na(t5.6007L$Speed_lag1)==F & is.na(t5.6007L$Hway_lag1)==F
prior.6007L$sp.occ<-with(t5.6007L[c.chld.spl ,],ns(Occ_lag1,df=4))
prior.6007L$sp.hway<-with(t5.6007L[c.chld.spl & t5.6007L$Hway_lag1>1.3 ,],ns(Hway_lag1,df=4))
prior.6007L$sp.spd<-with(t5.6007L[c.chld.spl,],ns(Speed_lag1,df=4))

aux.prop<-data.frame(Flow.idx=t5.6007L[c.pri.chld,]$Flow.idx,Dweek=t5.6007L[c.pri.chld,]$Dweek,
                     Flow.par=t5.6013B[c.pri.par,]$Flow,Flow.chd=t5.6007L[c.pri.chld,]$Flow,
                     Flow.prop=t5.6007L[c.pri.chld,]$Flow/t5.6013B[c.pri.par,]$Flow)
aux.prop<-aux.prop[is.na(aux.prop$Flow.prop)==F & aux.prop$Dweek==wkday ,]
prior.6007L$sp.prop<-with(aux.prop,ns(Flow.idx,knots=c(seq(60,120,by=10),seq(180,220,by=10))))

#--------------------------------------#
#          Spline prediction           #
#--------------------------------------#
spl.prop<-predict(prior.6007L$sp.prop,t5.6007L[c.chld,]$Flow.idx)
# removing missing (otherwise predict function won't work - obs won't be used anyway)
x.occ<-t5.6007L[c.chld,]$Occ_lag1
for (i in 1:length(x.occ)) { if (is.na(x.occ[i])==T) x.occ[i]=x.occ[i-1]}
spl.occ<-predict(prior.6007L$sp.occ,x.occ)
x.hway<-t5.6007L[c.chld & t5.6007L$Hway_lag1>1.3,]$Hway_lag1
for (i in 1:length(x.hway)) { if (is.na(x.hway[i])==T) x.hway[i]=x.hway[i-1]}
spl.hway<-predict(prior.6007L$sp.hway,x.hway)
x.spd<-t5.6007L[c.chld,]$Speed_lag1
for (i in 1:length(x.spd)) { if (is.na(x.spd[i])==T) x.spd[i]=x.spd[i-1]}
spl.spd<-predict(prior.6007L$sp.spd,x.spd)



#***************************#
#   6007L MODEL FITTING     #
#***************************#
#  S/S Model
MDM.6007L<-list()
occsp<-t5.6007L[c.chld,]$Occ_lag1*t5.6007L[c.chld,]$Speed_lag1
for (i in 1:length(occsp)) if (is.na(occsp[i])==T) occsp[i]<-occsp[i-1]
VDF=rep(0.88,3)
MDM.6007L[[1]]<-MDM.VLaw.v12(info=t5.6007L[c.chld,c("Flow","Date1","Hour")],nparents=1,beta.var=betas.ch.rnd[[2]],Vdf=VDF,
             X_t=t5.6013B[c.par,]$Flow*spl.prop[1:nrow(spl.prop),],F_t.star=spl.prop[1:nrow(spl.prop),],
             parent.f.mean=dlm.6013B[[1]]$f,parent.f.var=dlm.6013B[[1]]$Q,prior.obj=prior.6007L,n.test=n.test)
# S/F Model
MDM.6007L[[2]]<-MDM.VLaw.v12(info=t5.6007L[c.chld,c("Flow","Date1","Hour")],nparents=1,beta.var=betas.ch.rnd[[2]],Vdf=VDF,
             X_t=t5.6013B[c.par,]$Flow*spl.prop[1:nrow(spl.prop),],F_t.star=spl.prop[1:nrow(spl.prop),],
             parent.f.mean=dlm.6013B[[9]]$f,parent.f.var=dlm.6013B[[9]]$Q,prior.obj=prior.6007L,n.test=n.test)
# F/S Model
MDM.6007L[[3]]<-MDM.VLaw.v12(info=t5.6007L[c.chld,c("Flow","Date1","Hour")],nparents=1,beta.var=betas.ch.rnd[[2]],Vdf=VDF,
             X_t=as.matrix(cbind(t5.6013B[c.par,]$Flow*spl.prop[1:nrow(spl.prop),],spl.hway[1:nrow(spl.hway),],
             spl.occ[1:nrow(spl.occ),],spl.spd[1:nrow(spl.spd),],occsp)),
             F_t.star=spl.prop[1:nrow(spl.prop),],parent.f.mean=dlm.6013B[[1]]$f,parent.f.var=dlm.6013B[[1]]$Q,prior.obj=prior.6007L,n.test=n.test)
# F/F Model
MDM.6007L[[4]]<-MDM.VLaw.v12(info=t5.6007L[c.chld,c("Flow","Date1","Hour")],nparents=1,beta.var=betas.ch.rnd[[2]],Vdf=VDF,
               X_t=as.matrix(cbind(t5.6013B[c.par,]$Flow*spl.prop[1:nrow(spl.prop),],spl.hway[1:nrow(spl.hway),],
               spl.occ[1:nrow(spl.occ),],spl.spd[1:nrow(spl.spd),],occsp)),
               F_t.star=spl.prop[1:nrow(spl.prop),],parent.f.mean=dlm.6013B[[9]]$f,parent.f.var=dlm.6013B[[9]]$Q,prior.obj=prior.6007L,n.test=n.test)
# DLM F/F
sp.flow.6007L<-with(aux.prop,ns(Flow.idx,knots=c(seq(40,100,by=10),seq(180,220,by=10))))
spl.flow<-predict(sp.flow.6007L,t5.6007L[c.chld,]$Flow.idx)
MDM.6007L[[5]]<-MDM.VLaw.v12(info=t5.6007L[c.chld,c("Flow","Date1","Hour")],nparents=0,beta.var=betas.ch.rnd[[2]], Vdf=VDF,
     X_t=as.matrix(cbind(spl.flow[1:nrow(spl.flow),],spl.hway[1:nrow(spl.hway),],spl.occ[1:nrow(spl.occ),],
     spl.spd[1:nrow(spl.spd),],occsp)),
     prior.obj=prior.6007L,n.test=n.test)

r<--(1:n.test)
#-----------------------#
# Joint LPL 6013B/6007L:#
#-----------------------#
aux.l<-length(t5.6007L[c.chld,]$Date[r])
jlpl.6007L<-data.frame(Date=t5.6007L[c.chld,]$Date[r],Hour=t5.6007L[c.chld,]$Hour[r],lpl.DLM=rep(NA,aux.l),
                       lpl.m1=rep(NA,aux.l),lpl.m2=rep(NA,aux.l),lpl.m3=rep(NA,aux.l),lpl.m4=rep(NA,aux.l))

#"S/S"
jlpl.6007L$lpl.m1<-dnorm(t5.6013B[c.par,]$Flow[r],mean=dlm.6013B[[1]]$f[r],sd=sqrt(dlm.6013B[[1]]$Q[r]),log=TRUE)+
                   dnorm(t5.6007L[c.chld,]$Flow[r],mean=MDM.6007L[[1]]$f[r],sd=sqrt(MDM.6007L[[1]]$Q[r]),log=TRUE)
#"S/F"
jlpl.6007L$lpl.m2<-dnorm(t5.6013B[c.par,]$Flow[r],mean=dlm.6013B[[9]]$f[r],sd=sqrt(dlm.6013B[[9]]$Q[r]),log=TRUE)+
                   dnorm(t5.6007L[c.chld,]$Flow[r],mean=MDM.6007L[[2]]$f[r],sd=sqrt(MDM.6007L[[2]]$Q[r]),log=TRUE)
#"F/S"
jlpl.6007L$lpl.m3<-dnorm(t5.6013B[c.par,]$Flow[r],mean=dlm.6013B[[1]]$f[r],sd=sqrt(dlm.6013B[[1]]$Q[r]),log=TRUE)+
                   dnorm(t5.6007L[c.chld,]$Flow[r],mean=MDM.6007L[[3]]$f[r],sd=sqrt(MDM.6007L[[3]]$Q[r]),log=TRUE)
#"F/F"
jlpl.6007L$lpl.m4<-dnorm(t5.6013B[c.par,]$Flow[r],mean=dlm.6013B[[9]]$f[r],sd=sqrt(dlm.6013B[[9]]$Q[r]),log=TRUE)+
                   dnorm(t5.6007L[c.chld,]$Flow[r],mean=MDM.6007L[[4]]$f[r],sd=sqrt(MDM.6007L[[4]]$Q[r]),log=TRUE)

#"F/F DLM"
jlpl.6007L$lpl.DLM<-dnorm(t5.6013B[c.par,]$Flow[r],mean=dlm.6013B[[9]]$f[r],sd=sqrt(dlm.6013B[[9]]$Q[r]),log=TRUE)+
                   dnorm(t5.6007L[c.chld,]$Flow[r],mean=MDM.6007L[[5]]$f[r],sd=sqrt(MDM.6007L[[5]]$Q[r]),log=TRUE)

(jlpl.6007L<-data.frame(model=c("S/S","S/F","F/S","F/F","F/F DLM"),
                           lpl=c(sum(jlpl.6007L$lpl.m1,na.rm="T"),sum(jlpl.6007L$lpl.m2,na.rm="T"),
                                 sum(jlpl.6007L$lpl.m3,na.rm="T"),sum(jlpl.6007L$lpl.m4,na.rm="T"),
                                 sum(jlpl.6007L$lpl.DLM,na.rm="T"))))
                                 


save(MDM.6007L,file="/home/terrahome/oaj9/splines/MDM.6007L.Rdata")
save(prior.6007L,file="/home/terrahome/oaj9/splines/prior.6007L.Rdata")

#**************************************************************************************************#
#-------------------------------------6004L MODELLING----------------------------------------------#
#**************************************************************************************************#
#load("/home/terrahome/oaj9/splines/MDM.6007L.Rdata")
#-------------------------#
#     Data Selection      #
#-------------------------#
t5.6004L$Flow.idx<-rep(1:288,dim(table(t5.6004L$Date)))
# prior and learning data:
c.pri.gchld<-months(t5.6004L$Date1) %in% prior.months & !t5.6004L$Date %in% test.dates & t5.6004L$Hour %in% time.filter
c.gchld<-(months(t5.6004L$Date1) %in% learning.months | t5.6004L$Date %in% test.dates) & t5.6004L$Hour %in% time.filter  & t5.6004L$Dweek==wkday
period<-dim(table(t5.6004L[c.pri.gchld,]$Time5))
n.test<-period*length(test.dates)
with(t5.6004L[c.pri.gchld,],table(months(Date1)))
with(t5.6004L[c.gchld,],as.matrix(table(Date1)))
#----------#
#  priors  #
#----------#
#prior.6004L<-MDM.DLM.find.prior(t5.6004L[c.pri.gchld,]$Flow,
#                                     t5.6007L[c.pri.chld,]$Flow, period,0.001)
prior.6004L<-list()
prior.6004L$sigma2<-sum((t5.6004L[c.pri.gchld,]$Flow-t5.6007L[c.pri.chld,]$Flow)^2,na.rm=T)/
                  length(t5.6004L[c.pri.gchld,]$Flow)
prior.6004L$mean.prior<-with(t5.6004L[c.pri.gchld & t5.6004L$Time5==t5.6004L[c.pri.gchld,]$Time5[1],],
                      mean(Flow,na.rm=T))
# total flow per period (denominator of the vlaw):
prior.6004L$total.flow<-sum(t5.6004L[t5.6004L$Date =="04/10/2010",]$Flow,na.rm=T)

#--------------------------------------#
#          Spline fitting              #
#--------------------------------------#
c.gchld.spl<-c.pri.gchld & is.na(t5.6004L$Flow_lag1)==F  &  is.na(t5.6004L$Flow)==F & is.na(t5.6004L$Occ_lag1)==F &  is.na(t5.6004L$Hway_lag1)==F &  is.na(t5.6004L$Speed_lag1)==F
prior.6004L$sp.occ<-with(t5.6004L[c.gchld.spl ,],ns(Occ_lag1,df=4))
prior.6004L$sp.hway<-with(t5.6004L[c.gchld.spl & t5.6004L$Hway_lag1>1.3 ,],ns(Hway_lag1,df=4))
prior.6004L$sp.spd<-with(t5.6004L[c.gchld.spl,],ns(Speed_lag1,df=4))
aux.prop<-data.frame(Flow.idx=t5.6004L[c.pri.gchld,]$Flow.idx,Dweek=t5.6004L[c.pri.gchld,]$Dweek,
                     Flow.par=t5.6007L[c.pri.gchld,]$Flow,Flow.chd=t5.6004L[c.pri.gchld,]$Flow,
                     Flow.prop=t5.6004L[c.pri.gchld,]$Flow/t5.6007L[c.pri.chld,]$Flow)
aux.prop<-aux.prop[is.na(aux.prop$Flow.prop)==F & aux.prop$Dweek==wkday ,]
prior.6004L$sp.prop<-with(aux.prop,ns(Flow.idx,knots=c(seq(60,120,by=5),seq(180,220,by=5))))
#--------------------------------------#
#          Spline prediction           #
#--------------------------------------#
spl.prop<-predict(prior.6004L$sp.prop,t5.6004L[c.gchld,]$Flow.idx)
# removing missing (otherwise predict function won't work - obs won't be used anyway)
x.occ<-t5.6004L[c.gchld,]$Occ_lag1
for (i in 1:length(x.occ)) { if (is.na(x.occ[i])==T) x.occ[i]=x.occ[i-1]}
spl.occ<-predict(prior.6004L$sp.occ,x.occ)
t5.6004L[c.gchld,]$Hway_lag1[1]<-t5.6004L[c.gchld,]$Hway_lag1[2]
x.hway<-t5.6004L[c.gchld,]$Hway_lag1
for (i in 1:length(x.hway)) { if (is.na(x.hway[i])==T) x.hway[i]=x.hway[i-1]}
spl.hway<-predict(prior.6004L$sp.hway,x.hway)
x.spd<-t5.6004L[c.gchld ,]$Speed_lag1
for (i in 1:length(x.spd)) { if (is.na(x.spd[i])==T) x.spd[i]=x.spd[i-1]}
spl.spd<-predict(prior.6004L$sp.spd,x.spd)


#***************************#
#   6004L MODEL FITTING     #
#***************************#
#S/S/S
#source("/home/terrahome/oaj9/_codes/VLaw_Function-with_DLM_pack.r")
MDM.6004L<-list()
occsp<-t5.6004L[c.gchld,]$Occ_lag1*t5.6004L[c.gchld,]$Speed_lag1
for (i in 1:length(occsp)) if (is.na(occsp[i])==T) occsp[i]<-occsp[i-1]
VDF=rep(0.89,3)
MDM.6004L[[1]]<-MDM.VLaw.v12(info=t5.6004L[c.gchld,c("Flow","Date1","Hour")],nparents=1,beta.var=betas.ch.rnd[[2]],Vdf=VDF,
            X_t=t5.6007L[c.chld,]$Flow*spl.prop[1:nrow(spl.prop),],prior.obj=prior.6004L,n.test=n.test,
            parent.f.mean=MDM.6007L[[1]]$f.marg,parent.f.var=MDM.6007L[[1]]$marg.var,
            F_t.star=spl.prop[1:nrow(spl.prop),])

#S/S/F
MDM.6004L[[2]]<-MDM.VLaw.v12(info=t5.6004L[c.gchld,c("Flow","Date1","Hour")],nparents=1,beta.var=betas.ch.rnd[[2]],Vdf=VDF,
            X_t=t5.6007L[c.chld,]$Flow*spl.prop[1:nrow(spl.prop),],prior.obj=prior.6004L,n.test=n.test,
            parent.f.mean=MDM.6007L[[2]]$f.marg,parent.f.var=MDM.6007L[[2]]$marg.var,
            F_t.star=spl.prop[1:nrow(spl.prop),])
#S/F/F
MDM.6004L[[3]]<-MDM.VLaw.v12(info=t5.6004L[c.gchld,c("Flow","Date1","Hour")],nparents=1,beta.var=betas.ch.rnd[[2]],Vdf=VDF,
            X_t=t5.6007L[c.chld,]$Flow*spl.prop[1:nrow(spl.prop),],prior.obj=prior.6004L,n.test=n.test,
            parent.f.mean=MDM.6007L[[4]]$f.marg,parent.f.var=MDM.6007L[[4]]$marg.var,
            F_t.star=spl.prop[1:nrow(spl.prop),])
# F/F/F
MDM.6004L[[4]]<-MDM.VLaw.v12(info=t5.6004L[c.gchld,c("Flow","Date1","Hour")],nparents=1,beta.var=betas.ch.rnd[[2]],Vdf=VDF,
             X_t=as.matrix(cbind(t5.6007L[c.chld,]$Flow*spl.prop[1:nrow(spl.prop),],spl.hway[1:nrow(spl.hway),],
             spl.occ[1:nrow(spl.occ),],spl.spd[1:nrow(spl.spd),],occsp)),prior.obj=prior.6004L,n.test=n.test,
             parent.f.mean=MDM.6007L[[4]]$f.marg,parent.f.var=MDM.6007L[[4]]$marg.var,
             F_t.star=spl.prop[1:nrow(spl.prop),])
#single DLM with F:
sp.flow.6004L<-with(aux.prop,ns(Flow.idx,knots=c(seq(40,100,by=10),seq(180,220,by=10))))
spl.flow<-predict(sp.flow.6004L,t5.6004L[c.gchld,]$Flow.idx)
MDM.6004L[[5]]<-MDM.VLaw.v12(info=t5.6004L[c.gchld,c("Flow","Date1","Hour")],nparents=0,beta.var=betas.ch.rnd[[2]],Vdf=VDF,
   X_t=as.matrix(cbind(spl.flow[1:nrow(spl.flow),],spl.hway[1:nrow(spl.hway),],spl.occ[1:nrow(spl.occ),],
   spl.spd[1:nrow(spl.spd),],occsp)),prior.obj=prior.6007L,n.test=n.test)



r<--(1:n.test)
# Joint LPL 6013B/6007L/6004L:
aux.l<-length(t5.6007L[c.gchld,]$Date[r])
jlpl.6004L<-data.frame(Date=t5.6007L[c.chld,]$Date[r],Hour=t5.6007L[c.chld,]$Hour[r],lpl.DLM=rep(NA,aux.l),
                       lpl.m1=rep(NA,aux.l),lpl.m2=rep(NA,aux.l),lpl.m3=rep(NA,aux.l),lpl.m4=rep(NA,aux.l))
#"S/S/S"
jlpl.6004L$lpl.m1<-dnorm(t5.6013B[c.par,]$Flow[r],mean=dlm.6013B[[1]]$f[r],sd=sqrt(dlm.6013B[[1]]$Q[r]),log=TRUE)+
                   dnorm(t5.6007L[c.chld,]$Flow[r],mean=MDM.6007L[[1]]$f[r],sd=sqrt(MDM.6007L[[1]]$Q[r]),log=TRUE)+
                   dnorm(t5.6004L[c.gchld,]$Flow[r],mean=MDM.6004L[[1]]$f[r],sd=sqrt(MDM.6004L[[1]]$Q[r]),log=TRUE)
# "S/S/F"
jlpl.6004L$lpl.m2<-dnorm(t5.6013B[c.par,]$Flow[r],mean=dlm.6013B[[4]]$f[r],sd=sqrt(dlm.6013B[[4]]$Q[r]),log=TRUE)+
                   dnorm(t5.6007L[c.chld,]$Flow[r],mean=MDM.6007L[[2]]$f[r],sd=sqrt(MDM.6007L[[2]]$Q[r]),log=TRUE)+
                   dnorm(t5.6004L[c.gchld,]$Flow[r],mean=MDM.6004L[[2]]$f[r],sd=sqrt(MDM.6004L[[2]]$Q[r]),log=TRUE)
#"S/F/F"
jlpl.6004L$lpl.m3<-dnorm(t5.6013B[c.par,]$Flow[r],mean=dlm.6013B[[4]]$f[r],sd=sqrt(dlm.6013B[[4]]$Q[r]),log=TRUE)+
                   dnorm(t5.6007L[c.chld,]$Flow[r],mean=MDM.6007L[[4]]$f[r],sd=sqrt(MDM.6007L[[4]]$Q[r]),log=TRUE)+
                   dnorm(t5.6004L[c.gchld,]$Flow[r],mean=MDM.6004L[[3]]$f[r],sd=sqrt(MDM.6004L[[3]]$Q[r]),log=TRUE)
#"F/F/F"
jlpl.6004L$lpl.m4<-dnorm(t5.6013B[c.par,]$Flow[r],mean=dlm.6013B[[4]]$f[r],sd=sqrt(dlm.6013B[[4]]$Q[r]),log=TRUE)+
                   dnorm(t5.6007L[c.chld,]$Flow[r],mean=MDM.6007L[[4]]$f[r],sd=sqrt(MDM.6007L[[4]]$Q[r]),log=TRUE)+
                   dnorm(t5.6004L[c.gchld,]$Flow[r],mean=MDM.6004L[[4]]$f[r],sd=sqrt(MDM.6004L[[4]]$Q[r]),log=TRUE)

#"F/F/F DLM"
jlpl.6004L$lpl.DLM<-dnorm(t5.6013B[c.par,]$Flow[r],mean=dlm.6013B[[4]]$f[r],sd=sqrt(dlm.6013B[[4]]$Q[r]),log=TRUE)+
                   dnorm(t5.6007L[c.chld,]$Flow[r],mean=MDM.6007L[[5]]$f[r],sd=sqrt(MDM.6007L[[5]]$Q[r]),log=TRUE)+
                   dnorm(t5.6004L[c.gchld,]$Flow[r],mean=MDM.6004L[[5]]$f[r],sd=sqrt(MDM.6004L[[5]]$Q[r]),log=TRUE)

(jlpl.6004L<-data.frame(model=c("S/S/S","S/S/F","S/F/F","F/F/F","F/F/F DLM"),
                           lpl=c(sum(jlpl.6004L$lpl.m1,na.rm="T"),sum(jlpl.6004L$lpl.m2,na.rm="T"),
                                 sum(jlpl.6004L$lpl.m3,na.rm="T"),sum(jlpl.6004L$lpl.m4,na.rm="T"),
                                 sum(jlpl.6004L$lpl.DLM,na.rm="T"))))



save(MDM.6004L,file="/home/terrahome/oaj9/splines/MDM.6004L.Rdata")
save(prior.6004L,file="/home/terrahome/oaj9/splines/prior.6004L.Rdata")







