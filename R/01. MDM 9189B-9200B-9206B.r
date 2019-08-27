#+++++++++++++++++++++++++++++++++++++++++++++++++++++#
#     Spline MDM                                      #
#     DLM for 9206B                                   #
#     Set/2011                                        #
#+++++++++++++++++++++++++++++++++++++++++++++++++++++#
library(chron);library(splines);library(car)
source("/home/terrahome/oaj9/_codes/MDM-VLaw-v12.r")
load("/home/terrahome/oaj9/_data/all2010/t5.9206B.Rdata")
load("/home/terrahome/oaj9/_data/all2010/t5.9200B.Rdata")
load("/home/terrahome/oaj9/_data/all2010/t5.9189B.Rdata")
load(file="/home/terrahome/oaj9/paper 1 review/betas.rnd.Rdata")
load(file="/home/terrahome/oaj9/paper 1 review/betas.ch.rnd.Rdata")
load("/home/terrahome/oaj9/_data/all2010/stats/betas.9206B.Rdata")
time.filter<-6:20
test.dates<-c("18/08/2010","25/08/2010")
prior.months<- c("Jul","Aug")
#learning.months<-c("Sep","Oct","Nov")
learning.months<-c("Sep","Oct","Nov","Dec")
wkday<-"Wed"
#**************************************************************************************************#
#-------------------------------------9206B MODELLING----------------------------------------------#
#**************************************************************************************************#

#-------------------------#
#     Data Selection      #
#-------------------------#
# prior and learning data:
t5.9206B$Flow.idx<-rep(1:288,dim(table(t5.9206B$Date)))
time.filter<-6:20
c.pri.par<-months(t5.9206B$Date1) %in% prior.months & !t5.9206B$Date %in% test.dates & t5.9206B$Hour %in% time.filter
c.par<-((months(t5.9206B$Date1) %in% learning.months | t5.9206B$Date %in% test.dates) &  !t5.9206B$Date %in% c("24/11/2010"))  & t5.9206B$Hour %in% time.filter  & t5.9206B$Dweek==wkday
period<-dim(table(t5.9206B[c.pri.par,]$Time5))
t5.9206B[is.na(t5.9206B$Hway_lag1)==T,]$Flow<-NA
t5.9206B[is.na(t5.9206B$Occ_lag1)==T,]$Flow<-NA
n.test<-period*length(test.dates)

with(t5.9206B[c.pri.par,],table(months(Date1)))
with(t5.9206B[c.par,],as.matrix(table(Date1)))





#----------#
#  priors  #
#----------#
prior.9206B<-list()
prior.9206B$sigma2<-with(t5.9206B[c.pri.par,], prior.5min(Flow,t5.9206B[c.pri.par,],dim(table(Date1)),period)$S0)
prior.9206B$mean.prior<-with(t5.9206B[c.pri.par & t5.9206B$Time5==t5.9206B[c.pri.par,]$Time5[1],],
                      mean(Flow,na.rm=T))
prior.9206B$total.flow<-sum(t5.9206B[t5.9206B$Date =="04/10/2010",]$Flow,na.rm=T)
#--------------------------------------#
#          Spline fitting              #
#--------------------------------------#
c.spl.par<-c.pri.par  & is.na(t5.9206B$Flow)==F & is.na(t5.9206B$Occ_lag1)==F  & is.na(t5.9206B$Hway_lag1)==F  & is.na(t5.9206B$Speed_lag1)==F
c.spl.flow<-c.pri.par &  t5.9206B$Dweek==wkday & is.na(t5.9206B$Flow)==F
with(t5.9206B[c.spl.par,],table(months(Date1)))
with(t5.9206B[c.spl.flow,],table(months(Date1)))
prior.9206B$sp.flow<-with(t5.9206B[c.spl.flow,],ns(Flow.idx,knots=c(seq(40,130,by=5),seq(180,220,by=10))))
prior.9206B$sp.occ<-with(t5.9206B[c.spl.par,],ns(Occ_lag1,df=4))
prior.9206B$sp.hway<-with(t5.9206B[c.spl.par & t5.9206B$Hway_lag1>1,],ns(Hway_lag1,df=4))
prior.9206B$sp.spd<-with(t5.9206B[c.spl.par,],ns(Speed_lag1,df=4))


#--------------------------------------#
#          Spline prediction           #
#--------------------------------------#
# removing missing (otherwise predict function won't work - obs won't be used anyway)

x.flow<-t5.9206B[c.par,]$Flow.idx
for (i in 1:length(x.flow)) { if (is.na(x.flow[i])==T) x.flow[i]=x.flow[i-1]}
spl.flow<-predict(prior.9206B$sp.flow,x.flow)

x.occ<-t5.9206B[c.par,]$Occ_lag1
for (i in 1:length(x.occ)) { if (is.na(x.occ[i])==T) x.occ[i]=x.occ[i-1]}
spl.occ<-predict(prior.9206B$sp.occ,x.occ)

x.hway<-t5.9206B[c.par & t5.9206B$Hway_lag1>1,]$Hway_lag1
for (i in 1:length(x.hway)) { if (is.na(x.hway[i])==T) x.hway[i]=x.hway[i-1]}
spl.hway<-predict(prior.9206B$sp.hway,x.hway)

x.spd<-t5.9206B[c.par ,]$Speed_lag1
for (i in 1:length(x.spd)) { if (is.na(x.spd[i])==T) x.spd[i]=x.spd[i-1]}
spl.spd<-predict(prior.9206B$sp.spd,x.spd)

#***************************#
#   9206B MODEL FITTING     #
#***************************#
dlm.9206B<-list()
VDF=rep(0.99,3)



#   S:  Flow spline trend
source("/home/terrahome/oaj9/_codes/MDM-VLaw-v12.r")
ini<-proc.time()
dlm.9206B[[1]]<-MDM.VLaw.v12(info=t5.9206B[c.par,c("Flow","Date1","Hour")],nparents=0,beta.var=betas.rnd[[1]],Vdf=VDF,
                              X_t=as.matrix(spl.flow[1:nrow(spl.flow),]),prior.obj=prior.9206B,n.test=n.test)
#---------individual models----------------------------------#
#   SO  (S + occupancy)
source("/home/terrahome/oaj9/_codes/MDM-VLaw-v12.r")
dlm.9206B[[2]]<-MDM.VLaw.v12(info=t5.9206B[c.par,c("Flow","Date1","Hour")],nparents=0,beta.var=betas.rnd[[2]],Vdf=VDF,
                             X_t=as.matrix(cbind(spl.flow[1:nrow(spl.flow),],spl.occ[1:nrow(spl.occ),])),ind.vlaw=T,
                             prior.obj=prior.9206B,n.test=n.test)
#   SH (S + headway)
dlm.9206B[[3]]<-MDM.VLaw.v12(info=t5.9206B[c.par,c("Flow","Date1","Hour")],nparents=0,beta.var=betas.rnd[[1]],Vdf=VDF,
                             X_t=as.matrix(cbind(spl.flow[1:nrow(spl.flow),],spl.hway[1:nrow(spl.hway),])),
                             prior.obj=prior.9206B,n.test=n.test)
#   SHinv  (S + 1/headway)
dlm.9206B[[4]]<-MDM.VLaw.v12(info=t5.9206B[c.par,c("Flow","Date1","Hour")],nparents=0,beta.var=betas.rnd[[1]],Vdf=VDF,
                             X_t=as.matrix(cbind(spl.flow[1:nrow(spl.flow),],1/t5.9206B[c.par,]$Hway_lag1)),
                             prior.obj=prior.9206B,n.test=n.test)
# - SSp (S + Speed)
dlm.9206B[[5]]<-MDM.VLaw.v12(info=t5.9206B[c.par,c("Flow","Date1","Hour")],nparents=0,beta.var=betas.rnd[[1]],Vdf=VDF,
                             X_t=as.matrix(cbind(spl.flow[1:nrow(spl.flow),],spl.spd[1:nrow(spl.spd),])),
                             prior.obj=prior.9206B,n.test=n.test)
#------Interaction test---------------------------------------#
# - SOSp: (S + occupancy + speed)
dlm.9206B[[6]]<-MDM.VLaw.v12(info=t5.9206B[c.par,c("Flow","Date1","Hour")],nparents=0,beta.var=betas.rnd[[1]],Vdf=VDF,
                             X_t=as.matrix(cbind(spl.flow[1:nrow(spl.flow),],spl.occ[1:nrow(spl.occ),],spl.spd[1:nrow(spl.spd),])),
                             prior.obj=prior.9206B,n.test=n.test)
# - S(OSp) : (S + occupancy + speed + occ x speed)
dlm.9206B[[7]]<-MDM.VLaw.v12(info=t5.9206B[c.par,c("Flow","Date1","Hour")],nparents=0,beta.var=betas.rnd[[1]],Vdf=VDF,
                             X_t=as.matrix(cbind(spl.flow[1:nrow(spl.flow),],spl.occ[1:nrow(spl.occ),],spl.spd[1:nrow(spl.spd),],
                                                  t5.9206B[c.par,]$Occ_lag1*t5.9206B[c.par,]$Speed_lag1)),
                             prior.obj=prior.9206B,n.test=n.test)
#------models with all variables:------------------------------#
# - Spline Flow Level + O + H + Sp
dlm.9206B[[8]]<-MDM.VLaw.v12(info=t5.9206B[c.par,c("Flow","Date1","Hour")],nparents=0,beta.var=betas.rnd[[1]],Vdf=VDF,
                             X_t=as.matrix(cbind(spl.flow[1:nrow(spl.flow),],spl.occ[1:nrow(spl.occ),],
                                                 spl.hway[1:nrow(spl.hway),],spl.spd[1:nrow(spl.spd),])),
                             prior.obj=prior.9206B,n.test=n.test)
# - Spline Flow Level +  Sp + Occ + (OSp) +  Hway
dlm.9206B[[9]]<-MDM.VLaw.v12(info=t5.9206B[c.par,c("Flow","Date1","Hour")],nparents=0,beta.var=betas.rnd[[1]],Vdf=VDF,
                             X_t=as.matrix(cbind(spl.flow[1:nrow(spl.flow),],spl.hway[1:nrow(spl.hway),],spl.occ[1:nrow(spl.occ),],
                                                  spl.spd[1:nrow(spl.spd),],t5.9206B[c.par,]$Occ_lag1*t5.9206B[c.par,]$Speed_lag1)),
                             prior.obj=prior.9206B,n.test=n.test)


(res2.9206B<-data.frame(Model=c("S","SO","SH","SHinv","SSp","SOSp","S(OSp)","SOSpH","S(OSp)H"),
            LPL=c(dlm.9206B[[1]]$perf[5,2],dlm.9206B[[2]]$perf[5,2], dlm.9206B[[3]]$perf[5,2],
                  dlm.9206B[[4]]$perf[5,2],dlm.9206B[[5]]$perf[5,2], dlm.9206B[[6]]$perf[5,2],
                  dlm.9206B[[7]]$perf[5,2],dlm.9206B[[8]]$perf[5,2], dlm.9206B[[9]]$perf[5,2])))


save(dlm.9206B,file="/home/terrahome/oaj9/splines/dlm.9206B.Rdata")
save(prior.9206B,file="/home/terrahome/oaj9/splines/prior.9206B.Rdata")

#**************************************************************************************************#
#-------------------------------------9200B MODELLING----------------------------------------------#
#**************************************************************************************************#
#load("/home/terrahome/oaj9/splines/dlm.9206B.Rdata")
#-------------------------#
#     Data Selection      #
#-------------------------#
t5.9200B$Flow.idx<-rep(1:288,dim(table(t5.9200B$Date)))
# prior and learning data:
period<-dim(table(t5.9200B[c.pri.chld,]$Time5))
c.pri.chld<-months(t5.9200B$Date1) %in% prior.months & !t5.9200B$Date %in% test.dates & t5.9200B$Hour %in% time.filter
c.chld<-(months(t5.9200B$Date1) %in% learning.months | t5.9200B$Date %in% test.dates) & t5.9200B$Hour %in% time.filter  & t5.9200B$Dweek==wkday
n.test<-period*length(test.dates)

with(t5.9200B[c.pri.chld,],table(months(Date1)))
with(t5.9200B[c.chld,],as.matrix(table(Date1)))

#----------#
#  priors  #
#----------#
prior.9200B<-MDM.DLM.find.prior(t5.9200B[c.pri.chld,]$Flow,t5.9206B[c.pri.par,]$Flow,period,0.001)
prior.9200B$sigma2<-sum((t5.9200B[c.pri.chld,]$Flow-t5.9206B[c.pri.par,]$Flow)^2,na.rm=T)/
                         length(t5.9200B[c.pri.chld,]$Flow)
prior.9200B$mean.prior<-with(t5.9200B[c.pri.chld & t5.9200B$Time5==t5.9200B[c.pri.chld,]$Time5[1],],
                      mean(Flow,na.rm=T))

# total flow per period (denominator of the vlaw):
prior.9200B$total.flow<-sum(t5.9200B[t5.9200B$Date =="04/10/2010",]$Flow,na.rm=T)

#--------------------------------------#
#          Spline fitting              #
#--------------------------------------#
c.chld.spl<-c.pri.chld & is.na(t5.9200B$Flow_lag1)==F & is.na(t5.9200B$Flow)==F & is.na(t5.9200B$Occ_lag1)==F &  is.na(t5.9200B$Hway_lag1)==F &  is.na(t5.9200B$Speed_lag1)==F
prior.9200B$sp.occ<-with(t5.9200B[c.chld.spl ,],ns(Occ_lag1,df=3))
prior.9200B$sp.hway<-with(t5.9200B[c.chld.spl,],ns(Hway_lag1,df=3))
prior.9200B$sp.spd<-with(t5.9200B[c.chld.spl,],ns(Speed_lag1,df=3))

aux.prop<-data.frame(Flow.idx=t5.9200B[c.pri.chld,]$Flow.idx,Dweek=t5.9200B[c.pri.chld,]$Dweek,
                     Flow.par=t5.9206B[c.pri.par,]$Flow,Flow.chd=t5.9200B[c.pri.chld,]$Flow,
                     Flow.prop=t5.9200B[c.pri.chld,]$Flow/t5.9206B[c.pri.par,]$Flow)
aux.prop<-aux.prop[is.na(aux.prop$Flow.prop)==F & aux.prop$Dweek==wkday ,]
prior.9200B$sp.prop<-with(aux.prop,ns(Flow.idx,knots=c(seq(60,120,by=15),seq(180,220,by=15))))

#--------------------------------------#
#          Spline prediction           #
#--------------------------------------#
spl.prop<-predict(prior.9200B$sp.prop,t5.9200B[c.chld,]$Flow.idx)
# removing missing (otherwise predict function won't work - obs won't be used anyway)
x.occ<-t5.9200B[c.chld,]$Occ_lag1
for (i in 1:length(x.occ)) { if (is.na(x.occ[i])==T) x.occ[i]=x.occ[i-1]}
spl.occ<-predict(prior.9200B$sp.occ,x.occ)
x.hway<-t5.9200B[c.chld ,]$Hway_lag1
for (i in 1:length(x.hway)) { if (is.na(x.hway[i])==T) x.hway[i]=x.hway[i-1]}
spl.hway<-predict(prior.9200B$sp.hway,x.hway)
x.spd<-t5.9200B[c.par ,]$Speed_lag1
for (i in 1:length(x.spd)) { if (is.na(x.spd[i])==T) x.spd[i]=x.spd[i-1]}
spl.spd<-predict(prior.9200B$sp.spd,x.spd)


#***************************#
#   9200B MODEL FITTING     #
#***************************#
VDF=rep(0.99,3)
MDM.9200B<-list()
#  Spline proportion (S) and parent with spline flow (S) : S/S Model
MDM.9200B[[1]]<-MDM.VLaw.v12(info=t5.9200B[c.chld,c("Flow","Date1","Hour")],nparents=1,beta.var=betas.ch.rnd[[1]],Vdf=VDF,
             X_t=t5.9206B[c.par,]$Flow*spl.prop[1:nrow(spl.prop),],F_t.star=spl.prop[1:nrow(spl.prop),],ind.vlaw=F,
             parent.f.mean=dlm.9206B[[1]]$f,parent.f.var=dlm.9206B[[1]]$Q,prior.obj=prior.9200B,n.test=n.test)
table(MDM.9200B[[1]]$f<0)
             
             
             
             
             
#  Spline proportion (S) and parent with Full Model    S/F
MDM.9200B[[2]]<-MDM.VLaw.v12(info=t5.9200B[c.chld,c("Flow","Date1","Hour")],nparents=1,beta.var=betas.ch.rnd[[1]],Vdf=VDF,
             X_t=t5.9206B[c.par,]$Flow*spl.prop[1:nrow(spl.prop),],F_t.star=spl.prop[1:nrow(spl.prop),],
             parent.f.mean=dlm.9206B[[9]]$f,parent.f.var=dlm.9206B[[9]]$Q,prior.obj=prior.9200B,n.test=n.test)
# Full model and parent with spline flow, : F/S Model
MDM.9200B[[3]]<-MDM.VLaw.v12(info=t5.9200B[c.chld,c("Flow","Date1","Hour")],nparents=1,beta.var=betas.ch.rnd[[1]],Vdf=VDF,
             X_t=as.matrix(cbind(t5.9206B[c.par,]$Flow*spl.prop[1:nrow(spl.prop),],spl.hway[1:nrow(spl.hway),],
             spl.occ[1:nrow(spl.occ),],spl.spd[1:nrow(spl.spd),])),#t5.9200B[c.chld,]$Occ_lag1*t5.9200B[c.chld,]$Speed_lag1)),
             F_t.star=spl.prop[1:nrow(spl.prop),],parent.f.mean=dlm.9206B[[1]]$f,parent.f.var=dlm.9206B[[1]]$Q,prior.obj=prior.9200B,n.test=n.test)
# Full model and parent with Full model F/F Model
MDM.9200B[[4]]<-MDM.VLaw.v12(info=t5.9200B[c.chld,c("Flow","Date1","Hour")],nparents=1,beta.var=betas.ch.rnd[[1]],Vdf=VDF,
               X_t=as.matrix(cbind(t5.9206B[c.par,]$Flow*spl.prop[1:nrow(spl.prop),],spl.hway[1:nrow(spl.hway),],
               spl.occ[1:nrow(spl.occ),],spl.spd[1:nrow(spl.spd),])),#t5.9200B[c.chld,]$Occ_lag1*t5.9200B[c.chld,]$Speed_lag1)),
               F_t.star=spl.prop[1:nrow(spl.prop),],parent.f.mean=dlm.9206B[[9]]$f,parent.f.var=dlm.9206B[[9]]$Q,prior.obj=prior.9200B,n.test=n.test)
# DLM F/F
sp.flow.9200B<-with(aux.prop,ns(Flow.idx,knots=c(seq(40,100,by=10),seq(180,220,by=10))))
spl.flow<-predict(sp.flow.9200B,t5.9200B[c.chld,]$Flow.idx)
MDM.9200B[[5]]<-MDM.VLaw.v12(info=t5.9200B[c.chld,c("Flow","Date1","Hour")],nparents=0,beta.var=betas.ch.rnd[[1]],Vdf=VDF,
     X_t=as.matrix(cbind(spl.flow[1:nrow(spl.flow),],spl.hway[1:nrow(spl.hway),],spl.occ[1:nrow(spl.occ),],
     spl.spd[1:nrow(spl.spd),])),prior.obj=prior.9200B,n.test=n.test)
#,t5.9200B[c.par,]$Occ_lag1*t5.9200B[c.par,]$Speed_lag1

r<--(1:n.test)
#-----------------------#
# Joint LPL 9206B/9200B:#
#-----------------------#
aux.l<-length(t5.9200B[c.chld,]$Date[r])
jlpl.9200B<-data.frame(Date=t5.9200B[c.chld,]$Date[r],Hour=t5.9200B[c.chld,]$Hour[r],lpl.DLM=rep(NA,aux.l),
                       lpl.m1=rep(NA,aux.l),lpl.m2=rep(NA,aux.l),lpl.m3=rep(NA,aux.l),lpl.m4=rep(NA,aux.l))

#"S/S"
jlpl.9200B$lpl.m1<-dnorm(t5.9206B[c.par,]$Flow[r],mean=dlm.9206B[[1]]$f[r],sd=sqrt(dlm.9206B[[1]]$Q[r]),log=TRUE)+
                   dnorm(t5.9200B[c.chld,]$Flow[r],mean=MDM.9200B[[1]]$f[r],sd=sqrt(MDM.9200B[[1]]$Q[r]),log=TRUE)
#"S/F
jlpl.9200B$lpl.m2<-dnorm(t5.9206B[c.par,]$Flow[r],mean=dlm.9206B[[9]]$f[r],sd=sqrt(dlm.9206B[[9]]$Q[r]),log=TRUE)+
                   dnorm(t5.9200B[c.chld,]$Flow[r],mean=MDM.9200B[[2]]$f[r],sd=sqrt(MDM.9200B[[2]]$Q[r]),log=TRUE)
#"F/S"
jlpl.9200B$lpl.m3<-dnorm(t5.9206B[c.par,]$Flow[r],mean=dlm.9206B[[1]]$f[r],sd=sqrt(dlm.9206B[[1]]$Q[r]),log=TRUE)+
                   dnorm(t5.9200B[c.chld,]$Flow[r],mean=MDM.9200B[[3]]$f[r],sd=sqrt(MDM.9200B[[3]]$Q[r]),log=TRUE)
#"F/F"
jlpl.9200B$lpl.m4<-dnorm(t5.9206B[c.par,]$Flow[r],mean=dlm.9206B[[9]]$f[r],sd=sqrt(dlm.9206B[[9]]$Q[r]),log=TRUE)+
                   dnorm(t5.9200B[c.chld,]$Flow[r],mean=MDM.9200B[[4]]$f[r],sd=sqrt(MDM.9200B[[4]]$Q[r]),log=TRUE)

#"F/F DLM"
jlpl.9200B$lpl.DLM<-dnorm(t5.9206B[c.par,]$Flow[r],mean=dlm.9206B[[9]]$f[r],sd=sqrt(dlm.9206B[[9]]$Q[r]),log=TRUE)+
                   dnorm(t5.9200B[c.chld,]$Flow[r],mean=MDM.9200B[[5]]$f[r],sd=sqrt(MDM.9200B[[5]]$Q[r]),log=TRUE)

(jlpl.9200B<-data.frame(model=c("S/S","S/F","F/S","F/F","F/F DLM"),
                           lpl=c(sum(jlpl.9200B$lpl.m1,na.rm="T"),sum(jlpl.9200B$lpl.m2,na.rm="T"),
                                 sum(jlpl.9200B$lpl.m3,na.rm="T"),sum(jlpl.9200B$lpl.m4,na.rm="T"),
                                 sum(jlpl.9200B$lpl.DLM,na.rm="T"))))
                                 
save(MDM.9200B,file="/home/terrahome/oaj9/splines/MDM.9200B.Rdata")
save(prior.9200B,file="/home/terrahome/oaj9/splines/prior.9200B.Rdata")

#**************************************************************************************************#
#-------------------------------------9189B MODELLING----------------------------------------------#
#**************************************************************************************************#
#load("/home/terrahome/oaj9/splines/MDM.9200B.Rdata")
#-------------------------#
#     Data Selection      #
#-------------------------#
t5.9189B$Flow.idx<-rep(1:288,dim(table(t5.9189B$Date)))
# prior and learning data:
c.pri.gchld<-months(t5.9189B$Date1) %in% prior.months & !t5.9189B$Date %in% test.dates & t5.9189B$Hour %in% time.filter
c.gchld<-(months(t5.9189B$Date1) %in% learning.months | t5.9189B$Date %in% test.dates) & t5.9189B$Hour %in% time.filter  & t5.9189B$Dweek==wkday
period<-dim(table(t5.9189B[c.pri.gchld,]$Time5))
n.test<-period*length(test.dates)
with(t5.9189B[c.pri.gchld,],table(months(Date1)))
with(t5.9189B[c.gchld,],as.matrix(table(Date1)))
#----------#
#  priors  #
#----------#
#prior.9189B<-MDM.DLM.find.prior(t5.9189B[c.pri.gchld,]$Flow,
#                                     t5.9200B[c.pri.chld,]$Flow, period,0.001)
prior.9189B<-list()
prior.9189B$sigma2<-sum((t5.9189B[c.pri.gchld,]$Flow-t5.9200B[c.pri.chld,]$Flow)^2,na.rm=T)/
                  length(t5.9189B[c.pri.gchld,]$Flow)
prior.9189B$mean.prior<-with(t5.9189B[c.pri.gchld & t5.9189B$Time5==t5.9189B[c.pri.gchld,]$Time5[1],],
                      mean(Flow,na.rm=T))
# total flow per period (denominator of the vlaw):
prior.9189B$total.flow<-sum(t5.9189B[t5.9189B$Date =="04/10/2010",]$Flow,na.rm=T)

#--------------------------------------#
#          Spline fitting              #
#--------------------------------------#
c.gchld.spl<-c.pri.gchld & is.na(t5.9189B$Flow_lag1)==F  & is.na(t5.9189B$Flow)==F & is.na(t5.9189B$Occ_lag1)==F &  is.na(t5.9189B$Hway_lag1)==F &  is.na(t5.9189B$Speed_lag1)==F
prior.9189B$sp.occ<-with(t5.9189B[c.gchld.spl ,],ns(Occ_lag1,df=4))
prior.9189B$sp.hway<-with(t5.9189B[c.gchld.spl & t5.9189B$Hway_lag1>1.3 ,],ns(Hway_lag1,df=4))
prior.9189B$sp.spd<-with(t5.9189B[c.gchld.spl,],ns(Speed_lag1,df=4))

aux.prop<-data.frame(Flow.idx=t5.9189B[c.pri.gchld,]$Flow.idx,Dweek=t5.9189B[c.pri.gchld,]$Dweek,
                     Flow.par=t5.9200B[c.pri.chld,]$Flow,Flow.chd=t5.9189B[c.pri.chld,]$Flow,
                     Flow.prop=t5.9189B[c.pri.gchld,]$Flow/t5.9200B[c.pri.chld,]$Flow)
aux.prop<-aux.prop[is.na(aux.prop$Flow.prop)==F & aux.prop$Dweek==wkday ,]
prior.9189B$sp.prop<-with(aux.prop,ns(Flow.idx,knots=c(seq(60,120,by=10),seq(180,220,by=10))))
#--------------------------------------#
#          Spline prediction           #
#--------------------------------------#
spl.prop<-predict(prior.9189B$sp.prop,t5.9189B[c.gchld,]$Flow.idx)
# removing missing (otherwise predict function won't work - obs won't be used anyway)
x.occ<-t5.9189B[c.gchld,]$Occ_lag1
for (i in 1:length(x.occ)) { if (is.na(x.occ[i])==T) x.occ[i]=x.occ[i-1]}
spl.occ<-predict(prior.9189B$sp.occ,x.occ)
x.hway<-t5.9189B[c.gchld,]$Hway_lag1
for (i in 1:length(x.hway)) { if (is.na(x.hway[i])==T) x.hway[i]=x.hway[i-1]}
spl.hway<-predict(prior.9189B$sp.hway,x.hway)
x.spd<-t5.9189B[c.gchld,]$Speed_lag1
for (i in 1:length(x.spd)) { if (is.na(x.spd[i])==T) x.spd[i]=x.spd[i-1]}
spl.spd<-predict(prior.9189B$sp.spd,x.spd)



#***************************#
#   9189B MODEL FITTING     #
#***************************#
#S/S/S
VDF=rep(0.88,3)
MDM.9189B<-list()
#source("/home/terrahome/oaj9/_codes/VLaw_Function-with_DLM_pack.r")
#S/S/S
MDM.9189B<-list()
MDM.9189B[[1]]<-MDM.VLaw.v12(info=t5.9189B[c.gchld,c("Flow","Date1","Hour")],nparents=1,beta.var=betas.ch.rnd[[1]],
            X_t=t5.9200B[c.chld,]$Flow*spl.prop[1:nrow(spl.prop),],prior.obj=prior.9189B,n.test=n.test,
            parent.f.mean=MDM.9200B[[1]]$f.marg,parent.f.var=MDM.9200B[[1]]$marg.var,Vdf=VDF,
            F_t.star=spl.prop[1:nrow(spl.prop),])

#S/S/F
MDM.9189B[[2]]<-MDM.VLaw.v12(info=t5.9189B[c.gchld,c("Flow","Date1","Hour")],nparents=1,beta.var=betas.ch.rnd[[1]],
            X_t=t5.9200B[c.chld,]$Flow*spl.prop[1:nrow(spl.prop),],prior.obj=prior.9189B,n.test=n.test,
            parent.f.mean=MDM.9200B[[2]]$f.marg,parent.f.var=MDM.9200B[[2]]$marg.var,Vdf=VDF,
            F_t.star=spl.prop[1:nrow(spl.prop),])
#S/F/F
MDM.9189B[[3]]<-MDM.VLaw.v12(info=t5.9189B[c.gchld,c("Flow","Date1","Hour")],nparents=1,beta.var=betas.ch.rnd[[1]],
            X_t=t5.9200B[c.chld,]$Flow*spl.prop[1:nrow(spl.prop),],prior.obj=prior.9189B,n.test=n.test,
            parent.f.mean=MDM.9200B[[4]]$f.marg,parent.f.var=MDM.9200B[[4]]$marg.var,Vdf=VDF,
            F_t.star=spl.prop[1:nrow(spl.prop),])
# F/F/F
MDM.9189B[[4]]<-MDM.VLaw.v12(info=t5.9189B[c.gchld,c("Flow","Date1","Hour")],nparents=1,beta.var=betas.ch.rnd[[1]],
             X_t=as.matrix(cbind(t5.9200B[c.chld,]$Flow*spl.prop[1:nrow(spl.prop),],spl.hway[1:nrow(spl.hway),],
             spl.occ[1:nrow(spl.occ),],spl.spd[1:nrow(spl.spd),])), #t5.9189B[c.gchld,]$Occ_lag1*t5.9189B[c.gchld,]$Speed_lag1)),
             prior.obj=prior.9189B,n.test=n.test,F_t.star=spl.prop[1:nrow(spl.prop),],
             parent.f.mean=MDM.9200B[[4]]$f.marg,parent.f.var=MDM.9200B[[4]]$marg.var,Vdf=VDF)

#single DLM with SOH:
sp.flow.9189B<-with(aux.prop,ns(Flow.idx,knots=c(seq(40,100,by=10),seq(180,220,by=10))))
spl.flow<-predict(sp.flow.9189B,t5.9189B[c.gchld,]$Flow.idx)
MDM.9189B[[5]]<-MDM.VLaw.v12(info=t5.9189B[c.gchld,c("Flow","Date1","Hour")],nparents=0,beta.var=betas.ch.rnd[[1]],Vdf=VDF,
   X_t=as.matrix(cbind(spl.flow[1:nrow(spl.flow),],spl.occ[1:nrow(spl.occ),],spl.hway[1:nrow(spl.hway),],
   spl.occ[1:nrow(spl.occ),],spl.spd[1:nrow(spl.spd),])),#t5.9189B[c.gchld,]$Occ_lag1*t5.9189B[c.gchld,]$Speed_lag1)),
    prior.obj=prior.9200B,n.test=n.test)
                           

r<--(1:n.test)
# Joint LPL 9206B/9200B/9189B:
aux.l<-length(t5.9200B[c.gchld,]$Date[r])
jlpl.9189B<-data.frame(Date=t5.9200B[c.chld,]$Date[r],Hour=t5.9200B[c.chld,]$Hour[r],lpl.DLM=rep(NA,aux.l),
                       lpl.m1=rep(NA,aux.l),lpl.m2=rep(NA,aux.l),lpl.m3=rep(NA,aux.l),lpl.m4=rep(NA,aux.l))
#"S/S/S"
jlpl.9189B$lpl.m1<-dnorm(t5.9206B[c.par,]$Flow[r],mean=dlm.9206B[[1]]$f[r],sd=sqrt(dlm.9206B[[1]]$Q[r]),log=TRUE)+
                   dnorm(t5.9200B[c.chld,]$Flow[r],mean=MDM.9200B[[1]]$f[r],sd=sqrt(MDM.9200B[[1]]$Q[r]),log=TRUE)+
                   dnorm(t5.9189B[c.gchld,]$Flow[r],mean=MDM.9189B[[1]]$f[r],sd=sqrt(MDM.9189B[[1]]$Q[r]),log=TRUE)
# "S/S/F"
jlpl.9189B$lpl.m2<-dnorm(t5.9206B[c.par,]$Flow[r],mean=dlm.9206B[[9]]$f[r],sd=sqrt(dlm.9206B[[9]]$Q[r]),log=TRUE)+
                   dnorm(t5.9200B[c.chld,]$Flow[r],mean=MDM.9200B[[2]]$f[r],sd=sqrt(MDM.9200B[[2]]$Q[r]),log=TRUE)+
                   dnorm(t5.9189B[c.gchld,]$Flow[r],mean=MDM.9189B[[2]]$f[r],sd=sqrt(MDM.9189B[[2]]$Q[r]),log=TRUE)
#"S/F/F"
jlpl.9189B$lpl.m3<-dnorm(t5.9206B[c.par,]$Flow[r],mean=dlm.9206B[[9]]$f[r],sd=sqrt(dlm.9206B[[9]]$Q[r]),log=TRUE)+
                   dnorm(t5.9200B[c.chld,]$Flow[r],mean=MDM.9200B[[4]]$f[r],sd=sqrt(MDM.9200B[[4]]$Q[r]),log=TRUE)+
                   dnorm(t5.9189B[c.gchld,]$Flow[r],mean=MDM.9189B[[3]]$f[r],sd=sqrt(MDM.9189B[[3]]$Q[r]),log=TRUE)
#"F/F/F"
jlpl.9189B$lpl.m4<-dnorm(t5.9206B[c.par,]$Flow[r],mean=dlm.9206B[[9]]$f[r],sd=sqrt(dlm.9206B[[9]]$Q[r]),log=TRUE)+
                   dnorm(t5.9200B[c.chld,]$Flow[r],mean=MDM.9200B[[4]]$f[r],sd=sqrt(MDM.9200B[[4]]$Q[r]),log=TRUE)+
                   dnorm(t5.9189B[c.gchld,]$Flow[r],mean=MDM.9189B[[4]]$f[r],sd=sqrt(MDM.9189B[[4]]$Q[r]),log=TRUE)

#"F/F/F DLM"
jlpl.9189B$lpl.DLM<-dnorm(t5.9206B[c.par,]$Flow[r],mean=dlm.9206B[[9]]$f[r],sd=sqrt(dlm.9206B[[9]]$Q[r]),log=TRUE)+
                   dnorm(t5.9200B[c.chld,]$Flow[r],mean=MDM.9200B[[5]]$f[r],sd=sqrt(MDM.9200B[[5]]$Q[r]),log=TRUE)+
                   dnorm(t5.9189B[c.gchld,]$Flow[r],mean=MDM.9189B[[5]]$f[r],sd=sqrt(MDM.9189B[[5]]$Q[r]),log=TRUE)

(jlpl.9189B<-data.frame(model=c("S/S/S","S/S/F","S/F/F","F/F/F","F/F/F DLM"),
                           lpl=c(sum(jlpl.9189B$lpl.m1,na.rm="T"),sum(jlpl.9189B$lpl.m2,na.rm="T"),
                                 sum(jlpl.9189B$lpl.m3,na.rm="T"),sum(jlpl.9189B$lpl.m4,na.rm="T"),
                                 sum(jlpl.9189B$lpl.DLM,na.rm="T"))))



save(MDM.9189B,file="/home/terrahome/oaj9/splines/MDM.9189B.Rdata")
save(prior.9189B,file="/home/terrahome/oaj9/splines/prior.9189B.Rdata")





 table(is.na(MDM.9189B[[2]]$Q))