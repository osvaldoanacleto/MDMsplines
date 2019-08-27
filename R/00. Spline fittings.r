load("C:\\data\\all2010\\t5.9206B.Rdata")
load("C:\\data\\all2010\\t5.6013B.Rdata")
load("C:\\data\\all2010\\t5.9188A.Rdata")
load("C:\\data\\all2010\\t5.1431A.Rdata")
load("C:\\data\\all2010\\t5.9200B.Rdata")
load("C:\\data\\all2010\\t5.6007L.Rdata")
load("C:\\data\\all2010\\t5.9193J.Rdata")
load("C:\\data\\all2010\\t5.1437A.Rdata")
load("C:\\data\\all2010\\t5.9189B.Rdata")
load("C:\\data\\all2010\\t5.6004L.Rdata")
load("C:\\data\\all2010\\t5.1436M.Rdata")
load("C:\\data\\all2010\\t5.1441A.Rdata")




t5.9206B$Flow.idx<-rep(1:288,dim(table(t5.9206B$Date)))
t5.6013B$Flow.idx<-rep(1:288,dim(table(t5.6013B$Date)))
t5.9188A$Flow.idx<-rep(1:288,dim(table(t5.9188A$Date)))
t5.1431A$Flow.idx<-rep(1:288,dim(table(t5.1431A$Date)))


#--------------------------------------------------------------------------------------------------#
#-----------------------------------------SPLINE FLOWS---------------------------------------------#
#--------------------------------------------------------------------------------------------------#

# Data Selection:
cond.func<-function(obj){
#cond<-obj$Hour %in% 6:20 & obj$Dweek=="Mon" & months(obj$Date1) %in% c("Feb","Mar","Apr","May","Jun") & is.na(obj$Flow)==F
cond<-obj$Hour %in% 6:20 & obj$Dweek=="Wed" & months(obj$Date1) %in% c("Jul","Ago") & is.na(obj$Flow)==F
return(cond)}
cond1<-cond.func(t5.9206B);cond2<-cond.func(t5.6013B);cond3<-cond.func(t5.9188A);cond4<-cond.func(t5.1431A)

# Spline Function:
spl.flow<-function(flows,knts){
  sp.Flow<-with(flows,ns(Flow.idx,knots=knts))
  aux.plot<-data.frame(Flow.idx=flows$Flow.idx,
          Flow.pred=lm(Flow~sp.Flow[1:nrow(sp.Flow),],data=flows)$fitted.values)
  aux.plot<-aux.plot[order(aux.plot$Flow.idx),]
  axis(1,at=seq(6*12+1,23*12,by=12),lab=6:22)
  with(aux.plot,lines(Flow.idx,Flow.pred,col="red",lwd=2))}

# Flows + Splines:
par(mfrow=c(2,1))
with(t5.9206B[cond1,],plot(Flow~Flow.idx,type="p",col="dark grey",xaxt="n",pch=19,
                           main="Flows 9206B (prior data)",ylim=c(0,800)))
spl.flow(t5.9206B[cond1,],c(seq(40,100,by=5),seq(180,220,by=10)))

with(t5.6013B[cond2,],plot(Flow~Flow.idx,type="p",col="dark grey",xaxt="n",pch=19,
                           main="Flows 6013B (prior data)",ylim=c(0,800)))
spl.flow(t5.6013B[cond2,],c(seq(60,130,by=10),seq(180,220,by=10)))

with(t5.9188A[cond3,],plot(Flow~Flow.idx,type="p",col="dark grey",xaxt="n",pch=19,
                           main="Flows 9188A (prior data)",ylim=c(0,800)))
spl.flow(t5.9188A[cond3,],c(seq(60,130,by=10),seq(180,220,by=10)))

with(t5.1431A[cond4,],plot(Flow~Flow.idx,type="p",col="dark grey",xaxt="n",pch=19,
                           main="Flows 1431A (prior data)",ylim=c(0,800)))
spl.flow(t5.1431A[cond4,],c(seq(40,130,by=10),seq(180,220,by=10)))

#----------------------------------------------------------------------------------------------------#
#--------------------------------------SPLINE FLOW VS H/O/S------------------------------------------#
#----------------------------------------------------------------------------------------------------#
lag1<-function(vbl) c(NA,vbl[-length(vbl)])
t5.9206B$Speed_lag1<-with(t5.9206B,lag1(Speed));t5.6013B$Speed_lag1<-with(t5.6013B,lag1(Speed))
t5.9188A$Speed_lag1<-with(t5.9188A,lag1(Speed));t5.1431A$Speed_lag1<-with(t5.1431A,lag1(Speed))


cond.func1<-function(obj){
cond<-obj$Hour %in% 6:20 & obj$Dweek=="Fri" & months(obj$Date1) %in% c("Jul","Ago") & is.na(obj$Flow)==F & is.na(obj$Occ_lag1)==F & is.na(obj$Hway_lag1)==F & obj$Hway_lag1>1 &  is.na(obj$Speed_lag1)==F
return(cond)}
cond1<-cond.func1(t5.9206B);cond2<-cond.func1(t5.6013B);cond3<-cond.func1(t5.9188A);cond4<-cond.func1(t5.1431A)
cond1c<-cond.func1(t5.9200B);cond2c<-cond.func1(t5.6007L);cond3c<-cond.func1(t5.9193J);cond4c<-cond.func1(t5.1437A)
cond1g<-cond.func1(t5.9189B);cond2g<-cond.func1(t5.6004L);cond3g<-cond.func1(t5.1436M);cond4g<-cond.func1(t5.1441A)


#------------#
# Occupancy  #
#------------#


o.plot<-function(obj,node,ymax,xmax){
with(obj,plot(Occ_lag1,Flow,main=paste(node,"Flow vs Occupancy"),xlab="Occupancy at t-1",
              ylab="Flow",ylim=c(0,ymax),xlim=c(0,xmax)))


spl.occ<-with(obj,ns(Occ_lag1,df=4))
flow.pred.occ<-lm(Flow~spl.occ,data=obj)
aux<-data.frame(Occ=obj$Occ_lag1,Flow=flow.pred.occ$fitted.values)
aux<-aux[order(aux$Occ),]
lines(aux$Occ,aux$Flow,col="red",lwd=2) }

ymax<-400
xmax<-35
o.plot(t5.9200B[cond1c,],"9200B",ymax,xmax);o.plot(t5.6007L[cond2c,],"6007L",ymax)
o.plot(t5.9193J[cond3c,],"9193J",ymax); o.plot(t5.1437A[cond4c,],"1437A",ymax)


#ymax<-500
o.plot(t5.9189B[cond1g,],"9189B",ymax);o.plot(t5.6004L[cond2g,],"6004L",ymax)
o.plot(t5.1436M[cond3g,],"1436M",ymax); o.plot(t5.1441A[cond4g,],"1441A",ymax)
#test:


#----------#
# Headway  #
#----------#
h.plot<-function(obj,node){
with(obj,plot(Hway_lag1,Flow,xlim=c(0,15),ylim=c(0,800),main=paste(node,"Flow vs Headway"),xlab="Headway at t-1",ylab="Flow"))
spl.Hway<-with(obj,ns(Hway_lag1,df=4))
flow.pred.Hway<-lm(Flow~spl.Hway,data=obj)
aux<-data.frame(Hway=obj$Hway_lag1,Flow=flow.pred.Hway$fitted.values)
aux<-aux[order(aux$Hway),]
lines(aux$Hway,aux$Flow,col="red",lwd=2)

}

#----------#
# Speed    #
#----------#


s.plot<-function(obj,node){

with(obj,plot(Speed_lag1,Flow,xlim=c(0,150),ylim=c(0,800),main=paste(node,"Flow vs Speed"),xlab="Speed at t-1",ylab="Flow"))
spl.Speed<-with(obj,ns(Speed_lag1,df=4,intercept =TRUE))
flow.pred.Speed<-lm(Flow~spl.Speed,data=obj)
aux<-data.frame(Speed=obj$Speed_lag1,Flow=flow.pred.Speed$fitted.values)
aux<-aux[order(aux$Speed),]
lines(aux$Speed,aux$Flow,col="red",lwd=2)

}


par(mfrow=c(1,3))

ymax<-800
xmax<-27
#Fri
o.plot(t5.9206B[cond1,],"9206B",ymax,xmax);o.plot(t5.6013B[cond2,],"6013B",ymax,xmax)
o.plot(t5.9188A[cond3,],"9188A",ymax,xmax); o.plot(t5.1431A[cond4,],"1431A",ymax,xmax)
ymax<-450
o.plot(t5.9200B[cond1c,],"9200B",ymax,xmax);o.plot(t5.6007L[cond2c,],"6007L",ymax,xmax)
o.plot(t5.9193J[cond3c,],"9193J",ymax,xmax); o.plot(t5.1437A[cond4c,],"1437A",ymax,xmax)
ymax<-400
o.plot(t5.9189B[cond1g,],"9189B",ymax,xmax);o.plot(t5.6004L[cond2g,],"6004L",ymax,xmax)
o.plot(t5.1436M[cond3g,],"1436M",ymax,xmax); o.plot(t5.1441A[cond4g,],"1441A",ymax,xmax)
dev.copy2pdf(file="C:\\plots\\Occ_all.pdf")




#Fri
#par(mfrow=c(3,4))

h.plot(t5.9206B[cond1,],"9206B");h.plot(t5.6013B[cond2,],"6013B")
h.plot(t5.9188A[cond3,],"9188A");h.plot(t5.1431A[cond4,],"1431A")
ymax<-450
h.plot(t5.9200B[cond1c,],"9200B");h.plot(t5.6007L[cond2c,],"6007L")
h.plot(t5.9193J[cond3c,],"9193J"); h.plot(t5.1437A[cond4c,],"1437A")
ymax<-400
h.plot(t5.9189B[cond1g,],"9189B");h.plot(t5.6004L[cond2g,],"6004L")
h.plot(t5.1436M[cond3g,],"1436M"); h.plot(t5.1441A[cond4g,],"1441A")
dev.copy2pdf(file="C:\\plots\\Hway_all.pdf")

#Wed
s.plot(t5.9206B[cond1,],"9206B");s.plot(t5.6013B[cond2,],"6013B")
s.plot(t5.9188A[cond3,],"9188A");s.plot(t5.1431A[cond4,],"1431A")
s.plot(t5.9200B[cond1c,],"9200B");s.plot(t5.6007L[cond2c,],"6007L")
s.plot(t5.9193J[cond3c,],"9193J"); s.plot(t5.1437A[cond4c,],"1437A")
s.plot(t5.9189B[cond1g,],"9189B");s.plot(t5.6004L[cond2g,],"6004L")
s.plot(t5.1436M[cond3g,],"1436M"); s.plot(t5.1441A[cond4g,],"1441A")
dev.copy2pdf(file="C:\\plots\\Speed_all.pdf")


# Speed vs concentration:

par(mfrow=c(1,4))
with(t5.9206B[cond1,],plot(Occ_lag1, Speed,main="9206B Speed vs Occupancy", ylab="Speed @ t", xlab="Occupancy @ t",col="grey",pch=16,xlim=c(0,70)))
with(t5.6013B[cond2,],plot(Occ_lag1, Speed,main="6013B Speed vs Occupancy", ylab="Speed @ t", xlab="Occupancy @ t",col="grey",pch=16,xlim=c(0,70)))
with(t5.9188A[cond3,],plot(Occ_lag1, Speed,main="9188A Speed vs Occupancy", ylab="Speed @ t", xlab="Occupancy @ t",col="grey",pch=16,xlim=c(0,70)))
with(t5.1431A[cond4,],plot(Occ_lag1, Speed,main="1431A Speed vs Occupancy", ylab="Speed @ t", xlab="Occupancy @ t",col="grey",pch=16,xlim=c(0,70)))


#Flow vs Speed/occ


cnd<-t5.9188A$Hour %in% 6:20 & t5.9188A$Dweek=="Wed" & months(t5.9188A$Date1) %in% c("Apr","May","Jun")  &
    is.na(t5.9188A$Flow)==F & is.na(t5.9188A$Speed_lag1)==F & is.na(t5.9188A$Occ_lag1)==F

plot.SC<-function(obj,node,y.lim){
obj$Speed_lag1<-with(obj,lag1(Speed))
cnd<-obj$Hour %in% 6:20 & obj$Dweek=="Wed" & months(obj$Date1) %in% c("Apr","May","Jun") & is.na(obj$Flow)==F & is.na(obj$Flow_lag1)==F & is.na(obj$Speed_lag1)==F & is.na(obj$Occ_lag1)==F
spl.Speed.Occ<-with(obj[cnd,],ns(Speed_lag1*Occ_lag1,df=4))
flow.pred.Speed<-lm(Flow~spl.Speed.Occ,data=obj[cnd,])
aux<-data.frame(Speed.Occ=obj[cnd,]$Speed_lag1*obj[cnd,]$Occ_lag1,Flow=flow.pred.Speed$fitted.values)
aux<-aux[order(aux$Speed.Occ),]

#print(with(obj[cnd,],table(is.na(Occ_lag1*Speed_lag1))))
with(obj[cnd,],plot(Speed_lag1*Occ_lag1,Flow,ylim=y.lim,xlim=c(0,2500),
              main=paste(node,"Flow vs Speed*Occ"),xlab="Occ*Speed at t-1",ylab="Flow",pch=19,col="dark grey",type="n"))
              
cnd<-obj$Hour %in% 15:18 & obj$Dweek=="Wed" & months(obj$Date1) %in% c("Apr","May","Jun")
with(obj[cnd ,],points(Speed_lag1*Occ_lag1,Flow,pch=19,col="green"))
cnd<-obj$Hour %in% c(6:14,19:20) & obj$Dweek=="Wed" & months(obj$Date1) %in% c("Apr","May","Jun")
with(obj[cnd ,],points(Speed_lag1*Occ_lag1,Flow,pch=19,col="dark grey"))

#cnd<-obj$Hour %in% 6:10 & obj$Dweek=="Wed" & months(obj$Date1) %in% c("Apr","May","Jun")
#with(obj[cnd ,],points(Speed_lag1*Occ_lag1,Flow,pch=19,col="red"))
lines(aux$Speed.Occ,aux$Flow,col="blue",lwd=2)


}
par(mfrow=c(1,2))
plot.SC(t5.6013B,"6013B",y.lim=c(0,500));plot.SC(t5.6007L,"6007L",y.lim=c(0,500))



par(mfrow=c(2,4))
#par(mfrow=c(1,1))
plot.SC(t5.9206B,"9206B",y.lim=c(0,800));


plot.SC(t5.6013B,"6013B",y.lim=c(0,500));plot.SC(t5.9188A,"9188A",y.lim=c(0,600));plot.SC(t5.1431A,"1431A",y.lim=c(0,600))
plot.SC(t5.9200B,"9200B",y.lim=c(0,500));plot.SC(t5.6007L,"6007L",y.lim=c(0,250));plot.SC(t5.9193J,"9193J",y.lim=c(0,300));plot.SC(t5.1437A,"1437A",y.lim=c(0,400))


#xlim=c(0,150),

spl.Speed.Occ<-with(t5.9188A[cnd,],ns(Speed_lag1*Occ_lag1,df=4))
flow.pred.Speed<-lm(Flow~spl.Speed.Occ,data=t5.9188A[cnd,])
aux<-data.frame(Speed.Occ=t5.9188A[cnd,]$Speed_lag1*t5.9188A[cnd,]$Occ_lag1,Flow=flow.pred.Speed$fitted.values)
aux<-aux[order(aux$Speed.Occ),]
lines(aux$Speed.Occ,aux$Flow,col="blue",lwd=2)

#--------------------------#
# Comparing month by month #
#--------------------------#
library(lattice)
cond.func1<-function(obj,wkday="Wed"){
#months(obj$Date1) %in% mnth & obj$Dweek==wkday &
cond<-obj$Hour %in% 6:20 &  is.
na(obj$Flow)==F & is.na(obj$Occ_lag1)==F & is.na(obj$Hway_lag1)==F & obj$Hway_lag1>1 &  is.na(obj$Speed_lag1)==F
return(cond)}
cond1<-cond.func1(t5.9206B)
cond2<-cond.func1(t5.6013B)
cond3<-cond.func1(t5.9188A)
cond4<-cond.func1(t5.1431A)

with(t5.9206B[cond1,],xyplot(Flow~Hway_lag1|as.factor(months(Date1)),pch=19,col=1,as.table=T,main="Hway 9206B"))
dev.copy2pdf(file="C:\\plots\\Hway-9206B.pdf")

with(t5.6013B[cond1,],xyplot(Flow~Hway_lag1|as.factor(months(Date1)),pch=19,col=1,as.table=T,main="Hway 6013B"))
dev.copy2pdf(file="C:\\plots\\Hway-6013B.pdf")

with(t5.9188A[cond1,],xyplot(Flow~Hway_lag1|as.factor(months(Date1)),pch=19,col=1,as.table=T,main="Hway 9188A"))
dev.copy2pdf(file="C:\\plots\\Hway-9188A.pdf")

with(t5.1431A[cond1,],xyplot(Flow~Hway_lag1|as.factor(months(Date1)),pch=19,col=1,as.table=T,main="Hway 1431A"))
dev.copy2pdf(file="C:\\plots\\Hway-1431A.pdf")



#--------------------------------------------------------------------------------------------------#
#--------------------------------------SPLINE PROPORTIONS------------------------------------------#
#--------------------------------------------------------------------------------------------------#
load("C:\\data\\all2010\\t5.9200B.Rdata");load("C:\\data\\all2010\\t5.9189B.Rdata")
load("C:\\data\\all2010\\t5.6007L.Rdata");load("C:\\data\\all2010\\t5.6004L.Rdata")
load("C:\\data\\all2010\\t5.9193J.Rdata");load("C:\\data\\all2010\\t5.1436M.Rdata")
load("C:\\data\\all2010\\t5.1437A.Rdata");load("C:\\data\\all2010\\t5.1441A.Rdata")

t5.9200B$Flow.idx<-rep(1:288,dim(table(t5.9200B$Date)));t5.9189B$Flow.idx<-rep(1:288,dim(table(t5.9189B$Date)))
t5.6007L$Flow.idx<-rep(1:288,dim(table(t5.6007L$Date)));t5.6004L$Flow.idx<-rep(1:288,dim(table(t5.6004L$Date)))
t5.9193J$Flow.idx<-rep(1:288,dim(table(t5.9193J$Date)));t5.1436M$Flow.idx<-rep(1:288,dim(table(t5.1436M$Date)))
t5.1437A$Flow.idx<-rep(1:288,dim(table(t5.1437A$Date)));t5.1441A$Flow.idx<-rep(1:288,dim(table(t5.1441A$Date)))

# Data Selection:
c.func2<-function(obj){
cond<-obj$Hour %in% 6:22 & months(obj$Date1) %in% c("Feb","Mar","Apr","May","Jun") & obj$Dweek=="Mon"
return(cond)}
cond1<-c.func2(t5.9206B);cond2<-c.func2(t5.6013B);cond3<-c.func2(t5.9188A);cond4<-c.func2(t5.1431A)
c.chd1<-c.func2(t5.9200B);c.chd2<-c.func2(t5.6007L);c.chd3<-c.func2(t5.9193J);c.chd4<-c.func2(t5.1437A)
c.gchd1<-c.func2(t5.9189B);c.gchd2<-c.func2(t5.6004L);c.gchd3<-c.func2(t5.1436M);c.gchd4<-c.func2(t5.1441A)


# Proportion box plots:
n<-length(t5.9200B[c.chd1,]$Flow)
aux.boxplot<-data.frame(
rootnode=sort(rep(c("(1) 9206B to 9200B","(2) 9206B to 9189B",
                    "(3) 6013B to 6007L","(4) 6013B to 6004L",
                    "(5) 9188A to 9193J","(6) 9188A to 1436M",
                    "(7) 1431A to 1437A","(8) 1431A to 1441A"),n)),
prop=c(t5.9200B[c.chd1,]$Flow/t5.9206B[cond1,]$Flow,t5.9189B[c.gchd1,]$Flow/t5.9206B[cond1,]$Flow,
       t5.6007L[c.chd2,]$Flow/t5.6013B[cond2,]$Flow,t5.6004L[c.gchd2,]$Flow/t5.6013B[cond2,]$Flow,
       t5.9193J[c.chd3,]$Flow/t5.9188A[cond3,]$Flow,t5.1436M[c.gchd3,]$Flow/t5.9188A[cond3,]$Flow,
       t5.1437A[c.chd4,]$Flow/t5.1431A[cond4,]$Flow,t5.1441A[c.gchd4,]$Flow/t5.1431A[cond4,]$Flow))


par(mfrow=c(1,1))
par(mar=c(5.1, 4.1+5, 4.1, 2.1))
with(aux.boxplot,boxplot(prop~rootnode,col="dark grey",boxwex=0.5,outline=F,ylim=c(0,1),horizontal=T,yaxt="n",
              main="Flow proportions from root node to descendants",cex.axis=1.2,border="black"))
abline(h=c(2.5,4.5,6.5),col="dark grey")
lab.hr<-c("9206B to 9200B","9206B to 9189B","6013B to 6007L","6013B to 6004L",
          "9188A to 9193J","9188A to 1436M","1431A to 1437A","1431A to 1441A")
text(-0.05,1:8, adj = 1,          labels = lab.hr, xpd = TRUE,cex=1.2)
par(mar=c(5.1, 4.1, 4.1, 2.1))


# Function para TS proportion & spline

prop.fc<-function(rn.obj,rn,desc.obj,desc,y.lim){
  aux.prop<-data.frame(Flow.idx=desc.obj$Flow.idx,Flow.par=rn.obj$Flow,Flow.chd=desc.obj$Flow,Flow.prop=desc.obj$Flow/rn.obj$Flow)
  aux.prop<-aux.prop[is.na(aux.prop$Flow.prop)==F,]
  with(aux.prop,plot(Flow.idx,Flow.prop,xaxt="n",col="dark grey",pch=19,ylim=y.lim,
                     main=paste("Flow proportion:",rn," to", desc),ylab="proportion",xlab="Time"))
  axis(1,at=seq(73,276,by=12),lab=rep(6:22))
  sp.prop<-with(aux.prop,ns(Flow.idx,knots=c(seq(60,120,by=10),seq(180,220,by=10))))
  aux.plot<-data.frame(Flow.idx=aux.prop$Flow.idx,
            Flow.pred=lm(Flow.prop~sp.prop[1:nrow(sp.prop),],data=aux.prop)$fitted.values)
  aux.plot<-aux.plot[order(aux.plot$Flow.idx),]
  with(aux.plot,lines(Flow.idx,Flow.pred,col="red",lwd=2))
  }

#       ROOT NODE TO CHILDREN       #
par(mfrow=c(2,1))
prop.fc(rn.obj=t5.9206B[cond1,],rn="9206B",desc.obj=t5.9200B[c.chd1,],desc="9200B",y.lim=c(0.3,0.8))
prop.fc(rn.obj=t5.6013B[cond2,],rn="6013B",desc.obj=t5.6007L[c.chd2,],desc="6007L",y.lim=c(0.2,0.9))
prop.fc(rn.obj=t5.9188A[cond3,],rn="9188A",desc.obj=t5.9193J[c.chd3,],desc="9193J",y.lim=c(0,0.5))
prop.fc(rn.obj=t5.1431A[cond4,],rn="1431A",desc.obj=t5.1437A[c.chd4,],desc="1437A",y.lim=c(0.3,1))



#       CHILDREN TO grandchildren       #

par(mfrow=c(2,1))
prop.fc(rn.obj=t5.9200B[c.chd1,],rn="9200B",desc.obj=t5.9189B[c.gchd1,],desc="9189B",y.lim=c(0.7,1))
prop.fc(rn.obj=t5.6007L[c.chd2,],rn="6007L",desc.obj=t5.6004L[c.gchd2,],desc="6004L",y.lim=c(0.0,0.7))
prop.fc(rn.obj=t5.9193J[c.chd3,],rn="9193J",desc.obj=t5.1436M[c.gchd3,],desc="1436M",y.lim=c(0.2,1))
prop.fc(rn.obj=t5.1437A[c.chd4,],rn="1437A",desc.obj=t5.1441A[c.gchd4,],desc="1441A",y.lim=c(0.5,1))

