cond.func<-function(obj){
cond<-obj$Hour %in% 6:20 & obj$Dweek=="Fri"
return(cond)}
c1<-cond.func(t5.9206B);c2<-cond.func(t5.6013B);c3<-cond.func(t5.9188A);c4<-cond.func(t5.1431A)
c11<-cond.func(t5.9200B);c22<-cond.func(t5.6007L);c33<-cond.func(t5.9193J);c44<-cond.func(t5.1437A)
c111<-cond.func(t5.9189B);c222<-cond.func(t5.6004L);c333<-cond.func(t5.1436M);c444<-cond.func(t5.1441A)

plot.fc<-function(obj,node,x.lim){
with(obj,plot(1/Hway_lag1,Flow,
        main=paste(node,": 1/Hway vs Flow"),xlab="1/Hway @ t-1",ylab="Flow @ t",xlim=x.lim))
}

par(mfrow=c(2,2))
plot.fc(t5.9206B[c1,],"9206B",c(.05,.7))
plot.fc(t5.6013B[c2,],"6013B",c(.05,.5))
plot.fc(t5.9188A[c3,],"9188A",c(.05,.7))
plot.fc(t5.1431A[c4,],"1431A",c(.05,.5))



plot.fc(t5.9200B[c11,],"9200B",c(.05,.7))
plot.fc(t5.6007L[c22,],"6007L",c(.05,.35))
plot.fc(t5.9193J[c33,],"9193J",c(.05,.4))
plot.fc(t5.1437A[c4,],"1437A",c(.05,.6))
plot.fc(t5.9189B[c111,],"9189B",c(.05,.6))
plot.fc(t5.6004L[c222,],"6004L",c(.05,.25))
plot.fc(t5.1436M[c333,],"1436M",c(.05,.6))
plot.fc(t5.1441A[c444,],"1441A",c(.05,.5))




par(mfrow=c(1,1))
cond<-t5.1436M$Hour %in% 6:22 & t5.1436M$Dweek=="Mon"
with(t5.1436M[cond,],plot(Occupancy,Flow))



# Scatterplots

lag1<-function(vbl) c(NA,vbl[-length(vbl)])

par(mfrow=c(4,3))

kimelda<-function(node.obj,node){
node.obj$Speed_lag1<-with(node.obj,lag1(Speed))
par(mar=c(5.1-3, 4.1-2, 4.1-3, 2.1-2))
par(mgp=c(4,1,0))
cond.plot<-node.obj$Hour %in% 7:20 & !months(node.obj$Date1) %in% c("Mar","Apr","May","Jun")
           #node.obj$Dweek=="Mon"
with(node.obj[cond.plot,],plot(Occ_lag1,Flow,pch=19,main=paste(node,"F vs O")))
with(node.obj[cond.plot,],plot(Hway_lag1,Flow,pch=19,main=paste(node,"F vs H")))
with(node.obj[cond.plot,],plot(Speed_lag1,Flow,pch=19,main=paste(node,"F vs S")))
}
kimelda(t5.9206B,"9206B")      # Occ not good
kimelda(t5.6013B,"6013B")
kimelda(t5.9188A,"9188A")
kimelda(t5.1431A,"1431A")
dev.copy2pdf(file="C:\\plots\\Scatterplots A.pdf")


kimelda(t5.9200B,"9200B")      # Occ not good
kimelda(t5.6007L,"6007L")
kimelda(t5.9193J,"9193J")
kimelda(t5.1437A,"1437A")
dev.copy2pdf(file="C:\\plots\\Scatterplots B.pdf")


kimelda(t5.9189B,"9189B")      # Occ not good
kimelda(t5.6004L,"6004L")
kimelda(t5.1436M,"1436M")
kimelda(t5.1441A,"1441A")
dev.copy2pdf(file="C:\\plots\\Scatterplots C.pdf")
