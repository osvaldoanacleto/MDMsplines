load("/home/terrahome/oaj9/splines/dlm.9206B.Rdata")
load("/home/terrahome/oaj9/splines/dlm.6013B.Rdata")
load("/home/terrahome/oaj9/splines/dlm.9188A.Rdata")
load("/home/terrahome/oaj9/splines/dlm.1431A.Rdata")
load("/home/terrahome/oaj9/splines/MDM.9200B.Rdata")
load("/home/terrahome/oaj9/splines/MDM.6007L.Rdata")
load("/home/terrahome/oaj9/splines/MDM.9193J.Rdata")
load("/home/terrahome/oaj9/splines/MDM.1437A.Rdata")
load("/home/terrahome/oaj9/splines/MDM.9189B.Rdata")
load("/home/terrahome/oaj9/splines/MDM.6004L.Rdata")
load("/home/terrahome/oaj9/splines/MDM.1436M.Rdata")
load("/home/terrahome/oaj9/splines/MDM.1441A.Rdata")
load("/home/terrahome/oaj9/_data/all2010/t5.9206B.Rdata")
load("/home/terrahome/oaj9/_data/all2010/t5.6013B.Rdata")
load("/home/terrahome/oaj9/_data/all2010/t5.9188A.Rdata")
load("/home/terrahome/oaj9/_data/all2010/t5.1431A.Rdata")
load("/home/terrahome/oaj9/_data/all2010/t5.9200B.Rdata")
load("/home/terrahome/oaj9/_data/all2010/t5.6007L.Rdata")
load("/home/terrahome/oaj9/_data/all2010/t5.9193J.Rdata")
load("/home/terrahome/oaj9/_data/all2010/t5.1437A.Rdata")
load("/home/terrahome/oaj9/_data/all2010/t5.9189B.Rdata")
load("/home/terrahome/oaj9/_data/all2010/t5.6004L.Rdata")
load("/home/terrahome/oaj9/_data/all2010/t5.1436M.Rdata")
load("/home/terrahome/oaj9/_data/all2010/t5.1441A.Rdata")


   time.filter<-6:20
  period<-96
  learning.months<-c("Sep","Oct","Nov")
# Data Selection:
cond.func<-function(obj){
cond<-months(obj$Date1) %in% learning.months  & obj$Hour %in% time.filter  & obj$Dweek=="Wed"  & obj$Hour %in% 6:20
return(cond)}
cond1<-cond.func(t5.9206B);cond3<-cond.func(t5.9188A)
cond2<-cond.func(t5.6013B);cond4<-cond.func(t5.1431A)

ndays<-with(t5.9206B[cond1,],dim(table(Date1)))

plot.fc<-function(node,Flow.obj,model,y.lim,model.type){
 pdf(file=paste("/home/terrahome/oaj9/plots/",node,".pdf"),width=200)
 rr<--(1:360)
 plot(Flow.obj$Flow,pch=19,ylim=y.lim,xaxt="n")
  abline(v=seq(1,period*ndays,by=period/length(time.filter)),lty=1,col=rgb(0,0,0,alpha=0.1))
  abline(h=seq(0,y.lim[2],by=50),col=rgb(0,0,0,alpha=0.1),lty=2)
  abline(v=seq(1,period*ndays,by=period),lty=2,col="grey")
  axis(1,at=seq(1,period*ndays,by=period/length(time.filter)),lab=rep(time.filter,ndays))
  
  if (model.type=="DLM"){
  lines(model[[1]]$f[rr],col="red")
  lines(model[[1]]$f[rr]+1.64*sqrt(model[[1]]$Q[rr]),col="red",lty=2)
  lines(model[[1]]$f[rr]-1.64*sqrt(model[[1]]$Q[rr]),col="red",lty=2)
  lines(model[[9]]$f[rr],col="blue")
  lines(model[[9]]$f[rr]+1.64*sqrt(model[[9]]$Q[rr]),col="blue",lty=2)
  lines(model[[9]]$f[rr]-1.64*sqrt(model[[9]]$Q[rr]),col="blue",lty=2)}

  else
 {lines(model[[1]]$f.marg[rr],col="red")
  lines(model[[1]]$f.marg[rr]+1.64*sqrt(model[[1]]$marg.var[rr]),col="red",lty=2)
  lines(model[[1]]$f.marg[rr]-1.64*sqrt(model[[1]]$marg.var[rr]),col="red",lty=2)
  lines(model[[4]]$f.marg[rr],col="blue")
  lines(model[[4]]$f.marg[rr]+1.64*sqrt(model[[4]]$marg.var[rr]),col="blue",lty=2)
  lines(model[[4]]$f.marg[rr]-1.64*sqrt(model[[4]]$marg.var[rr]),col="blue",lty=2)}
  dev.off()
  }
  
#DLM:
plot.fc("9206B",t5.9206B[cond1,],dlm.9206B,c(0,1000),"DLM")
plot.fc("6013B",t5.6013B[cond2,],dlm.6013B,c(0,500),"DLM")
plot.fc("9188A",t5.9188A[cond3,],dlm.9188A,c(0,700),"DLM")
plot.fc("1431A",t5.1431A[cond4,],dlm.1431A,c(0,600),"DLM")



#MDM:

cond1c<-cond.func(t5.9200B);cond3c<-cond.func(t5.9193J)
cond2c<-cond.func(t5.6007L);cond4c<-cond.func(t5.1437A)
cond1g<-cond.func(t5.9189B);cond3g<-cond.func(t5.1436M)
cond2g<-cond.func(t5.6004L);cond4g<-cond.func(t5.1441A)


plot.fc("9200B",t5.9200B[cond1c,],MDM.9200B,c(0,500),"MDM")
plot.fc("6007L",t5.6007L[cond2c,],MDM.6007L,c(0,300),"MDM")
plot.fc("9193J",t5.9193J[cond3c,],MDM.9193J,c(0,300),"MDM")
plot.fc("1437A",t5.1437A[cond4c,],MDM.1437A,c(0,500),"MDM")

plot.fc("9189B",t5.9189B[cond1g,],MDM.9189B,c(0,500),"MDM")
plot.fc("6004L",t5.6004L[cond2g,],MDM.6004L,c(0,150),"MDM")
plot.fc("1436M",t5.1436M[cond3g,],MDM.1436M,c(0,700),"MDM")
plot.fc("1441A",t5.1441A[cond4g,],MDM.1441A,c(0,400),"MDM")


with(t5.1441A[cond4g,],as.matrix(table(Date1)))


rr<--(1:360)
summary(2*sqrt(MDM.1441A[[4]]$marg.var[rr]))






lines(MDM.9200B[[1]]$f.marg,col="red")
lines(MDM.9200B[[1]]$f.marg+2*sqrt(MDM.9200B[[1]]$marg.var),col="red",lty=2)
lines(MDM.9200B[[1]]$f.marg-2*sqrt(MDM.9200B[[1]]$marg.var),col="red",lty=2)

lines(MDM.9200B[[2]]$f.marg,col="blue")
lines(MDM.9200B[[2]]$f.marg+2*sqrt(MDM.9200B[[2]]$marg.var),col="blue",lty=2)
lines(MDM.9200B[[2]]$f.marg-2*sqrt(MDM.9200B[[2]]$marg.var),col="blue",lty=2)

#lines(MDM.9200B[[3]]$f.marg,col="green")
#lines(MDM.9200B[[3]]$f.marg+2*sqrt(MDM.9200B[[3]]$marg.var),col="green",lty=2)
#lines(MDM.9200B[[3]]$f.marg-2*sqrt(MDM.9200B[[3]]$marg.var),col="green",lty=2)


#lines(MDM.9200B[[4]]$f.marg,col="blue")
#lines(MDM.9200B[[4]]$f.marg+2*sqrt(MDM.9200B[[4]]$marg.var),col="blue",lty=2)
#lines(MDM.9200B[[4]]$f.marg-2*sqrt(MDM.9200B[[4]]$marg.var),col="blue",lty=2)

legend(15,160,bty="n",col="red",lty=1,"ML and parent with ML")
legend(15,130,bty="n",col="blue",lty=1,"ML and parent with ML+occ+way")
#legend(15,100,bty="n",col="green",lty=1,"ML+occ+way and parent with ML")
#legend(15,70,bty="n",col="orange",lty=1,"ML+occ+way and parent with ML+occ+way")

