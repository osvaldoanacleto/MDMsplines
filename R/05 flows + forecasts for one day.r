plot.day.fc<-function(node,Flow.obj,model,y.lim,model.type,Date,ind.A,mA.leg,ind.B,mB.leg){
  time.filter<-6:20
  test.dates<-c("18/08/2010","25/08/2010")
  learning.months<-c("Sep","Oct")
  wkday<-"Wed"
  d<-substr(Date,1,2)
  m<-substr(Date,4,5)
  y<-substr(Date,7,10)

  aux.cond<-(months(Flow.obj$Date1) %in% learning.months | Flow.obj$Date %in% test.dates)  & Flow.obj$Hour %in% time.filter  & Flow.obj$Dweek==wkday
  if (model.type=="DLM"){
   aux<-data.frame(Date=Flow.obj[aux.cond,]$Date,Flow=Flow.obj[aux.cond,]$Flow,modA.f=model[[ind.A]]$f,modB.f=model[[ind.B]]$f,
                  modA.var=model[[ind.A]]$Q,modB.var=model[[ind.B]]$Q)}
  else
   {aux<-data.frame(Date=Flow.obj[aux.cond,]$Date,Flow=Flow.obj[aux.cond,]$Flow,modA.f=model[[ind.A]]$f.marg,modB.f=model[[ind.B]]$f.marg,
                    modA.var=model[[ind.A]]$marg.var,modB.var=model[[ind.B]]$marg.var)}
 cond<-aux$Date==Date
 pdf(file=paste("/home/terrahome/oaj9/plots/",node,"-",paste(d,"-",m,sep=""),".pdf"),width=15)
 with(aux[cond,],plot(Flow,pch=19,ylim=y.lim,xaxt="n",ylab="",xlab="Time",
                      main=paste(node,": flows and forecasts on ",paste(d,"-",m,"-",y,sep=""),sep="")))

 abline(h=seq(0,y.lim[2],by=50),col=rgb(0,0,0,alpha=0.1),lty=2)
 abline(v=seq(1,period,by=12),col=rgb(0,0,0,alpha=0.1),lty=2)
 axis(1,at=seq(1,period,by=12),lab=time.filter)
 
 with(aux[cond,],lines(modA.f,col="red"))
 with(aux[cond,], lines(modA.f+1.64*sqrt(modA.var),col="red",lty=2,lwd=2))
 with(aux[cond,], lines(modA.f-1.64*sqrt(modA.var),col="red",lty=2,lwd=2))

 with(aux[cond,],lines(modB.f,col="Blue"))
 with(aux[cond,], lines(modB.f+1.64*sqrt(modB.var),col="Blue",lty=2,lwd=2))
 with(aux[cond,], lines(modB.f-1.64*sqrt(modB.var),col="Blue",lty=2,lwd=2))
 inc<-30
 ini<-11
 max.leg<-150
 legend((ini-6)*12,max.leg,bty="n",pch=19,legend="observed flows")
 legend((ini-6)*12,max.leg-inc,bty="n",col="red",lty=1,legend=paste("forecast means: ",mA.leg),cex=1,lwd=2)
 legend((ini-6)*12,max.leg-2*inc,bty="n",col="red",lty=2,legend=paste("forecast limits:",mA.leg),cex=1,lwd=2)
 legend((ini-6)*12,max.leg-3*inc,bty="n",col="blue",lty=1,legend=paste("forecast means:",mB.leg),cex=1,lwd=2)
 legend((ini-6)*12,max.leg-4*inc,bty="n",col="blue",lty=2,legend=paste("forecast limits:",mB.leg),cex=1,lwd=2)

 dev.off()

  }
  

#9200B
plot.day.fc("9200B",t5.9200B,MDM.9200B,c(0,500),"MDM",Date="27/10/2010",
            ind.A=1,mA.leg="daily cycle model",ind.B=4,mB.leg="full model")


  

#1431A
plot.day.fc("1431A",t5.1431A,dlm.1431A,c(0,800),"DLM",Date="20/10/2010",
            ind.A=1,mA.leg="daily cycle model",ind.B=9,mB.leg="full model")

  
#9206B
plot.day.fc("9206B",t5.9206B,dlm.9206B,c(0,1000),"DLM",Date="27/10/2010",
            ind.A=1,mA.leg="daily cycle model",ind.B=4,mB.leg="full model")
            
#9193J
plot.day.fc("9193J",t5.9193J,MDM.9193J,c(0,300),"MDM",Date="03/11/2010",
            ind.A=1,mA.leg="daily cycle model",ind.B=4,mB.leg="full model")



#9193J
plot.day.fc("9193J",t5.9193J,MDM.9193J,c(0,300),"MDM",Date="03/11/2010",
            ind.A=1,mA.leg="daily cycle model",ind.B=4,mB.leg="full model")



#1441A
plot.day.fc("1441A",t5.1441A,MDM.1441A,c(0,500),"MDM",Date="17/11/2010",
            ind.A=1,mA.leg="daily cycle model",ind.B=4,mB.leg="full model")
