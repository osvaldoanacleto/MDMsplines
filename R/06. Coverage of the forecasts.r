
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







r<--c(1:360)


c.par1<-(months(t5.9206B$Date1) %in% learning.months | t5.9206B$Date %in% test.dates)  &  t5.9206B$Hour %in% time.filter  & t5.9206B$Dweek==wkday
c.par2<-(months(t5.6013B$Date1) %in% learning.months | t5.6013B$Date %in% test.dates)  & t5.6013B$Hour %in% time.filter  & t5.6013B$Dweek==wkday
c.par3<-(months(t5.9188A$Date1) %in% learning.months | t5.9188A$Date %in% test.dates)  & t5.9188A$Hour %in% time.filter  & t5.9188A$Dweek==wkday
c.par4<-(months(t5.1431A$Date1) %in% learning.months | t5.1431A$Date %in% test.dates)  & t5.1431A$Hour %in% time.filter  & t5.1431A$Dweek==wkday


c.chld1<-(months(t5.9200B$Date1) %in% learning.months | t5.9200B$Date %in% test.dates) & t5.9200B$Hour %in% time.filter  & t5.9200B$Dweek==wkday
c.chld2<-(months(t5.6007L$Date1) %in% learning.months | t5.6007L$Date %in% test.dates) & t5.6007L$Hour %in% time.filter  & t5.6007L$Dweek==wkday
c.chld3<-(months(t5.9193J$Date1) %in% learning.months | t5.9193J$Date %in% test.dates) & t5.9193J$Hour %in% time.filter  & t5.9193J$Dweek==wkday
c.chld4<-(months(t5.1437A$Date1) %in% learning.months | t5.1437A$Date %in% test.dates) & t5.1437A$Hour %in% time.filter  & t5.1437A$Dweek==wkday


c.gchld1<-(months(t5.9189B$Date1) %in% learning.months | t5.9189B$Date %in% test.dates) & t5.9189B$Hour %in% time.filter  & t5.9189B$Dweek==wkday
c.gchld2<-(months(t5.6004L$Date1) %in% learning.months | t5.6004L$Date %in% test.dates) & t5.6004L$Hour %in% time.filter  & t5.6004L$Dweek==wkday
c.gchld3<-(months(t5.1436M$Date1) %in% learning.months | t5.1436M$Date %in% test.dates) & t5.1436M$Hour %in% time.filter  & t5.1436M$Dweek==wkday
c.gchld4<-(months(t5.1441A$Date1) %in% learning.months | t5.1441A$Date %in% test.dates) & t5.1441A$Hour %in% time.filter  & t5.1441A$Dweek==wkday
cov.func.spl<-function(node,flow,model){
  inf.n<-model[["f"]][r]-1.96*sqrt(model[["Q"]][r])
  sup.n<-model[["f"]][r]+1.96*sqrt(model[["Q"]][r])
  aux.new<-data.frame(out.inf=flow<inf.n,out.sup=flow>sup.n)
  aux.new$out.total=aux.new$out.inf+aux.new$out.sup
  total<-table(is.na(flow))[1]
    ans<-rep(NA,4)
  ans[1:3]<-with(aux.new,c(sum(out.inf,na.rm=T),sum(out.sup,na.rm=T),sum(out.total,na.rm=T)))
  ans[4]<-round(1-ans[3]/total,2)
  return(ans)
}

mean(cov.func.spl("9206B",t5.9206B[c.par1,]$Flow[r],dlm.9206B[[9]])[4],
cov.func.spl("6013B",t5.6013B[c.par2,]$Flow[r],dlm.6013B[[9]])[4],
cov.func.spl("9188A",t5.9188A[c.par3,]$Flow[r],dlm.9188A[[9]])[4],
cov.func.spl("1431A",t5.1431A[c.par4,]$Flow[r],dlm.1431A[[9]])[4])



mean(cov.func.spl("9200B",t5.9200B[c.chld1,]$Flow[r],MDM.9200B[[5]])[4],
cov.func.spl("6007L",t5.6007L[c.chld2,]$Flow[r],MDM.6007L[[5]])[4],
cov.func.spl("9193J",t5.9193J[c.chld3,]$Flow[r],MDM.9193J[[5]])[4],
cov.func.spl("1437A",t5.1437A[c.chld4,]$Flow[r],MDM.1437A[[5]])[4])




mean(cov.func.spl("9189B",t5.9189B[c.gchld1,]$Flow[r],MDM.9189B[[5]])[4],
cov.func.spl("6004L",t5.6004L[c.gchld2,]$Flow[r],MDM.6004L[[5]])[4],
cov.func.spl("1436M",t5.1436M[c.gchld3,]$Flow[r],MDM.1436M[[5]])[4],
cov.func.spl("1441A",t5.1441A[c.gchld4,]$Flow[r],MDM.1441A[[5]])[4])



covs<-list()
covs[["RN"]]<-data.frame(daily=rep(NA,4),full=rep(NA,4))
row.names(covs[["RN"]])<-c("9206B","6013B","9188A","1431A")

covs[["child"]]<-data.frame(daily=rep(NA,4),full=rep(NA,4))
row.names(covs[["child"]])<-c("9200B","6007L","9193J","1437A")

covs[["Gchild"]]<-data.frame(daily=rep(NA,4),full=rep(NA,4))
row.names(covs[["Gchild"]])<-c("9189B","6004L","1436M","1441A")


covs[["RN"]]$full<-c(cov.func.spl("9206B",t5.9206B[c.par1,]$Flow[r],dlm.9206B[[9]])[4],
cov.func.spl("6013B",t5.6013B[c.par2,]$Flow[r],dlm.6013B[[9]])[4],
cov.func.spl("9188A",t5.9188A[c.par3,]$Flow[r],dlm.9188A[[9]])[4],
cov.func.spl("1431A",t5.1431A[c.par4,]$Flow[r],dlm.1431A[[9]])[4])

covs[["RN"]]$daily<-c(cov.func.spl("9206B",t5.9206B[c.par1,]$Flow[r],dlm.9206B[[1]])[4],
cov.func.spl("6013B",t5.6013B[c.par2,]$Flow[r],dlm.6013B[[1]])[4],
cov.func.spl("9188A",t5.9188A[c.par3,]$Flow[r],dlm.9188A[[1]])[4],
cov.func.spl("1431A",t5.1431A[c.par4,]$Flow[r],dlm.1431A[[1]])[4])




covs[["child"]]$full<-c(cov.func.spl("9200B",t5.9200B[c.chld1,]$Flow[r],MDM.9200B[[5]])[4],
cov.func.spl("6007L",t5.6007L[c.chld2,]$Flow[r],MDM.6007L[[5]])[4],
cov.func.spl("9193J",t5.9193J[c.chld3,]$Flow[r],MDM.9193J[[5]])[4],
cov.func.spl("1437A",t5.1437A[c.chld4,]$Flow[r],MDM.1437A[[5]])[4])

covs[["child"]]$daily<-c(cov.func.spl("9200B",t5.9200B[c.chld1,]$Flow[r],MDM.9200B[[1]])[4],
cov.func.spl("6007L",t5.6007L[c.chld2,]$Flow[r],MDM.6007L[[1]])[4],
cov.func.spl("9193J",t5.9193J[c.chld3,]$Flow[r],MDM.9193J[[1]])[4],
cov.func.spl("1437A",t5.1437A[c.chld4,]$Flow[r],MDM.1437A[[1]])[4])



covs[["Gchild"]]$full<-c(cov.func.spl("9189B",t5.9189B[c.gchld1,]$Flow[r],MDM.9189B[[5]])[4],
cov.func.spl("6004L",t5.6004L[c.gchld2,]$Flow[r],MDM.6004L[[5]])[4],
cov.func.spl("1436M",t5.1436M[c.gchld3,]$Flow[r],MDM.1436M[[5]])[4],
cov.func.spl("1441A",t5.1441A[c.gchld4,]$Flow[r],MDM.1441A[[5]])[4])

covs[["Gchild"]]$daily<-c(cov.func.spl("9189B",t5.9189B[c.gchld1,]$Flow[r],MDM.9189B[[1]])[4],
cov.func.spl("6004L",t5.6004L[c.gchld2,]$Flow[r],MDM.6004L[[1]])[4],
cov.func.spl("1436M",t5.1436M[c.gchld3,]$Flow[r],MDM.1436M[[1]])[4],
cov.func.spl("1441A",t5.1441A[c.gchld4,]$Flow[r],MDM.1441A[[1]])[4])


