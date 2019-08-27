

RN<-data.frame(measure=names(lpl.9206B)[-1],N9206B=as.numeric(lpl.9206B[5,-1]),N6013B=as.numeric(lpl.6013B[5,-1]),
      N9188A=as.numeric(lpl.9188A[5,-1]),N1431A=as.numeric(lpl.1431A[5,-1]))
write.csv(RN,file="/home/terrahome/oaj9/splines/RN.csv")

#RB-chd
RN.chd<-cbind(jlpl.9200B,jlpl.6007L[,2],jlpl.9193J[,2],jlpl.1437A[,2])
write.csv(RN.chd,file="/home/terrahome/oaj9/splines/RN-chd.csv")
#RB-chd-gchld
RN.chd.gchld<-cbind(jlpl.9189B,jlpl.6004L[,2],jlpl.1436M[,2],jlpl.1441A[,2])
write.csv(RN.chd.gchld,file="/home/terrahome/oaj9/splines/RN-chd-gchld.csv")