#-------------------------------------------------------------#
#     Manchester Network                                      #
#     Reading Data from 2010 (just flows)                     #
#     Nodes:                                                  #
#     M62/1438B, M62/1436B, M602/6013B, M602/6002B, M602/6007L#
#     "M60/9206B","M60/9188A","M62/1431A"                     #
# (variance law/error analysis/root mod nodes                 #
#     (details on metadata.xls)                               #
#     Based on Casper's program convert_raw_data.R            #
#-------------------------------------------------------------#
library(car);library(chron);library(lattice)

# Reading all data:
all.2010.raw.a<-read.csv("/home/terrahome/oaj9/_data/all2010/01-03.txt", header = TRUE, as.is = c(1,6,7))
all.2010.raw.b<-read.csv("/home/terrahome/oaj9/_data/all2010/04-06.txt", header = TRUE, as.is = c(1,6,7))
all.2010.raw.c<-read.csv("/home/terrahome/oaj9/_data//all2010/07-09.txt", header = TRUE, as.is = c(1,6,7))
all.2010.raw.d<-read.csv("/home/terrahome/oaj9/_data/all2010/10-12.txt", header = TRUE, as.is = c(1,6,7))
all.2010.add<-read.csv("/home/terrahome/oaj9/_data/all2010/addnodes2010.txt", header = TRUE, as.is = c(1,6,7))
raw.9193J<-read.csv("/home/terrahome/oaj9/_data/all2010/9193J_not04-06.txt", header = TRUE, as.is = c(1,6,7))
all.2010.raw<-rbind(all.2010.raw.a,all.2010.raw.b,all.2010.raw.c,all.2010.raw.d,all.2010.add,raw.9193J)



library(chron)
all.2010.raw$Date1<-chron(dates = all.2010.raw$Date, format = "d/m/y",out.format="d/mon/y")

#all.2010.raw<-all.2010.raw[months(all.2010.raw$Date1)=="May",]

#checking if all nodes and dates are in the file (ok 26/01)
write.csv(table(all.2010.raw$Date1,all.2010.raw$GeogAddress),
          file="/home/terrahome/oaj9/delme/delme.csv")
#Excluding weekends:
#all.2010.raw<-all.2010.raw[is.weekend(all.2010.raw$Date1)==F,]
#Excluding Holidays:
all.2010.raw<-all.2010.raw[!all.2010.raw$Date %in% c("01/01/2010", "02/04/2010","05/04/2010","03/05/2010","31/05/2010",
                           "30/08/2010","25/12/2010","26/12/2010","31/12/2010"),]
all.2010.raw<-all.2010.raw[order(all.2010.raw$GeogAddress,all.2010.raw$Date,all.2010.raw$Time),]

all.2010.raw$Flow <-with(all.2010.raw,(Category_1_flow+Category_2_flow+Category_3_flow+Category_4_flow)/60)
all.2010.raw$Occupancy<-rowSums(all.2010.raw[,c("Lane_1_occupancy","Lane_2_occupancy","Lane_3_occupancy","Lane_4_occupancy")],na.rm=T)/
                     rowSums(is.na(all.2010.raw[,c("Lane_1_occupancy","Lane_2_occupancy",
                                                "Lane_3_occupancy","Lane_4_occupancy")])==F)
all.2010.raw$Headway<-rowSums(all.2010.raw[,c("Lane_1_headway","Lane_2_headway","Lane_3_headway","Lane_4_headway")],na.rm=T)/
                     rowSums(is.na(all.2010.raw[,c("Lane_1_headway","Lane_2_headway",
                                                "Lane_3_headway","Lane_4_headway")])==F)
#speed:
#sfi=speed x flow @ lane i

all.2010.raw$sf1<-with(all.2010.raw,Lane_1_average_speed*Lane_1_flow)
all.2010.raw$sf2<-with(all.2010.raw,Lane_2_average_speed*Lane_2_flow)
all.2010.raw$sf3<-with(all.2010.raw,Lane_3_average_speed*Lane_3_flow)
all.2010.raw$sf4<-with(all.2010.raw,Lane_4_average_speed*Lane_4_flow)

all.2010.raw$Speed<-rowSums(all.2010.raw[,c("sf1","sf2","sf3","sf4")],na.rm=T)/
                 rowSums(all.2010.raw[,c("Lane_1_flow","Lane_2_flow","Lane_3_flow","Lane_4_flow")],na.rm=T)                   
#table(all.2010.raw$GeogAddress)        
dim(table(all.2010.raw$Date))

length(all.2010.raw$Date)-254*19*1440-61*4*1440

all.2010.raw<-all.2010.raw
#save(all.2010.raw,file="/home/terrahome/oaj9/_data/all2010/all.2010.raw-WK.Rdata")
#removing aux files

load("/home/terrahome/oaj9/_data/all2010/all.2010.raw.Rdata")


delme<-with(all.2010.raw[all.2010.raw$GeogAddress=="M60/9200B",],aggregate(list(Flow),
                                             list(Hour), summary,na.rm=T))


summary(all.2010.raw$Speed)
summary(all.2010.raw$Occupancy)
summary(all.2010.raw$Headway)
summary(all.2010.raw$Flow)


table(is.na(all.2010.raw$Speed))/length(all.2010.raw$Date)
table(is.na(all.2010.raw$Occupancy))/length(all.2010.raw$Date)
table(is.na(all.2010.raw$Headway))/length(all.2010.raw$Date)
table(is.na(all.2010.raw$Flow))/length(all.2010.raw$Date)


#---------------------- ----------------------#
#    aggregating data into 5-min int          #
#---------------------------------------------#
load("/home/terrahome/oaj9/_data/all2010/all.2010.raw-WK.Rdata")
all.2010.raw<-all.2010.raw[months(all.2010.raw$Date1)=="May",]
str(all.2010.raw)
time.aux<-data.frame(Time=with(all.2010.raw,as.data.frame(table(Time))$Time))
for (i in 1:1440){
if (as.numeric(substr(time.aux$Time[i],5,5))<5)
  { time.aux$Time5[i]<-paste(substr(time.aux$Time[i],1,2),":",substr(time.aux$Time[i],4,4),0,"-",
                             substr(time.aux$Time[i],1,2),":",substr(time.aux$Time[i],4,4),4,sep="")}
    else {time.aux$Time5[i]<-paste(substr(time.aux$Time[i],1,2),":",substr(time.aux$Time[i],4,4),5,"-",
                            substr(time.aux$Time[i],1,2),":",substr(time.aux$Time[i],4,4),9,sep="") }}
all.2010.raw<-merge(all.2010.raw, time.aux, by="Time")
# convert to 5-min data
mymean <- function(X){
 if(sum(is.na(X))>2){
   mymean <- NA
 }else{
   mymean <-mean(X,na.rm=TRUE)
 }
}


all.2010.raw1<-all.2010.raw

all.2010.raw1 <- all.2010.raw1[with(all.2010.raw1,order(GeogAddress,Date1,Time5)),]
all.2010.Time5 <- with(all.2010.raw1,aggregate(list(Flow,Occupancy,Headway,Speed),
                                             list(GeogAddress,Date,Time5), mymean))
names(all.2010.Time5)<-c("GeogAddress","Date","Time5","Flow","Occupancy","Headway","Speed")
nrow(all.2010.Time5)-19*254*288-61*4*288
all.2010.Time5$Date1<-chron(dates = all.2010.Time5$Date, format = "d/m/y",out.format="d/mon/y")
#sort
all.2010.Time5 <- all.2010.Time5[with(all.2010.Time5,order(GeogAddress,Date1,Time5)),]
all.2010.Time5$Flow<-all.2010.Time5$Flow*5
all.2010.Time5$Dweek<-factor(weekdays(all.2010.Time5$Date1),exclude=as.factor(c("Sat","Sun"))) 
all.2010.Time5$Hour<-as.numeric(substr(all.2010.Time5$Time,1,2))
all.2010.t5<-all.2010.Time5
rm(all.2010.Time5)
#save(all.2010.t5,file="/home/terrahome/oaj9/_data/all2010/all.2010.t5.Rdata")
with(all.2010.t5,table(GeogAddress))

load("/home/terrahome/oaj9/_data/all2010/all.2010.t5.Rdata")
table(all.2010.t5$GeogAddress)


#save(all.2010.t5,file="/home/terrahome/oaj9/_data/all2010/all.2010.t5-wk.Rdata")

# Selecting Data site and lagging variables:
#t5.9206B<-all.2010.t30[all.2010.t30$GeogAddress=="M60/9206B" & all.2010.t30$Dweek=="Mon",]
t5.9200B<-all.2010.t5[all.2010.t5$GeogAddress=="M60/9200B",]

# lags
lag1<-function(vbl) c(NA,vbl[-length(vbl)])

t5.9200B$Hway_lag1<-with(t5.9200B,lag1(Headway))
t5.9200B$Occ_lag1<-with(t5.9200B,lag1(Occupancy))
t5.9200B$Flow_lag1<-with(t5.9200B,lag1(Flow))
#save(t5.9200B,file="/home/terrahome/oaj9/_data/all2010/t5.9200B.Rdata")

t5.9200B.May.wk<-t5.9200B
save(t5.9200B.May.wk,file="/home/terrahome/oaj9/_data/all2010/t5.9200B.May.wk.Rdata")



t5.9196B<-all.2010.t5[all.2010.t5$GeogAddress=="M60/9196B",]

# lags
lag1<-function(vbl) c(NA,vbl[-length(vbl)])

t5.9196B$Hway_lag1<-with(t5.9196B,lag1(Headway))
t5.9196B$Occ_lag1<-with(t5.9196B,lag1(Occupancy))
t5.9196B$Flow_lag1<-with(t5.9196B,lag1(Flow))
#save(t5.9196B,file="/home/terrahome/oaj9/_data/all2010/t5.9196B.Rdata")








#---------------------------------------------#
#    aggregating data into 30-min int         #
#---------------------------------------------#
load("/home/terrahome/oaj9/_data/all2010/all.2010.raw.Rdata")
str(all.2010.raw)
time.aux<-data.frame(Time=with(all.2010.raw,as.data.frame(table(Time))$Time))


# Creating interval scale
for (i in 1:1440){
    if (as.numeric(substr(time.aux$Time[i],4,4)) %in%  0:2)
  { time.aux$Time30[i]<-paste(substr(time.aux$Time[i],1,2),":","00","-",
                             substr(time.aux$Time[i],1,2),":",29,sep="")}
  else
   { time.aux$Time30[i]<-paste(substr(time.aux$Time[i],1,2),":",30,"-",
                             substr(time.aux$Time[i],1,2),":",59,sep="")}}

# Checking new interval scale
write.csv(with(time.aux,table(Time, Time30)),
          file="/home/terrahome/oaj9/delme/delme.csv")

all.2010.raw1<-merge(all.2010.raw, time.aux, by="Time")


# Aggregating to 30-min data
mymean <- function(X){
 if(sum(is.na(X))>2){
   mymean <- NA
 }else{
   mymean <-mean(X,na.rm=TRUE)
 }
}
# ordering
all.2010.raw1 <- all.2010.raw1[with(all.2010.raw1,order(GeogAddress,Date1,Time30)),]

all.2010.Time30 <- with(all.2010.raw1,aggregate(list(Flow,Occupancy,Headway,Speed),
                                             list(GeogAddress,Date,Time30), mymean))
names(all.2010.Time30)<-c("GeogAddress","Date","Time30","Flow","Occupancy","Headway","Speed")
nrow(all.2010.Time30)-18*254*48-61*5*48
library(chron)
all.2010.Time30$Date1<-chron(dates = all.2010.Time30$Date, format = "d/m/y",out.format="d/mon/y")
#sort
all.2010.Time30 <- all.2010.Time30[with(all.2010.Time30,order(GeogAddress,Date1,Time30)),]
all.2010.Time30$Flow<-all.2010.Time30$Flow*30
all.2010.Time30$Dweek<-factor(weekdays(all.2010.Time30$Date1),exclude=as.factor(c("Sat","Sun")))
all.2010.Time30$Hour<-as.numeric(substr(all.2010.Time30$Time,1,2))
all.2010.t30<-all.2010.Time30
#save(all.2010.t30,file="/home/terrahome/oaj9/_data/all2010/all.2010.t30.Rdata")
rm(all.2010.Time30)

nrows(all.2010.t30)



#*********************************#
####### 9189B      ################
#*********************************#

# Reading all data:
all.9189B.raw<-read.csv("/home/terrahome/oaj9/_data/all2010/9189B_2010.txt", header = TRUE, as.is = c(1,6,7))
all.9189B.raw$Date1<-chron(dates = all.9189B.raw$Date, format = "d/m/y",out.format="d/mon/y")
#checking if all and dates are in the file (ok 26/01)
write.csv(table(all.9189B.raw$Date1,all.9189B.raw$GeogAddress),
          file="/home/terrahome/oaj9/delme/delme.csv")
#Excluding weekends:
all.9189B.raw<-all.9189B.raw[is.weekend(all.9189B.raw$Date1)==F,]
#Excluding Holidays:
all.9189B.raw<-all.9189B.raw[!all.9189B.raw$Date %in% c("01/01/2010", "02/04/2010","05/04/2010","03/05/2010","31/05/2010",
                           "30/08/2010","25/12/2010","26/12/2010","31/12/2010"),]
all.9189B.raw<-all.9189B.raw[order(all.9189B.raw$GeogAddress,all.9189B.raw$Date,all.9189B.raw$Time),]

all.9189B.raw$Flow <-with(all.9189B.raw,(Category_1_flow+Category_2_flow+Category_3_flow+Category_4_flow)/60)
all.9189B.raw$Occupancy<-rowSums(all.9189B.raw[,c("Lane_1_occupancy","Lane_2_occupancy","Lane_3_occupancy","Lane_4_occupancy")],na.rm=T)/
                     rowSums(is.na(all.9189B.raw[,c("Lane_1_occupancy","Lane_2_occupancy",
                                                "Lane_3_occupancy","Lane_4_occupancy")])==F)
all.9189B.raw$Headway<-rowSums(all.9189B.raw[,c("Lane_1_headway","Lane_2_headway","Lane_3_headway","Lane_4_headway")],na.rm=T)/
                     rowSums(is.na(all.9189B.raw[,c("Lane_1_headway","Lane_2_headway",
                                                "Lane_3_headway","Lane_4_headway")])==F)
#speed:
#sfi=speed x flow @ lane i

all.9189B.raw$sf1<-with(all.9189B.raw,Lane_1_average_speed*Lane_1_flow)
all.9189B.raw$sf2<-with(all.9189B.raw,Lane_2_average_speed*Lane_2_flow)
all.9189B.raw$sf3<-with(all.9189B.raw,Lane_3_average_speed*Lane_3_flow)
all.9189B.raw$sf4<-with(all.9189B.raw,Lane_4_average_speed*Lane_4_flow)

all.9189B.raw$Speed<-rowSums(all.9189B.raw[,c("sf1","sf2","sf3","sf4")],na.rm=T)/
                 rowSums(all.9189B.raw[,c("Lane_1_flow","Lane_2_flow","Lane_3_flow","Lane_4_flow")],na.rm=T)
#table(all.9189B.raw$GeogAddress)
dim(table(all.9189B.raw$Date))

length(all.9189B.raw$Date)-254*1440
#save(all.9189B.raw,file="/home/terrahome/oaj9/_data/all2010/all.9189B.raw.Rdata")

summary(all.9189B.raw$Speed)
summary(all.9189B.raw$Occupancy)
summary(all.9189B.raw$Headway)
summary(all.9189B.raw$Flow)

table(is.na(all.9189B.raw$Speed))/length(all.9189B.raw$Date)
table(is.na(all.9189B.raw$Occupancy))/length(all.9189B.raw$Date)
table(is.na(all.9189B.raw$Headway))/length(all.9189B.raw$Date)
table(is.na(all.9189B.raw$Flow))/length(all.9189B.raw$Date)


#---------------------- ----------------------#
#    aggregating data into 5-min int          #
#---------------------------------------------#
load("/home/terrahome/oaj9/_data/all2010/all.9189B.raw.Rdata")
str(all.9189B.raw)
time.aux<-data.frame(Time=with(all.9189B.raw,as.data.frame(table(Time))$Time))
for (i in 1:1440){
if (as.numeric(substr(time.aux$Time[i],5,5))<5)
  { time.aux$Time5[i]<-paste(substr(time.aux$Time[i],1,2),":",substr(time.aux$Time[i],4,4),0,"-",
                             substr(time.aux$Time[i],1,2),":",substr(time.aux$Time[i],4,4),4,sep="")}
    else {time.aux$Time5[i]<-paste(substr(time.aux$Time[i],1,2),":",substr(time.aux$Time[i],4,4),5,"-",
                            substr(time.aux$Time[i],1,2),":",substr(time.aux$Time[i],4,4),9,sep="") }}
all.9189B.raw1<-merge(all.9189B.raw, time.aux, by="Time")
# convert to 5-min data
mymean <- function(X){
 if(sum(is.na(X))>2){
   mymean <- NA
 }else{
   mymean <-mean(X,na.rm=TRUE)
 }
}
all.9189B.raw1 <- all.9189B.raw1[with(all.9189B.raw1,order(GeogAddress,Date1,Time5)),]
t5.9189B <- with(all.9189B.raw1,aggregate(list(Flow,Occupancy,Headway,Speed),
                                             list(GeogAddress,Date,Time5), mymean))
names(t5.9189B)<-c("GeogAddress","Date","Time5","Flow","Occupancy","Headway","Speed")
nrow(t5.9189B)-254*288
t5.9189B$Date1<-chron(dates = t5.9189B$Date, format = "d/m/y",out.format="d/mon/y")
#sort
t5.9189B <- t5.9189B[with(t5.9189B,order(GeogAddress,Date1,Time5)),]
t5.9189B$Flow<-t5.9189B$Flow*5
t5.9189B$Dweek<-factor(weekdays(t5.9189B$Date1),exclude=as.factor(c("Sat","Sun")))
t5.9189B$Hour<-as.numeric(substr(t5.9189B$Time,1,2))
#save(t5.9189B,file="/home/terrahome/oaj9/_data/all2010/t5.9189B.Rdata")
with(t5.9189B,table(GeogAddress))

# lags
lag1<-function(vbl) c(NA,vbl[-length(vbl)])

t5.9189B$Hway_lag1<-with(t5.9189B,lag1(Headway))
t5.9189B$Occ_lag1<-with(t5.9189B,lag1(Occupancy))
t5.9189B$Flow_lag1<-with(t5.9189B,lag1(Flow))
#save(t5.9189B,file="/home/terrahome/oaj9/_data/all2010/t5.9189B.Rdata")


#*********************************#
####### 1445B and 1438B     ################
#*********************************#

# Reading all data:
add.spl.raw<-read.csv("/home/terrahome/oaj9/_data/all2010/addspl2010.txt", header = TRUE, as.is = c(1,6,7))
add.spl.raw$Date1<-chron(dates = add.spl.raw$Date, format = "d/m/y",out.format="d/mon/y")
#checking if all and dates are in the file (ok 26/01)
write.csv(table(add.spl.raw$Date1,add.spl.raw$GeogAddress),
          file="/home/terrahome/oaj9/delme/delme.csv")
#Excluding weekends:
add.spl.raw<-add.spl.raw[is.weekend(add.spl.raw$Date1)==F,]
#Excluding Holidays:
add.spl.raw<-add.spl.raw[!add.spl.raw$Date %in% c("01/01/2010", "02/04/2010","05/04/2010","03/05/2010","31/05/2010",
                           "30/08/2010","25/12/2010","26/12/2010","31/12/2010"),]
add.spl.raw<-add.spl.raw[order(add.spl.raw$GeogAddress,add.spl.raw$Date,add.spl.raw$Time),]

add.spl.raw$Flow <-with(add.spl.raw,(Category_1_flow+Category_2_flow+Category_3_flow+Category_4_flow)/60)
add.spl.raw$Occupancy<-rowSums(add.spl.raw[,c("Lane_1_occupancy","Lane_2_occupancy","Lane_3_occupancy","Lane_4_occupancy")],na.rm=T)/
                     rowSums(is.na(add.spl.raw[,c("Lane_1_occupancy","Lane_2_occupancy",
                                                "Lane_3_occupancy","Lane_4_occupancy")])==F)
add.spl.raw$Headway<-rowSums(add.spl.raw[,c("Lane_1_headway","Lane_2_headway","Lane_3_headway","Lane_4_headway")],na.rm=T)/
                     rowSums(is.na(add.spl.raw[,c("Lane_1_headway","Lane_2_headway",
                                                "Lane_3_headway","Lane_4_headway")])==F)
#speed:
#sfi=speed x flow @ lane i

add.spl.raw$sf1<-with(add.spl.raw,Lane_1_average_speed*Lane_1_flow)
add.spl.raw$sf2<-with(add.spl.raw,Lane_2_average_speed*Lane_2_flow)
add.spl.raw$sf3<-with(add.spl.raw,Lane_3_average_speed*Lane_3_flow)
add.spl.raw$sf4<-with(add.spl.raw,Lane_4_average_speed*Lane_4_flow)

add.spl.raw$Speed<-rowSums(add.spl.raw[,c("sf1","sf2","sf3","sf4")],na.rm=T)/
                 rowSums(add.spl.raw[,c("Lane_1_flow","Lane_2_flow","Lane_3_flow","Lane_4_flow")],na.rm=T)
#table(add.spl.raw$GeogAddress)
dim(table(add.spl.raw$Date))

length(add.spl.raw$Date)-2*254*1440
#save(add.spl.raw,file="/home/terrahome/oaj9/_data/all2010/add.spl.raw.Rdata")

summary(add.spl.raw$Speed)
summary(add.spl.raw$Occupancy)
summary(add.spl.raw$Headway)
summary(add.spl.raw$Flow)

table(is.na(add.spl.raw$Speed))/length(add.spl.raw$Date)
table(is.na(add.spl.raw$Occupancy))/length(add.spl.raw$Date)
table(is.na(add.spl.raw$Headway))/length(add.spl.raw$Date)
table(is.na(add.spl.raw$Flow))/length(add.spl.raw$Date)


#---------------------- ----------------------#
#    aggregating data into 5-min int          #
#---------------------------------------------#
load("/home/terrahome/oaj9/_data/all2010/add.spl.raw.Rdata")
str(add.spl.raw)
time.aux<-data.frame(Time=with(add.spl.raw,as.data.frame(table(Time))$Time))
for (i in 1:1440){
if (as.numeric(substr(time.aux$Time[i],5,5))<5)
  { time.aux$Time5[i]<-paste(substr(time.aux$Time[i],1,2),":",substr(time.aux$Time[i],4,4),0,"-",
                             substr(time.aux$Time[i],1,2),":",substr(time.aux$Time[i],4,4),4,sep="")}
    else {time.aux$Time5[i]<-paste(substr(time.aux$Time[i],1,2),":",substr(time.aux$Time[i],4,4),5,"-",
                            substr(time.aux$Time[i],1,2),":",substr(time.aux$Time[i],4,4),9,sep="") }}
add.spl.raw1<-merge(add.spl.raw, time.aux, by="Time")
# convert to 5-min data
mymean <- function(X){
 if(sum(is.na(X))>2){
   mymean <- NA
 }else{
   mymean <-mean(X,na.rm=TRUE)
 }
}
add.spl.raw1 <- add.spl.raw1[with(add.spl.raw1,order(GeogAddress,Date1,Time5)),]
t5.add.spl <- with(add.spl.raw1,aggregate(list(Flow,Occupancy,Headway,Speed),
                                             list(GeogAddress,Date,Time5), mymean))
names(t5.add.spl)<-c("GeogAddress","Date","Time5","Flow","Occupancy","Headway","Speed")
nrow(t5.add.spl)-2*254*288
t5.add.spl$Date1<-chron(dates = t5.add.spl$Date, format = "d/m/y",out.format="d/mon/y")
#sort
t5.add.spl <- t5.add.spl[with(t5.add.spl,order(GeogAddress,Date1,Time5)),]
t5.add.spl$Flow<-t5.add.spl$Flow*5
t5.add.spl$Dweek<-factor(weekdays(t5.add.spl$Date1),exclude=as.factor(c("Sat","Sun")))
t5.add.spl$Hour<-as.numeric(substr(t5.add.spl$Time,1,2))
#save(t5.add.spl,file="/home/terrahome/oaj9/_data/all2010/t5.add.spl.Rdata")
with(t5.add.spl,table(GeogAddress))




t5.1445B<-all.2010.t5[all.2010.t5$GeogAddress=="M62/1445B",]

# lags
lag1<-function(vbl) c(NA,vbl[-length(vbl)])

t5.1445B$Hway_lag1<-with(t5.1445B,lag1(Headway))
t5.1445B$Occ_lag1<-with(t5.1445B,lag1(Occupancy))
t5.1445B$Flow_lag1<-with(t5.1445B,lag1(Flow))
#save(t5.1445B,file="/home/terrahome/oaj9/_data/all2010/t5.1445B.Rdata")





t5.1438B<-all.2010.t5[all.2010.t5$GeogAddress=="M62/1438B",]

# lags
lag1<-function(vbl) c(NA,vbl[-length(vbl)])

t5.1438B$Hway_lag1<-with(t5.1438B,lag1(Headway))
t5.1438B$Occ_lag1<-with(t5.1438B,lag1(Occupancy))
t5.1438B$Flow_lag1<-with(t5.1438B,lag1(Flow))
#save(t5.1438B,file="/home/terrahome/oaj9/_data/all2010/t5.1438B.Rdata")





#1437A-1437A-1441A
load(file="/home/terrahome/oaj9/_data/all2010/all.2010.t5.Rdata")

with(all.2010.t5,table(GeogAddress))

#1431A-1437A-1441A
t5.1431A<-all.2010.t5[all.2010.t5$GeogAddress=="M62/1431A",]
# lags
lag1<-function(vbl) c(NA,vbl[-length(vbl)])
t5.1431A$Hway_lag1<-with(t5.1431A,lag1(Headway))
t5.1431A$Occ_lag1<-with(t5.1431A,lag1(Occupancy))
t5.1431A$Flow_lag1<-with(t5.1431A,lag1(Flow))
#save(t5.1431A,file="/home/terrahome/oaj9/_data/all2010/t5.1431A.Rdata")


t5.1437A<-all.2010.t5[all.2010.t5$GeogAddress=="M62/1437A",]
# lags
lag1<-function(vbl) c(NA,vbl[-length(vbl)])
t5.1437A$Hway_lag1<-with(t5.1437A,lag1(Headway))
t5.1437A$Occ_lag1<-with(t5.1437A,lag1(Occupancy))
t5.1437A$Flow_lag1<-with(t5.1437A,lag1(Flow))
#save(t5.1437A,file="/home/terrahome/oaj9/_data/all2010/t5.1437A.Rdata")

t5.1441A<-all.2010.t5[all.2010.t5$GeogAddress=="M62/1441A",]
# lags
lag1<-function(vbl) c(NA,vbl[-length(vbl)])
t5.1441A$Hway_lag1<-with(t5.1441A,lag1(Headway))
t5.1441A$Occ_lag1<-with(t5.1441A,lag1(Occupancy))
t5.1441A$Flow_lag1<-with(t5.1441A,lag1(Flow))
#save(t5.1441A,file="/home/terrahome/oaj9/_data/all2010/t5.1441A.Rdata")




#6013B-6007L-6004L
load(file="/home/terrahome/oaj9/_data/all2010/all.2010.t5.Rdata")

t5.6013B<-all.2010.t5[all.2010.t5$GeogAddress=="M602/6013B",]
# lags
lag1<-function(vbl) c(NA,vbl[-length(vbl)])
t5.6013B$Hway_lag1<-with(t5.6013B,lag1(Headway))
t5.6013B$Occ_lag1<-with(t5.6013B,lag1(Occupancy))
t5.6013B$Flow_lag1<-with(t5.6013B,lag1(Flow))
#save(t5.6013B,file="/home/terrahome/oaj9/_data/all2010/t5.6013B.Rdata")


t5.6007L<-all.2010.t5[all.2010.t5$GeogAddress=="M602/6007L",]
# lags
lag1<-function(vbl) c(NA,vbl[-length(vbl)])
t5.6007L$Hway_lag1<-with(t5.6007L,lag1(Headway))
t5.6007L$Occ_lag1<-with(t5.6007L,lag1(Occupancy))
t5.6007L$Flow_lag1<-with(t5.6007L,lag1(Flow))
#save(t5.6007L,file="/home/terrahome/oaj9/_data/all2010/t5.6007L.Rdata")

t5.6004L<-all.2010.t5[all.2010.t5$GeogAddress=="M602/6004L",]
# lags
lag1<-function(vbl) c(NA,vbl[-length(vbl)])
t5.6004L$Hway_lag1<-with(t5.6004L,lag1(Headway))
t5.6004L$Occ_lag1<-with(t5.6004L,lag1(Occupancy))
t5.6004L$Flow_lag1<-with(t5.6004L,lag1(Flow))
#save(t5.6004L,file="/home/terrahome/oaj9/_data/all2010/t5.6004L.Rdata")

with(t5.6013B,table(months(Date1)))
with(t5.6007L,table(months(Date1)))
with(t5.6004L,table(months(Date1)))

summary(t5.6007L$Flow/t5.6013B$Flow)
summary(t5.6004L$Flow/t5.6013B$Flow)




