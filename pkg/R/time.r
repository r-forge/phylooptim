setwd("/home/michels/Hallowed/repsitory/phylooptim/pkg/R")

require(gplots)

load("geigerctsresults.RData")
full <- list()
full[[1]] <- geiger_geo
full[[2]] <- geiger_aqui
full[[3]] <- geiger_mono
full[[4]] <- geiger_mam
rm(list=ls()[-1*c(which(ls()=="full"))])

load("geigerdiscresults.RData")
full[[5]] <- geiger_geo
full[[6]] <- geiger_aqui
full[[7]] <- geiger_mono
full[[8]] <- geiger_mam
rm(list=ls()[-1*c(which(ls()=="full"))])

load("acectsresults.RData")
full[[9]] <- ace_geo
full[[10]] <- ace_aqui
full[[11]] <- ace_mono
full[[12]] <- ace_mam
rm(list=ls()[-1*c(which(ls()=="full"))])

load("acediscresults.RData")
full[[13]] <- ace_geo
full[[14]] <- ace_aqui
full[[15]] <- ace_mono
full[[16]] <- ace_mam
rm(list=ls()[-1*c(which(ls()=="full"))])

load("ouchctsresults.RData")
full[[17]] <- ouch_geo
full[[18]] <- ouch_aqui
full[[19]] <- ouch_mono
full[[20]] <- ouch_mam
rm(list=ls()[-1*c(which(ls()=="full"))])

well <- c("spg", "Rcgmin", "Rvmmin", "bobyqa","L-BFGS-B","nlminb","ucminf","Nelder-Mead","nlm","CG","BFGS","newuoa")

for (i in c(1:length(full))){for(j in c(1:length(well))){full[[i]][[j]] <- data.frame(full[[i]][[j]])}}

for (j in c(1:length(full))){
  w <- full[[j]]$'Overall MLE'
  for (i in c(1:length(well))){
    for (k in c(1:length(full[[j]][[i]]$time))){
      if (is.na(full[[j]][[i]]$L[k])){full[[j]][[i]]$time[k] <- NA}else{
        if (w <= 0){
          if (as.numeric(full[[j]][[i]]$L[k]) > 2*w & as.numeric(full[[j]][[i]]$L[k]) < -2*w){
            full[[j]][[i]]$time[k] <- full[[j]][[i]]$time[k]}else{
            full[[j]][[i]]$time[k] <- NA}}else{
          if (as.numeric(full[[j]][[i]]$L[k]) > -2*w & as.numeric(full[[j]][[i]]$L[k]) < 2*w){
            full[[j]][[i]]$time[k] <- full[[j]][[i]]$time[k]}else{
            full[[j]][[i]]$time[k] <- NA
          }
        }
      }
    }
  }
}

#By Ace (Ace Continuous)
for (j in c(9:12)){
  for (i in c(1:length(well))){
    a <- full[[j]][[i]][-c(2:50),]
    a[1,6] <- mean(na.omit(full[[j]][[i]]$time))
    full[[j]][[i]] <- a
  }
}

#By Geospiza
t1 <- c(1,5,9,13,17)
time1 <- vector("list", length(well))
prop1 <- vector("list", length(well))
names(time1) <- well
names(prop1) <- well

for (j in c(1:length(t1))){for (i in c(1:length(well))){time1[[i]] <- c(time1[[i]],as.numeric(na.omit(full[[t1[j]]][[i]]$time)));prop1[[i]] <- c(prop1[[i]],as.numeric(full[[t1[j]]][[i]]$prptn))}}

time.tab1 <- matrix(NA,ncol=2,nrow=length(well))
MIN1 <- rep(NA,length(well))
MAX1 <- rep(NA,length(well))
UI1 <- rep(NA,length(well))
LI1 <- rep(NA,length(well))
for (i in c(1:length(well))){
  time.tab1[i,] <- c(mean(time1[[i]]),mean(prop1[[i]]))
  MIN1[i] <- min(time1[[i]])
  MAX1[i] <- max(time1[[i]])
  UI1[i] <- mean(time1[[i]])+qnorm(.975)*sd(time1[[i]])
  LI1[i] <- mean(time1[[i]])+qnorm(.025)*sd(time1[[i]])
}

time.tab1[,1] <- jitter(time.tab1[,1])
time.tab1[,2] <- jitter(time.tab1[,2])
png(file="Geospizatime.png")
#plotCI(time.tab[,1],time.tab[,2],ui=UI,li=LI,xlim=c(min(MIN),max(MAX)),err="x",xlab="Time (Sec)",ylab="Proportion",main="Geospiza",pch=NA,gap=3)
plotCI(time.tab1[,1],time.tab1[,2],ui=UI1,li=LI1,xlim=c(0,550),ylim=c(0,1),err="x",xlab="Time (Min)",ylab="Proportion",main="Geospiza",pch=NA,gap=3)
for (i in c(1:length(well)))
  {
  text(x=time.tab1[,1][i],y=time.tab1[,2][i],labels=well[i])
}
dev.off()

#By Aquilegia
t2 <- c(2,6,10,14,18)
time2 <- vector("list", length(well))
prop2 <- vector("list", length(well))
names(time2) <- well
names(prop2) <- well

for (j in c(1:length(t2))){for (i in c(1:length(well))){time2[[i]] <- c(time2[[i]],as.numeric(na.omit(full[[t2[j]]][[i]]$time)));prop2[[i]] <- c(prop2[[i]],as.numeric(full[[t2[j]]][[i]]$prptn))}}

time.tab2 <- matrix(NA,ncol=2,nrow=length(well))
MIN2 <- rep(NA,length(well))
MAX2 <- rep(NA,length(well))
UI2 <- rep(NA,length(well))
LI2 <- rep(NA,length(well))
for (i in c(1:length(well))){
  time.tab2[i,] <- c(mean(time2[[i]]),mean(prop2[[i]]))
  MIN2[i] <- min(time2[[i]])
  MAX2[i] <- max(time2[[i]])
  UI2[i] <- mean(time2[[i]])+qnorm(.975)*sd(time2[[i]])
  LI2[i] <- mean(time2[[i]])+qnorm(.025)*sd(time2[[i]])
}

time.tab2[,1] <- jitter(time.tab2[,1])
time.tab2[,2] <- jitter(time.tab2[,2])
png(file="Aquilegiatime.png")
plotCI(time.tab2[,1],time.tab2[,2],ui=UI2,li=LI2,xlim=c(0,550),ylim=c(0,1),err="x",xlab="Time (Min)",ylab="Proportion",main="Aquilegia",pch=NA,gap=3)
for (i in c(1:length(well)))
  {
  text(x=time.tab2[,1][i],y=time.tab2[,2][i],labels=well[i])
}
dev.off()

#By Monocot
t3 <- c(3,7,11,15,19)
time3 <- vector("list", length(well))
prop3 <- vector("list", length(well))
names(time3) <- well
names(prop3) <- well

for (j in c(1:length(t3))){for (i in c(1:length(well))){time3[[i]] <- c(time3[[i]],as.numeric(na.omit(full[[t3[j]]][[i]]$time)));prop3[[i]] <- c(prop3[[i]],as.numeric(full[[t3[j]]][[i]]$prptn))}}

time.tab3 <- matrix(NA,ncol=2,nrow=length(well))
MIN3 <- rep(NA,length(well))
MAX3 <- rep(NA,length(well))
UI3 <- rep(NA,length(well))
LI3 <- rep(NA,length(well))
for (i in c(1:length(well))){
  time.tab3[i,] <- c(mean(time3[[i]]),mean(prop3[[i]]))
  MIN3[i] <- min(time3[[i]])
  MAX3[i] <- max(time3[[i]])
  UI3[i] <- mean(time3[[i]])+qnorm(.975)*sd(time3[[i]])
  LI3[i] <- mean(time3[[i]])+qnorm(.025)*sd(time3[[i]])
}

time.tab3[,1] <- jitter(time.tab3[,1])
time.tab3[,2] <- jitter(time.tab3[,2])
png(file="Monocottime.png")
plotCI(time.tab3[,1],time.tab3[,2],ui=UI3,li=LI3,xlim=c(0,550),ylim=c(0,1),err="x",xlab="Time (Min)",ylab="Proportion",main="Monocot",pch=NA,gap=3)
for (i in c(1:length(well)))
  {
  text(x=time.tab3[,1][i],y=time.tab3[,2][i],labels=well[i])
}
dev.off()

#By Mammal
t4 <- c(4,8,12,16,20)
time4 <- vector("list", length(well))
prop4 <- vector("list", length(well))
names(time4) <- well
names(prop4) <- well

for (j in c(1:length(t4))){for (i in c(1:length(well))){time4[[i]] <- c(time4[[i]],as.numeric(na.omit(full[[t4[j]]][[i]]$time)));prop4[[i]] <- c(prop4[[i]],as.numeric(full[[t4[j]]][[i]]$prptn))}}

time.tab4 <- matrix(NA,ncol=2,nrow=length(well))
MIN4 <- rep(NA,length(well))
MAX4 <- rep(NA,length(well))
UI4 <- rep(NA,length(well))
LI4 <- rep(NA,length(well))
for (i in c(1:length(well))){
  time.tab4[i,] <- c(mean(time4[[i]]),mean(prop4[[i]]))
  MIN4[i] <- min(time4[[i]])
  MAX4[i] <- max(time4[[i]])
  UI4[i] <- mean(time4[[i]])+qnorm(.975)*sd(time4[[i]])
  LI4[i] <- mean(time4[[i]])+qnorm(.025)*sd(time4[[i]])
}

time.tab4[,1] <- jitter(time.tab4[,1])
time.tab4[,2] <- jitter(time.tab4[,2])
png(file="Mammaltime.png")
plotCI(time.tab4[,1],time.tab4[,2],ui=UI4,li=LI4,xlim=c(0,550),ylim=c(0,1),err="x",xlab="Time (Min)",ylab="Proportion",main="Mammal",pch=NA,gap=3)
for (i in c(1:length(well)))
  {
  text(x=time.tab4[,1][i],y=time.tab4[,2][i],labels=well[i])
}
dev.off()

#Summary
png(file="datasettime.png")
par(mfrow = c(2,2))
plotCI(time.tab1[,1],time.tab1[,2],ui=UI1,li=LI1,xlim=c(0,1100),ylim=c(0,1),err="x",xlab="Time (Min)",ylab="Proportion",main="Geospiza",pch=NA,gap=3)
for (i in c(1:length(well)))
  {
  text(x=time.tab1[,1][i],y=time.tab1[,2][i],labels=well[i])
}
plotCI(time.tab2[,1],time.tab2[,2],ui=UI2,li=LI2,xlim=c(0,1100),ylim=c(0,1),err="x",xlab="Time (Min)",ylab="Proportion",main="Aquilegia",pch=NA,gap=3)
for (i in c(1:length(well)))
  {
  text(x=time.tab2[,1][i],y=time.tab2[,2][i],labels=well[i])
}
plotCI(time.tab3[,1],time.tab3[,2],ui=UI3,li=LI3,xlim=c(0,1100),ylim=c(0,1),err="x",xlab="Time (Min)",ylab="Proportion",main="Monocot",pch=NA,gap=3)
for (i in c(1:length(well)))
  {
  text(x=time.tab3[,1][i],y=time.tab3[,2][i],labels=well[i])
}
plotCI(time.tab4[,1],time.tab4[,2],ui=UI4,li=LI4,xlim=c(0,1100),ylim=c(0,1),err="x",xlab="Time (Min)",ylab="Proportion",main="Mammal",pch=NA,gap=3)
for (i in c(1:length(well)))
  {
  text(x=time.tab4[,1][i],y=time.tab4[,2][i],labels=well[i])
}
dev.off()

#By geiger (fitContinuous)
time5 <- vector("list", length(well))
prop5 <- vector("list", length(well))
names(time5) <- well
names(prop5) <- well

for (j in c(1:4)){for (i in c(1:length(well))){time5[[i]] <- c(time5[[i]],as.numeric(na.omit(full[[j]][[i]]$time)));prop5[[i]] <- c(prop5[[i]],as.numeric(full[[j]][[i]]$prptn))}}

time.tab5 <- matrix(NA,ncol=2,nrow=length(well))
MIN5 <- rep(NA,length(well))
MAX5 <- rep(NA,length(well))
UI5 <- rep(NA,length(well))
LI5 <- rep(NA,length(well))
for (i in c(1:length(well))){
  time.tab5[i,] <- c(mean(time5[[i]]),mean(prop5[[i]]))
  MIN5[i] <- min(time5[[i]])
  MAX5[i] <- max(time5[[i]])
  UI5[i] <- mean(time5[[i]])+qnorm(.975)*sd(time5[[i]])
  LI5[i] <- mean(time5[[i]])+qnorm(.025)*sd(time5[[i]])
}

time.tab5[,1] <- jitter(time.tab5[,1])
time.tab5[,2] <- jitter(time.tab5[,2])
png(file="fitContinuoustime.png")
plotCI(time.tab5[,1],time.tab5[,2],ui=UI5,li=LI5,xlim=c(0,900),ylim=c(0,1),err="x",xlab="Time (Min)",ylab="Proportion",main="fitContinuous",pch=NA,gap=3)
for (i in c(1:length(well)))
  {
  text(x=time.tab5[,1][i],y=time.tab5[,2][i],labels=well[i])
}
dev.off()

#By geiger (fitDiscrete)
time6 <- vector("list", length(well))
prop6 <- vector("list", length(well))
names(time6) <- well
names(prop6) <- well

for (j in c(5:8)){for (i in c(1:length(well))){time6[[i]] <- c(time6[[i]],as.numeric(na.omit(full[[j]][[i]]$time)));prop6[[i]] <- c(prop6[[i]],as.numeric(full[[j]][[i]]$prptn))}}

time.tab6 <- matrix(NA,ncol=2,nrow=length(well))
MIN6 <- rep(NA,length(well))
MAX6 <- rep(NA,length(well))
UI6 <- rep(NA,length(well))
LI6 <- rep(NA,length(well))
for (i in c(1:length(well))){
  time.tab6[i,] <- c(mean(time6[[i]]),mean(prop6[[i]]))
  MIN6[i] <- min(time6[[i]])
  MAX6[i] <- max(time6[[i]])
  UI6[i] <- mean(time6[[i]])+qnorm(.975)*sd(time6[[i]])
  LI6[i] <- mean(time6[[i]])+qnorm(.025)*sd(time6[[i]])
}

time.tab6[,1] <- jitter(time.tab6[,1])
time.tab6[,2] <- jitter(time.tab6[,2])
png(file="fitDiscretetime.png")
plotCI(time.tab6[,1],time.tab6[,2],ui=UI6,li=LI6,xlim=c(0,900),ylim=c(0,1),err="x",xlab="Time (Min)",ylab="Proportion",main="fitDiscrete",pch=NA,gap=3)
for (i in c(1:length(well)))
  {
  text(x=time.tab6[,1][i],y=time.tab6[,2][i],labels=well[i])
}
dev.off()

#By Ace (Ace Discrete)
time7 <- vector("list", length(well))
prop7 <- vector("list", length(well))
names(time7) <- well
names(prop7) <- well

for (j in c(13:16)){for (i in c(1:length(well))){time7[[i]] <- c(time7[[i]],as.numeric(na.omit(full[[j]][[i]]$time)));prop7[[i]] <- c(prop7[[i]],as.numeric(full[[j]][[i]]$prptn))}}

time.tab7 <- matrix(NA,ncol=2,nrow=length(well))
MIN7 <- rep(NA,length(well))
MAX7 <- rep(NA,length(well))
UI7 <- rep(NA,length(well))
LI7 <- rep(NA,length(well))
for (i in c(1:length(well))){
  time.tab7[i,] <- c(mean(time7[[i]]),mean(prop7[[i]]))
  MIN7[i] <- min(time7[[i]])
  MAX7[i] <- max(time7[[i]])
  UI7[i] <- mean(time7[[i]])+qnorm(.975)*sd(time7[[i]])
  LI7[i] <- mean(time7[[i]])+qnorm(.025)*sd(time7[[i]])
}

time.tab7[,1] <- jitter(time.tab7[,1])
time.tab7[,2] <- jitter(time.tab7[,2])
png(file="AceDisctime.png")
plotCI(time.tab7[,1],time.tab7[,2],ui=UI7,li=LI7,xlim=c(0,900),ylim=c(0,1),err="x",xlab="Time (Min)",ylab="Proportion",main="Ace Discrete",pch=NA,gap=3)
for (i in c(1:length(well)))
  {
  text(x=time.tab7[,1][i],y=time.tab7[,2][i],labels=well[i])
}
dev.off()

#By Ouch
time8 <- vector("list", length(well))
prop8 <- vector("list", length(well))
names(time8) <- well
names(prop8) <- well

for (j in c(17:20)){for (i in c(1:length(well))){time8[[i]] <- c(time8[[i]],as.numeric(na.omit(full[[j]][[i]]$time)));prop8[[i]] <- c(prop8[[i]],as.numeric(full[[j]][[i]]$prptn))}}

time.tab8 <- matrix(NA,ncol=2,nrow=length(well))
MIN8 <- rep(NA,length(well))
MAX8 <- rep(NA,length(well))
UI8 <- rep(NA,length(well))
LI8 <- rep(NA,length(well))
for (i in c(1:length(well))){
  time.tab8[i,] <- c(mean(time8[[i]]),mean(prop8[[i]]))
  MIN8[i] <- min(time8[[i]])
  MAX8[i] <- max(time8[[i]])
  UI8[i] <- mean(time8[[i]])+qnorm(.975)*sd(time8[[i]])
  LI8[i] <- mean(time8[[i]])+qnorm(.025)*sd(time8[[i]])
}

time.tab8[,1] <- jitter(time.tab8[,1])
time.tab8[,2] <- jitter(time.tab8[,2])
png(file="Ouchtime.png")
plotCI(time.tab8[,1],time.tab8[,2],ui=UI8,li=LI8,xlim=c(0,900),ylim=c(0,1),err="x",xlab="Time (Min)",ylab="Proportion",main="Hansen",pch=NA,gap=3)
for (i in c(1:length(well)))
  {
  text(x=time.tab8[,1][i],y=time.tab8[,2][i],labels=well[i])
}
dev.off()

#Everything
t9 <- seq(1,length=22,by=2)
for (j in c(1:length(full))){
  for (i in c(1:length(well))){
    w <- length(full[[j]][[i]][,1])
    if (w==1){next}else{
      if (w > 22){
       ww <- c(1:w)
        tt <- rep(0,w)
        for (k in c(1:length(t9))){tt[t9[k]] <- t9[k]}
        r <- tt-ww
        for (k in c(1:w)){if (r[k]==0){r[k] <- NA}}
        full[[j]][[i]] <- full[[j]][[i]][na.omit(r),]
      }
    }
  }
}

#True time (90%)
time9 <- vector("list", length(well))
prop9 <- vector("list", length(well))
names(time9) <- well
names(prop9) <- well

for (j in c(1:20)){for (i in c(1:length(well))){time9[[i]] <- c(time9[[i]],as.numeric(na.omit(full[[j]][[i]]$time)));prop9[[i]] <- c(prop9[[i]],as.numeric(full[[j]][[i]]$prptn))}}

time.tab9 <- matrix(NA,ncol=5,nrow=length(well))
#MIN9 <- rep(NA,length(well))
#MAX9 <- rep(NA,length(well))

#for (i in c(1:length(well))){
#  time.tab9[i,] <- c(mean(time9[[i]]),mean(prop9[[i]]))
#  MIN9[i] <- min(time9[[i]])
#  MAX9[i] <- max(time9[[i]])
#  UI9[i] <- mean(time9[[i]])+qnorm(.975)*sd(time9[[i]])
#  LI9[i] <- mean(time9[[i]])+qnorm(.025)*sd(time9[[i]])
#}

#LI9[1] <- -85
#LI9[12] <- -82

#time.tab9[,1] <- jitter(time.tab9[,1])
#time.tab9[,2] <- jitter(time.tab9[,2])

for (i in c(1:length(well))){
  time.tab9[i,] <- c(mean(time9[[i]]),mean(prop9[[i]]),quantile(time9[[i]],.05),(quantile(time9[[i]],.95)+quantile(time9[[i]],.05))/2,quantile(time9[[i]],.95))
}

row.names(time.tab9) <- well

time.tab9 <- time.tab9[order(time.tab9[,2]) , ]
t <- time.tab9[order(time.tab9[,2]) , ]

for (i in c(1:(length(well)-1))){if (as.numeric(t[,2][i+1])-as.numeric(t[,2][i]) <=.02){t[,2][i+1] <- t[,2][i]+.02}}

time.tab9 <- cbind(time.tab9,t[,2])
colnames(time.tab9) <- c("time","prop","LI","MI","UI","position")

png(file="Time.png")
plotCI(time.tab9[,1],time.tab9[,6],ui=time.tab9[,5],li=time.tab9[,3],xlim=c(-60,650),ylim=c(0.38,0.82),err="x",xlab="Time (Min)",ylab="Proportion",main=NA,pch=NA,gap=0,axes=FALSE)
points(time.tab9[,1],time.tab9[,6],pch=16,col=2)
for (i in c(1:length(well)))
  {
  text(x=time.tab9[,4][i],y=time.tab9[,6][i]+.008,labels=row.names(time.tab9)[i])
}
x.labels<-c(0,100,200,300,400,500,600)
axis(side = 1, at = x.labels)
y.labels<-round(time.tab9[,6],3)
axis(side = 2, at = y.labels,labels=FALSE,pos=-20)
text(-60, time.tab9[,6], paste(round(time.tab9[,2],3)),cex = .8)
dev.off()

#Skip some x-axis (95%)
for (i in c(1:length(well))){
  time.tab9[i,] <- c(mean(time9[[i]]),mean(prop9[[i]]),quantile(time9[[i]],.025),(quantile(time9[[i]],.975)+quantile(time9[[i]],.025))/2,quantile(time9[[i]],.975))
}

row.names(time.tab9) <- well

time.tab9 <- time.tab9[order(time.tab9[,2]) , ]
t <- time.tab9[order(time.tab9[,2]) , ]

for (i in c(1:(length(well)-1))){if (as.numeric(t[,2][i+1])-as.numeric(t[,2][i]) <=.02){t[,2][i+1] <- t[,2][i]+.02}}

time.tab9 <- cbind(time.tab9,t[,2])
colnames(time.tab9) <- c("time","prop","LI","MI","UI","position")
time.tab9[8,5] <- time.tab9[8,5]-670
time.tab9[8,4] <- (time.tab9[8,5]+time.tab9[8,3])/2

png(file="Time2.png")
plotCI(time.tab9[,1],time.tab9[,6],ui=time.tab9[,5],li=time.tab9[,3],xlim=c(-60,550),ylim=c(0.38,0.82),err="x",xlab="Time (Min)",ylab="Proportion",main=NA,pch=NA,gap=0,axes=FALSE)
points(time.tab9[,1],time.tab9[,6],pch=16,col=2)
points(c(300,300+50/4,300+3*(50/4),350),c(.38,.39,.37,.38),type="l")
for (i in c(1:length(well)))
  {
  text(x=time.tab9[,4][i],y=time.tab9[,6][i]+.008,labels=row.names(time.tab9)[i])
}
x.labels1<-c(0,100,200,300)
x.labels2<-c(350,450,550)
axis(side = 1, at = x.labels1,labels=FALSE,pos=.38)
axis(side = 1, at = x.labels2,labels=FALSE,pos=.38)
mtext(at=c(0,100,200,300,350,450,550), line=-25.3, text=c(0,100,200,300,1000,1100,1200),cex = 1)
y.labels<-round(time.tab9[,6],3)
axis(side = 2, at = y.labels,labels=FALSE,pos=-20)
text(-60, time.tab9[,6], paste(round(time.tab9[,2],3)),cex = .8)
dev.off()

#Log time
time9 <- vector("list", length(well))
prop9 <- vector("list", length(well))
names(time9) <- well
names(prop9) <- well

for (j in c(1:20)){for (i in c(1:length(well))){time9[[i]] <- c(time9[[i]],log(as.numeric(na.omit(full[[j]][[i]]$time))));prop9[[i]] <- c(prop9[[i]],as.numeric(full[[j]][[i]]$prptn))}}

time.tab9 <- matrix(NA,ncol=5,nrow=length(well))

for (i in c(1:length(well))){
  time.tab9[i,] <- c(mean(time9[[i]]),mean(prop9[[i]]),quantile(time9[[i]],.025),(quantile(time9[[i]],.975)+quantile(time9[[i]],.025))/2,quantile(time9[[i]],.975))
}

row.names(time.tab9) <- well

time.tab9 <- time.tab9[order(time.tab9[,2]) , ]
t <- time.tab9[order(time.tab9[,2]) , ]

for (i in c(1:(length(well)-1))){if (as.numeric(t[,2][i+1])-as.numeric(t[,2][i]) <=.02){t[,2][i+1] <- t[,2][i]+.02}}

time.tab9 <- cbind(time.tab9,t[,2])
colnames(time.tab9) <- c("time","prop","LI","MI","UI","position")

png(file="logTime.png")
plotCI(time.tab9[,1],time.tab9[,6],ui=time.tab9[,5],li=time.tab9[,3],xlim=c(-8,8),ylim=c(0.38,0.82),err="x",xlab="Time (Min)",ylab="Proportion",main=NA,pch=NA,gap=4,axes=FALSE)
for (i in c(1:length(well)))
  {
  text(x=time.tab9[,1][i],y=time.tab9[,6][i],labels=row.names(time.tab9)[i])
}
x.ticks <- c(-7,log(10), log(60), log(180), log(600))
x.labels <- c(0, 10, 60, 180, 600)
axis(side = 1, at = x.ticks,labels=x.labels)
#text(x=x.ticks[1], y=.34, label=x.labels[1])
mtext(text=x.labels[4], at = 5.2, line = -26.1)
y.labels<-round(time.tab9[,6],3)
axis(side = 2, at = y.labels,labels=FALSE,pos=-7)
text(-8, time.tab9[,6], paste(round(time.tab9[,2],3)),cex = .8)
dev.off()

#Summary
png(file="functiontime.png")
par(mfrow = c(2,2))
plotCI(time.tab5[,1],time.tab5[,2],ui=UI5,li=LI5,xlim=c(0,900),ylim=c(0,1),err="x",xlab="Time (Min)",ylab="Proportion",main="fitContinuous",pch=NA,gap=3)
for (i in c(1:length(well)))
  {
  text(x=time.tab5[,1][i],y=time.tab5[,2][i],labels=well[i])
}
plotCI(time.tab6[,1],time.tab6[,2],ui=UI6,li=LI6,xlim=c(0,900),ylim=c(0,1),err="x",xlab="Time (Min)",ylab="Proportion",main="fitDiscrete",pch=NA,gap=3)
for (i in c(1:length(well)))
  {
  text(x=time.tab6[,1][i],y=time.tab6[,2][i],labels=well[i])
}
plotCI(time.tab7[,1],time.tab7[,2],ui=UI7,li=LI7,xlim=c(0,900),ylim=c(0,1),err="x",xlab="Time (Min)",ylab="Proportion",main="Ace Discrete",pch=NA,gap=3)
for (i in c(1:length(well)))
  {
  text(x=time.tab7[,1][i],y=time.tab7[,2][i],labels=well[i])
}
plotCI(time.tab8[,1],time.tab8[,2],ui=UI8,li=LI8,xlim=c(0,900),ylim=c(0,1),err="x",xlab="Time (Min)",ylab="Proportion",main="Hansen",pch=NA,gap=3)
for (i in c(1:length(well)))
  {
  text(x=time.tab8[,1][i],y=time.tab8[,2][i],labels=well[i])
}
dev.off()
