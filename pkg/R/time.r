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
t <- c(1,5,9,13,17)
time <- vector("list", length(well))
prop <- vector("list", length(well))
names(time) <- well
names(prop) <- well

for (j in c(1:length(t))){for (i in c(1:length(well))){time[[i]] <- c(time[[i]],as.numeric(na.omit(full[[t[j]]][[i]]$time))*60);prop[[i]] <- c(prop[[i]],as.numeric(full[[t[j]]][[i]]$prptn))}}

time.tab <- matrix(NA,ncol=2,nrow=length(well))
MIN <- rep(NA,length(well))
MAX <- rep(NA,length(well))
UI <- rep(NA,length(well))
LI <- rep(NA,length(well))
for (i in c(1:length(well))){
  time.tab[i,] <- c(mean(time[[i]]),mean(prop[[i]]))
  MIN[i] <- min(time[[i]])
  MAX[i] <- max(time[[i]])
  UI[i] <- mean(time[[i]])+qnorm(.975)*sd(time[[i]])
  LI[i] <- mean(time[[i]])+qnorm(.025)*sd(time[[i]])
}

time.tab[,1] <- jitter(time.tab[,1])
time.tab[,2] <- jitter(time.tab[,2])
png(file="Geospizatime.png")
#plotCI(time.tab[,1],time.tab[,2],ui=UI,li=LI,xlim=c(min(MIN),max(MAX)),err="x",xlab="Time (Sec)",ylab="Proportion",main="Geospiza",pch=NA,gap=3)
plotCI(time.tab[,1],time.tab[,2],ui=UI,li=LI,xlim=c(0,35),err="x",xlab="Time (Sec)",ylab="Proportion",main="Geospiza",pch=NA,gap=3)
for (i in c(1:length(well)))
  {
  text(x=time.tab[,1][i],y=time.tab[,2][i],labels=well[i])
}
dev.off()

#By Aquilegia
t <- c(2,6,10,14,18)
time <- vector("list", length(well))
prop <- vector("list", length(well))
names(time) <- well
names(prop) <- well

for (j in c(1:length(t))){for (i in c(1:length(well))){time[[i]] <- c(time[[i]],as.numeric(na.omit(full[[t[j]]][[i]]$time))*60);prop[[i]] <- c(prop[[i]],as.numeric(full[[t[j]]][[i]]$prptn))}}

time.tab <- matrix(NA,ncol=2,nrow=length(well))
MIN <- rep(NA,length(well))
MAX <- rep(NA,length(well))
UI <- rep(NA,length(well))
LI <- rep(NA,length(well))
for (i in c(1:length(well))){
  time.tab[i,] <- c(mean(time[[i]]),mean(prop[[i]]))
  MIN[i] <- min(time[[i]])
  MAX[i] <- max(time[[i]])
  UI[i] <- mean(time[[i]])+qnorm(.975)*sd(time[[i]])
  LI[i] <- mean(time[[i]])+qnorm(.025)*sd(time[[i]])
}

time.tab[,1] <- jitter(time.tab[,1])
time.tab[,2] <- jitter(time.tab[,2])
png(file="Aquilegiatime.png")
plotCI(time.tab[,1],time.tab[,2],ui=UI,li=LI,xlim=c(0,40),err="x",xlab="Time (Sec)",ylab="Proportion",main="Aquilegia",pch=NA,gap=3)
for (i in c(1:length(well)))
  {
  text(x=time.tab[,1][i],y=time.tab[,2][i],labels=well[i])
}
dev.off()

#By Monocot
t <- c(3,7,11,15,19)
time <- vector("list", length(well))
prop <- vector("list", length(well))
names(time) <- well
names(prop) <- well

for (j in c(1:length(t))){for (i in c(1:length(well))){time[[i]] <- c(time[[i]],as.numeric(na.omit(full[[t[j]]][[i]]$time)));prop[[i]] <- c(prop[[i]],as.numeric(full[[t[j]]][[i]]$prptn))}}

time.tab <- matrix(NA,ncol=2,nrow=length(well))
MIN <- rep(NA,length(well))
MAX <- rep(NA,length(well))
UI <- rep(NA,length(well))
LI <- rep(NA,length(well))
for (i in c(1:length(well))){
  time.tab[i,] <- c(mean(time[[i]]),mean(prop[[i]]))
  MIN[i] <- min(time[[i]])
  MAX[i] <- max(time[[i]])
  UI[i] <- mean(time[[i]])+qnorm(.975)*sd(time[[i]])
  LI[i] <- mean(time[[i]])+qnorm(.025)*sd(time[[i]])
}

time.tab[,1] <- jitter(time.tab[,1])
time.tab[,2] <- jitter(time.tab[,2])
png(file="Monocottime.png")
plotCI(time.tab[,1],time.tab[,2],ui=UI,li=LI,xlim=c(0,175),err="x",xlab="Time (Min)",ylab="Proportion",main="Monocot",pch=NA,gap=3)
for (i in c(1:length(well)))
  {
  text(x=time.tab[,1][i],y=time.tab[,2][i],labels=well[i])
}
dev.off()

#By Mammal
t <- c(4,8,12,16,20)
time <- vector("list", length(well))
prop <- vector("list", length(well))
names(time) <- well
names(prop) <- well

for (j in c(1:length(t))){for (i in c(1:length(well))){time[[i]] <- c(time[[i]],as.numeric(na.omit(full[[t[j]]][[i]]$time)));prop[[i]] <- c(prop[[i]],as.numeric(full[[t[j]]][[i]]$prptn))}}

time.tab <- matrix(NA,ncol=2,nrow=length(well))
MIN <- rep(NA,length(well))
MAX <- rep(NA,length(well))
UI <- rep(NA,length(well))
LI <- rep(NA,length(well))
for (i in c(1:length(well))){
  time.tab[i,] <- c(mean(time[[i]]),mean(prop[[i]]))
  MIN[i] <- min(time[[i]])
  MAX[i] <- max(time[[i]])
  UI[i] <- mean(time[[i]])+qnorm(.975)*sd(time[[i]])
  LI[i] <- mean(time[[i]])+qnorm(.025)*sd(time[[i]])
}

time.tab[,1] <- jitter(time.tab[,1])
time.tab[,2] <- jitter(time.tab[,2])
png(file="Mammaltime.png")
plotCI(time.tab[,1],time.tab[,2],ui=UI,li=LI,xlim=c(0,1200),err="x",xlab="Time (Min)",ylab="Proportion",main="Mammal",pch=NA,gap=3)
for (i in c(1:length(well)))
  {
  text(x=time.tab[,1][i],y=time.tab[,2][i],labels=well[i])
}
dev.off()

#By geiger (fitContinuous)
time <- vector("list", length(well))
prop <- vector("list", length(well))
names(time) <- well
names(prop) <- well

for (j in c(1:4)){for (i in c(1:length(well))){time[[i]] <- c(time[[i]],as.numeric(na.omit(full[[j]][[i]]$time)));prop[[i]] <- c(prop[[i]],as.numeric(full[[j]][[i]]$prptn))}}

time.tab <- matrix(NA,ncol=2,nrow=length(well))
MIN <- rep(NA,length(well))
MAX <- rep(NA,length(well))
UI <- rep(NA,length(well))
LI <- rep(NA,length(well))
for (i in c(1:length(well))){
  time.tab[i,] <- c(mean(time[[i]]),mean(prop[[i]]))
  MIN[i] <- min(time[[i]])
  MAX[i] <- max(time[[i]])
  UI[i] <- mean(time[[i]])+qnorm(.975)*sd(time[[i]])
  LI[i] <- mean(time[[i]])+qnorm(.025)*sd(time[[i]])
}

time.tab[,1] <- jitter(time.tab[,1])
time.tab[,2] <- jitter(time.tab[,2])
png(file="fitContinuoustime.png")
plotCI(time.tab[,1],time.tab[,2],ui=UI,li=LI,xlim=c(0,900),err="x",xlab="Time (Min)",ylab="Proportion",main="fitContinuous",pch=NA,gap=3)
for (i in c(1:length(well)))
  {
  text(x=time.tab[,1][i],y=time.tab[,2][i],labels=well[i])
}
dev.off()

#By geiger (fitDiscrete)
time <- vector("list", length(well))
prop <- vector("list", length(well))
names(time) <- well
names(prop) <- well

for (j in c(5:8)){for (i in c(1:length(well))){time[[i]] <- c(time[[i]],as.numeric(na.omit(full[[j]][[i]]$time)));prop[[i]] <- c(prop[[i]],as.numeric(full[[j]][[i]]$prptn))}}

time.tab <- matrix(NA,ncol=2,nrow=length(well))
MIN <- rep(NA,length(well))
MAX <- rep(NA,length(well))
UI <- rep(NA,length(well))
LI <- rep(NA,length(well))
for (i in c(1:length(well))){
  time.tab[i,] <- c(mean(time[[i]]),mean(prop[[i]]))
  MIN[i] <- min(time[[i]])
  MAX[i] <- max(time[[i]])
  UI[i] <- mean(time[[i]])+qnorm(.975)*sd(time[[i]])
  LI[i] <- mean(time[[i]])+qnorm(.025)*sd(time[[i]])
}

time.tab[,1] <- jitter(time.tab[,1])
time.tab[,2] <- jitter(time.tab[,2])
png(file="fitDiscretetime.png")
plotCI(time.tab[,1],time.tab[,2],ui=UI,li=LI,xlim=c(0,20),err="x",xlab="Time (Min)",ylab="Proportion",main="fitDiscrete",pch=NA,gap=3)
for (i in c(1:length(well)))
  {
  text(x=time.tab[,1][i],y=time.tab[,2][i],labels=well[i])
}
dev.off()

#By Ace (Ace Discrete)
time <- vector("list", length(well))
prop <- vector("list", length(well))
names(time) <- well
names(prop) <- well

for (j in c(13:16)){for (i in c(1:length(well))){time[[i]] <- c(time[[i]],as.numeric(na.omit(full[[j]][[i]]$time))*60);prop[[i]] <- c(prop[[i]],as.numeric(full[[j]][[i]]$prptn))}}

time.tab <- matrix(NA,ncol=2,nrow=length(well))
MIN <- rep(NA,length(well))
MAX <- rep(NA,length(well))
UI <- rep(NA,length(well))
LI <- rep(NA,length(well))
for (i in c(1:length(well))){
  time.tab[i,] <- c(mean(time[[i]]),mean(prop[[i]]))
  MIN[i] <- min(time[[i]])
  MAX[i] <- max(time[[i]])
  UI[i] <- mean(time[[i]])+qnorm(.975)*sd(time[[i]])
  LI[i] <- mean(time[[i]])+qnorm(.025)*sd(time[[i]])
}

time.tab[,1] <- jitter(time.tab[,1])
time.tab[,2] <- jitter(time.tab[,2])
png(file="AceDisctime.png")
plotCI(time.tab[,1],time.tab[,2],ui=UI,li=LI,xlim=c(0,50),err="x",xlab="Time (Sec)",ylab="Proportion",main="Ace Discrete",pch=NA,gap=3)
for (i in c(1:length(well)))
  {
  text(x=time.tab[,1][i],y=time.tab[,2][i],labels=well[i])
}
dev.off()

#By Ouch
time <- vector("list", length(well))
prop <- vector("list", length(well))
names(time) <- well
names(prop) <- well

for (j in c(17:20)){for (i in c(1:length(well))){time[[i]] <- c(time[[i]],as.numeric(na.omit(full[[j]][[i]]$time)));prop[[i]] <- c(prop[[i]],as.numeric(full[[j]][[i]]$prptn))}}

time.tab <- matrix(NA,ncol=2,nrow=length(well))
MIN <- rep(NA,length(well))
MAX <- rep(NA,length(well))
UI <- rep(NA,length(well))
LI <- rep(NA,length(well))
for (i in c(1:length(well))){
  time.tab[i,] <- c(mean(time[[i]]),mean(prop[[i]]))
  MIN[i] <- min(time[[i]])
  MAX[i] <- max(time[[i]])
  UI[i] <- mean(time[[i]])+qnorm(.975)*sd(time[[i]])
  LI[i] <- mean(time[[i]])+qnorm(.025)*sd(time[[i]])
}

time.tab[,1] <- jitter(time.tab[,1])
time.tab[,2] <- jitter(time.tab[,2])
png(file="Ouchtime.png")
plotCI(time.tab[,1],time.tab[,2],ui=UI,li=LI,xlim=c(0,50),err="x",xlab="Time (Min)",ylab="Proportion",main="Hansen",pch=NA,gap=3)
for (i in c(1:length(well)))
  {
  text(x=time.tab[,1][i],y=time.tab[,2][i],labels=well[i])
}
dev.off()

#Everything
t <- seq(1,length=22,by=2)
for (j in c(1:length(full))){
  for (i in c(1:length(well))){
    w <- length(full[[j]][[i]][,1])
    if (w==1){next}else{
      if (w > 22){
       ww <- c(1:w)
        tt <- rep(0,w)
        for (k in c(1:length(t))){tt[t[k]] <- t[k]}
        r <- tt-ww
        for (k in c(1:w)){if (r[k]==0){r[k] <- NA}}
        full[[j]][[i]] <- full[[j]][[i]][na.omit(r),]
      }
    }
  }
}

time <- vector("list", length(well))
prop <- vector("list", length(well))
names(time) <- well
names(prop) <- well

for (j in c(1:20)){for (i in c(1:length(well))){time[[i]] <- c(time[[i]],as.numeric(na.omit(full[[j]][[i]]$time)));prop[[i]] <- c(prop[[i]],as.numeric(full[[j]][[i]]$prptn))}}

time.tab <- matrix(NA,ncol=2,nrow=length(well))
MIN <- rep(NA,length(well))
MAX <- rep(NA,length(well))
UI <- rep(NA,length(well))
LI <- rep(NA,length(well))
for (i in c(1:length(well))){
  time.tab[i,] <- c(mean(time[[i]]),mean(prop[[i]]))
  MIN[i] <- min(time[[i]])
  MAX[i] <- max(time[[i]])
  UI[i] <- mean(time[[i]])+qnorm(.975)*sd(time[[i]])
  LI[i] <- mean(time[[i]])+qnorm(.025)*sd(time[[i]])
}

time.tab[,1] <- jitter(time.tab[,1])
time.tab[,2] <- jitter(time.tab[,2])
png(file="Time.png")
plotCI(time.tab[,1],time.tab[,2],ui=UI,li=LI,xlim=c(0,550),err="x",xlab="Time (Min)",ylab="Proportion",main="Everything",pch=NA,gap=3)
for (i in c(1:length(well)))
  {
  text(x=time.tab[,1][i],y=time.tab[,2][i],labels=well[i])
}
dev.off()

