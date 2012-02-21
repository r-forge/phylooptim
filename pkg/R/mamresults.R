require(gplots)

#Results for mammal data
load("/home/michels/repository/phylooptim/pkg/R/geiger/mamgeigererror.RData")
geiger <- l
l <- NULL
geiger$mean_lik <- round(c(mean(as.numeric(geiger[[1]][,1])),mean(as.numeric(geiger[[2]][,1])),mean(as.numeric(geiger[[3]][,1])),mean(as.numeric(geiger[[4]][,1])),mean(as.numeric(geiger[[5]][,1])),mean(as.numeric(geiger[[6]][,1])),mean(as.numeric(geiger[[7]][,1])),mean(as.numeric(geiger[[8]][,1])),mean(as.numeric(geiger[[9]][,1])),mean(as.numeric(geiger[[10]][,1])),mean(as.numeric(geiger[[11]][,1])),mean(as.numeric(geiger[[12]][,1]))),5)
names(geiger$mean_lik) <- names(geiger)[1:12]
geiger$MLE <- round(c(max(as.numeric(geiger[[1]][,1])),max(as.numeric(geiger[[2]][,1])),max(as.numeric(geiger[[3]][,1])),max(as.numeric(geiger[[4]][,1])),max(as.numeric(geiger[[5]][,1])),max(as.numeric(geiger[[6]][,1])),max(as.numeric(geiger[[7]][,1])),max(as.numeric(geiger[[8]][,1])),max(as.numeric(geiger[[9]][,1])),max(as.numeric(geiger[[10]][,1])),max(as.numeric(geiger[[11]][,1])),max(as.numeric(geiger[[12]][,1]))),5)
names(geiger$MLE) <- names(geiger)[1:12]

load("/home/michels//repository/phylooptim/pkg/R/ace/mamaceerror.RData")
ace <- l
l <- NULL
ace$mean_lik <- round(c(mean(ace[[1]][,5]),mean(ace[[2]][,5]),mean(ace[[3]][,5]),mean(ace[[4]][,5]),mean(ace[[5]][,5]),mean(ace[[6]][,5]),mean(ace[[7]][,5]),mean(ace[[8]][,5]),mean(ace[[9]][,5]),mean(ace[[10]][,5]),mean(ace[[11]][,5]),mean(ace[[12]][,5])),5)
names(ace$mean_lik) <- names(ace)[1:12]
ace$MLE <- c(max(ace[[1]]$L),max(ace[[2]]$L),max(ace[[3]]$L),max(ace[[4]]$L),max(ace[[5]]$L),max(ace[[6]]$L),max(ace[[7]]$L),max(ace[[8]]$L),max(ace[[9]]$L),max(ace[[10]]$L),max(ace[[11]]$L),max(ace[[12]]$L))
names(ace$MLE) <- names(ace)[1:12]

load("/home/michels/repository/phylooptim/pkg/R/ouch/mamoucherror.RData")
ouch <- l
l <- NULL
ouch$mean_lik <- round(c(mean(ouch[[1]][,5]),mean(ouch[[2]][,5]),mean(ouch[[3]][,5]),mean(ouch[[4]][,5]),mean(ouch[[5]][,5]),mean(ouch[[6]][,5]),mean(ouch[[7]][,5]),mean(ouch[[8]][,5]),mean(ouch[[9]][,5]),mean(ouch[[10]][,5]),mean(ouch[[11]][,5]),mean(ouch[[12]][,5])),5)
names(ouch$mean_lik) <- names(ouch)[1:12]
ouch$MLE <- round(c(max(ouch[[1]][,5]),max(ouch[[2]][,5]),max(ouch[[3]][,5]),max(ouch[[4]][,5]),max(ouch[[5]][,5]),max(ouch[[6]][,5]),max(ouch[[7]][,5]),max(ouch[[8]][,5]),max(ouch[[9]][,5]),max(ouch[[10]][,5]),max(ouch[[11]][,5]),max(ouch[[12]][,5])),5)
names(ouch$MLE) <- names(ouch)[1:12]

#True Liklihood
##NOTE: I AM FIXING THIS VALUE BECAUSE RIGHT NOW IT IS WRONG!
#tglik <- round(max(na.omit(as.numeric(geiger$MLE))),4)
tglik <- round(min(na.omit(as.numeric(geiger$MLE))),4)
gliktf <- rep(NA,length(geiger$mean_lik))
for (i in c(1:length(gliktf))){if (is.na(geiger$mean_lik[i])){gliktf[i] <- "NA"}else{if (as.numeric(geiger$mean_lik[i])==tglik){gliktf[i] <- "Y"}else{if (tglik-.05 < as.numeric(geiger$mean_lik[i]) && tglik+.05 > as.numeric(geiger$mean_lik[i])){gliktf[i] <- "O"}else{gliktf[i] <- "--"}}}}

for (i in c(1:length(well))){
  prptn <- rep(NA,length(geiger[[i]][,1]))
    for (j in c(1:length(prptn))){if (is.na(geiger[[i]][,5][j])){prptn[j] <- 0}else{if (round(as.numeric(geiger[[i]][,1]),4)[j]==tglik){prptn[j] <- 1}else{prptn[j] <- 0}}}
  geiger[[i]] <- cbind(geiger[[i]],prptn)
}
gprptn <- c(round(tglik,3),round(mean(geiger[[1]][,9]),3),round(mean(geiger[[2]][,9]),3),round(mean(geiger[[3]][,9]),3),round(mean(geiger[[4]][,9]),3),round(mean(geiger[[5]][,9]),3),round(mean(geiger[[6]][,9]),3),round(mean(geiger[[7]][,9]),3),round(mean(geiger[[8]][,9]),3),round(mean(geiger[[9]][,9]),3),round(mean(geiger[[10]][,9]),3),round(mean(geiger[[11]][,9]),3),round(mean(geiger[[12]][,9]),3))
            
talik <- round(max(na.omit(as.numeric(ace$MLE))),4)
aliktf <- rep(NA,length(ace$mean_lik))
for (i in c(1:length(aliktf))){if (is.na(ace$mean_lik[i])){aliktf[i] <- "NA"}else{if (as.numeric(ace$mean_lik[i])==talik){aliktf[i] <- "Y"}else{if (talik-.05 < as.numeric(ace$mean_lik[i]) && talik+.05 > as.numeric(ace$mean_lik[i])){aliktf[i] <- "O"}else{aliktf[i] <- "--"}}}}

for (i in c(1:length(well))){
  prptn <- rep(NA,length(ace[[i]][,1]))
    for (j in c(1:length(prptn))){if (is.na(ace[[i]][,5][j])){prptn[j] <- 0}else{if (round(ace[[i]][,5],4)[j]==talik){prptn[j] <- 1}else{prptn[j] <- 0}}}
  ace[[i]] <- cbind(ace[[i]],prptn)
}
aprptn <- c(round(talik,3),round(mean(ace[[1]][,7]),3),round(mean(ace[[2]][,7]),3),round(mean(ace[[3]][,7]),3),round(mean(ace[[4]][,7]),3),round(mean(ace[[5]][,7]),3),round(mean(ace[[6]][,7]),3),round(mean(ace[[7]][,7]),3),round(mean(ace[[8]][,7]),3),round(mean(ace[[9]][,7]),3),round(mean(ace[[10]][,7]),3),round(mean(ace[[11]][,7]),3),round(mean(ace[[12]][,7]),3))
            
tolik <- round(max(na.omit(as.numeric(ouch$MLE))),4)
oliktf <- rep(NA,length(ouch$mean_lik))
for (i in c(1:length(oliktf))){if (is.na(ouch$mean_lik[i])){oliktf[i] <- "NA"}else{if (as.numeric(ouch$mean_lik[i])==tolik){oliktf[i] <- "Y"}else{if (tolik-.05 < as.numeric(ouch$mean_lik[i]) && tolik+.05 > as.numeric(ouch$mean_lik[i])){oliktf[i] <- "O"}else{oliktf[i] <- "--"}}}}

for (i in c(1:length(well))){
  prptn <- rep(NA,length(ouch[[i]][,1]))
    for (j in c(1:length(prptn))){if (is.na(ouch[[i]][,5][j])){prptn[j] <- 0}else{if (round(ouch[[i]][,5],4)[j]==tolik){prptn[j] <- 1}else{prptn[j] <- 0}}}
  ouch[[i]] <- cbind(ouch[[i]],prptn)
}
oprptn <- c(round(tolik,3),round(mean(ouch[[1]][,9]),3),round(mean(ouch[[2]][,9]),3),round(mean(ouch[[3]][,9]),3),round(mean(ouch[[4]][,9]),3),round(mean(ouch[[5]][,9]),3),round(mean(ouch[[6]][,9]),3),round(mean(ouch[[7]][,9],3)),round(mean(ouch[[8]][,9]),3),round(mean(ouch[[9]][,9]),3),round(mean(ouch[[10]][,9]),3),round(mean(ouch[[11]][,9]),3),round(mean(ouch[[12]][,9]),3))

smoke <- matrix(c(gliktf,aliktf,oliktf),ncol=12,byrow=TRUE)
colnames(smoke) <- names(ace)[1:12]
rownames(smoke) <- c("geiger","ace","ouch")
mammal <- as.table(smoke)
mammal

smoke1 <- matrix(c(gprptn,aprptn,oprptn),ncol=13,byrow=TRUE)
colnames(smoke1) <- c("Liklihood",names(ace)[1:12])
rownames(smoke1) <- c("geiger","ace","ouch")
mamprptn <- as.table(smoke1)
mamprptn

#Plot geiger
time.table <- vector("list", length(well))
for (i in c(1:length(well))){
  time.table[[i]] <- unlist(as.numeric(geiger[[i]]$time))}
names(time.table) <- well

proportion <- c(mean(na.omit(geiger[[1]][,9])),mean(na.omit(geiger[[2]][,9])),mean(na.omit(geiger[[3]][,9])),mean(na.omit(geiger[[4]][,9])),mean(na.omit(geiger[[5]][,9])),mean(na.omit(geiger[[6]][,9])),mean(na.omit(geiger[[7]][,9])),mean(na.omit(geiger[[8]][,9])),mean(na.omit(geiger[[9]][,9])),mean(na.omit(geiger[[10]][,9])),mean(na.omit(geiger[[11]][,9])),mean(na.omit(geiger[[12]][,9])))

time.tab <- matrix(NA,ncol=2,nrow=length(well))
MIN <- rep(NA,length(well))
MAX <- rep(NA,length(well))
UI <- rep(NA,length(well))
LI <- rep(NA,length(well))
for (i in c(1:length(well))){
  time.tab[i,] <- c(mean(time.table[[i]]),proportion[i])
  MIN[i] <- min(time.table[[i]])
  MAX[i] <- max(time.table[[i]])
  UI[i] <- mean(time.table[[i]])+qnorm(.975)*sd(time.table[[i]])
  LI[i] <- mean(time.table[[i]])+qnorm(.025)*sd(time.table[[i]])
}

png(file="geiger_mam_propvstime.png")
plotCI(time.tab[,1],time.tab[,2],ui=UI,li=LI,xlim=c(min(MIN),max(MAX)),err="x",xlab="Time (Minutes)",ylab="Proportion",main="fitContinuous Mammal",pch=NA,gap=3)
for (i in c(1:length(well)))
  {
  text(x=time.tab[,1][i],y=time.tab[,2][i],labels=well[i])
}
dev.off()

#Plot ace
time.table <- vector("list", length(well))
for (i in c(1:length(well))){
  time.table[[i]] <- unlist(ace[[i]]$time)*60}
names(time.table) <- well

proportion <- c(mean(na.omit(ace[[1]][,7])),mean(na.omit(ace[[2]][,7])),mean(na.omit(ace[[3]][,7])),mean(na.omit(ace[[4]][,7])),mean(na.omit(ace[[5]][,7])),mean(na.omit(ace[[6]][,7])),mean(na.omit(ace[[7]][,7])),mean(na.omit(ace[[8]][,7])),mean(na.omit(ace[[9]][,7])),mean(na.omit(ace[[10]][,7])),mean(na.omit(ace[[11]][,7])),mean(na.omit(ace[[12]][,7])))

time.tab <- matrix(NA,ncol=2,nrow=length(well))
MIN <- rep(NA,length(well))
MAX <- rep(NA,length(well))
UI <- rep(NA,length(well))
LI <- rep(NA,length(well))
for (i in c(1:length(well))){
  time.tab[i,] <- c(mean(time.table[[i]]),proportion[i])
  MIN[i] <- min(time.table[[i]])
  MAX[i] <- max(time.table[[i]])
  UI[i] <- mean(time.table[[i]])+qnorm(.975)*sd(time.table[[i]])
  LI[i] <- mean(time.table[[i]])+qnorm(.025)*sd(time.table[[i]])
}

png(file="ace_mam_propvstime.png")
plotCI(time.tab[,1],time.tab[,2],ui=UI,li=LI,xlim=c(min(MIN),max(MAX)),err="x",xlab="Time (Seconds)",ylab="Proportion",main="Ace Mammal",pch=NA,gap=3)
for (i in c(1:length(well)))
  {
  text(x=time.tab[,1][i],y=time.tab[,2][i],labels=well[i])
}
dev.off()

#Plot ouch
time.table <- vector("list", length(well))
for (i in c(1:length(well))){
  time.table[[i]] <- unlist(ouch[[i]][,8])}
names(time.table) <- well

proportion <- c(mean(na.omit(ouch[[1]][,9])),mean(na.omit(ouch[[2]][,9])),mean(na.omit(ouch[[3]][,9])),mean(na.omit(ouch[[4]][,9])),mean(na.omit(ouch[[5]][,9])),mean(na.omit(ouch[[6]][,9])),mean(na.omit(ouch[[7]][,9])),mean(na.omit(ouch[[8]][,9])),mean(na.omit(ouch[[9]][,9])),mean(na.omit(ouch[[10]][,9])),mean(na.omit(ouch[[11]][,9])),mean(na.omit(ouch[[12]][,9])))

time.tab <- matrix(NA,ncol=2,nrow=length(well))
MIN <- rep(NA,length(well))
MAX <- rep(NA,length(well))
UI <- rep(NA,length(well))
LI <- rep(NA,length(well))
for (i in c(1:length(well))){
  time.tab[i,] <- c(mean(time.table[[i]]),proportion[i])
  MIN[i] <- min(time.table[[i]])
  MAX[i] <- max(time.table[[i]])
  UI[i] <- mean(time.table[[i]])+qnorm(.975)*sd(time.table[[i]])
  LI[i] <- mean(time.table[[i]])+qnorm(.025)*sd(time.table[[i]])
}

png(file="ouch_mam_propvstime.png")
plotCI(time.tab[,1],time.tab[,2],ui=UI,li=LI,xlim=c(min(MIN),max(MAX)),err="x",xlab="Time (Minutes)",ylab="Proportion",main="Ouch Mammal",pch=NA,gap=3)
for (i in c(1:length(well)))
  {
  text(x=time.tab[,1][i],y=time.tab[,2][i],labels=well[i])
}
dev.off()
