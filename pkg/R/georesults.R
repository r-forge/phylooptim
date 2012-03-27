rm(list = ls())

require(gplots)

p <- function(x, digits=4, prefix="", cex.cor, ...)
{
    r <- x$MLE
    txt <- format(c(r, 0.123456789), digits=digits)
    txt <- paste(prefix, txt, sep="")[1:12]
#    if(missing(cex.cor)) cex.cor <- 0.5/strwidth(txt)
#    text(0.5, 0.5, txt, cex = cex.cor )
    return(txt)
}

#Results for aquilegia data
load("/home/michels/Hallowed/repository/phylooptim/pkg/R/geiger/geogeigercts.RData")
geiger_cts <- l
l <- NULL
geiger_cts$mean_lik <- round(c(mean(as.numeric(geiger_cts[[1]][,1])),mean(as.numeric(geiger_cts[[2]][,1])),mean(as.numeric(geiger_cts[[3]][,1])),mean(as.numeric(geiger_cts[[4]][,1])),mean(as.numeric(geiger_cts[[5]][,1])),mean(as.numeric(geiger_cts[[6]][,1])),mean(as.numeric(geiger_cts[[7]][,1])),mean(as.numeric(geiger_cts[[8]][,1])),mean(as.numeric(geiger_cts[[9]][,1])),mean(as.numeric(geiger_cts[[10]][,1])),mean(as.numeric(geiger_cts[[11]][,1])),mean(as.numeric(geiger_cts[[12]][,1]))),4)
names(geiger_cts$mean_lik) <- names(geiger_cts)[1:12]
geiger_cts$MLE <- round(c(max(as.numeric(geiger_cts[[1]][,1])),max(as.numeric(geiger_cts[[2]][,1])),max(as.numeric(geiger_cts[[3]][,1])),max(as.numeric(geiger_cts[[4]][,1])),max(as.numeric(geiger_cts[[5]][,1])),max(as.numeric(geiger_cts[[6]][,1])),max(as.numeric(geiger_cts[[7]][,1])),max(as.numeric(geiger_cts[[8]][,1])),max(as.numeric(geiger_cts[[9]][,1])),max(as.numeric(geiger_cts[[10]][,1])),max(as.numeric(geiger_cts[[11]][,1])),max(as.numeric(geiger_cts[[12]][,1]))),4)
names(geiger_cts$MLE) <- names(geiger_cts)[1:12]

load("/home/michels/Hallowed/repository/phylooptim/pkg/R/geiger/geogeigerdisc.RData")
geiger_disc <- l
l <- NULL
geiger_disc$mean_lik <- round(c(mean(as.numeric(geiger_disc[[1]][,7])),mean(as.numeric(geiger_disc[[2]][,7])),mean(as.numeric(geiger_disc[[3]][,7])),mean(as.numeric(geiger_disc[[4]][,7])),mean(as.numeric(geiger_disc[[5]][,7])),mean(as.numeric(geiger_disc[[6]][,7])),mean(as.numeric(geiger_disc[[7]][,7])),mean(as.numeric(geiger_disc[[8]][,7])),mean(as.numeric(geiger_disc[[9]][,7])),mean(as.numeric(geiger_disc[[10]][,7])),mean(as.numeric(geiger_disc[[11]][,7])),mean(as.numeric(geiger_disc[[12]][,7]))),4)
names(geiger_disc$mean_lik) <- names(geiger_disc)[1:12]
geiger_disc$MLE <- round(c(max(as.numeric(geiger_disc[[1]][,7])),max(as.numeric(geiger_disc[[2]][,7])),max(as.numeric(geiger_disc[[3]][,7])),max(as.numeric(geiger_disc[[4]][,7])),max(as.numeric(geiger_disc[[5]][,7])),max(as.numeric(geiger_disc[[6]][,7])),max(as.numeric(geiger_disc[[7]][,7])),max(as.numeric(geiger_disc[[8]][c(-28,-31),7])),max(as.numeric(geiger_disc[[9]][,7])),max(as.numeric(geiger_disc[[10]][,7])),max(as.numeric(geiger_disc[[11]][,7])),max(as.numeric(geiger_disc[[12]][c(-(15:20),-(22:25)),7]))),4)
names(geiger_disc$MLE) <- names(geiger_disc)[1:12]

load("/home/michels/Hallowed/repository/phylooptim/pkg/R/ace/geoacects.RData")
ace_cts <- l
l <- NULL
ace_cts$mean_lik <- round(c(mean(na.omit(ace_cts[[1]][,5])),mean(na.omit(ace_cts[[2]][,5])),NA,mean(na.omit(ace_cts[[4]][,5])),mean(na.omit(ace_cts[[5]][,5])),mean(na.omit(ace_cts[[6]][,5])),mean(na.omit(ace_cts[[7]][,5])),mean(na.omit(ace_cts[[8]][,5])),mean(na.omit(ace_cts[[9]][,5])),mean(na.omit(ace_cts[[10]][,5])),mean(na.omit(ace_cts[[11]][,5])),mean(na.omit(ace_cts[[12]][,5]))),4)
names(ace_cts$mean_lik) <- names(ace_cts)[1:12]
ace_cts$MLE <- round(c(max(na.omit(ace_cts[[1]][,5])),max(na.omit(ace_cts[[2]][,5])),NA,max(na.omit(ace_cts[[4]][,5])),max(na.omit(ace_cts[[5]][,5])),max(na.omit(ace_cts[[6]][,5])),max(na.omit(ace_cts[[7]][,5])),max(na.omit(ace_cts[[8]][,5])),max(na.omit(ace_cts[[9]][,5])),max(na.omit(ace_cts[[10]][,5])),max(na.omit(ace_cts[[11]][,5])),max(na.omit(ace_cts[[12]][,5]))),4)
names(ace_cts$MLE) <- names(ace_cts)[1:12]

load("/home/michels/Hallowed/repository/phylooptim/pkg/R/geiger/aquigeigerdisc.RData")
ace_disc <- l
ace_disc[[2]] <- ace_disc[[2]][c(-3,-11,-19,-27,-35,-43),]
l <- NULL
ace_disc$mean_lik <- round(c(mean(na.omit(ace_disc[[1]][,5])),mean(na.omit(ace_disc[[2]][,5])),mean(na.omit(ace_disc[[3]][,5])),mean(na.omit(ace_disc[[4]][,5])),mean(na.omit(ace_disc[[5]][,5])),mean(na.omit(ace_disc[[6]][,5])),mean(na.omit(ace_disc[[7]][,5])),mean(na.omit(ace_disc[[8]][,5])),mean(na.omit(ace_disc[[9]][,5])),mean(na.omit(ace_disc[[10]][,5])),mean(na.omit(ace_disc[[11]][,5])),mean(na.omit(ace_disc[[12]][,5]))),4)
names(ace_disc$mean_lik) <- names(ace_disc)[1:12]
ace_disc$MLE <- round(c(max(na.omit(ace_disc[[1]][,5])),max(na.omit(ace_disc[[2]][,5])),max(na.omit(ace_disc[[3]][,5])),max(na.omit(ace_disc[[4]][,5])),max(na.omit(ace_disc[[5]][,5])),max(na.omit(ace_disc[[6]][,5])),max(na.omit(ace_disc[[7]][,5])),max(na.omit(ace_disc[[8]][,5])),max(na.omit(ace_disc[[9]][,5])),max(na.omit(ace_disc[[10]][,5])),max(na.omit(ace_disc[[11]][,5])),max(na.omit(ace_disc[[12]][,5]))),4)
names(ace_disc$MLE) <- names(ace_disc)[1:12]

load("/home/michels/Hallowed/repository/phylooptim/pkg/R/ouch/geoouchcts.RData")
ouch <- l
l <- NULL
ouch$mean_lik <- round(c(mean(ouch[[1]][,5]),mean(ouch[[2]][,5]),mean(ouch[[3]][,5]),mean(ouch[[4]][,5]),mean(ouch[[5]][,5]),mean(ouch[[6]][,5]),mean(ouch[[7]][,5]),mean(ouch[[8]][,5]),mean(ouch[[9]][,5]),mean(ouch[[10]][,5]),mean(ouch[[11]][,5]),mean(ouch[[12]][,5])),4)
names(ouch$mean_lik) <- names(ouch)[1:12]
ouch$MLE <- round(c(max(ouch[[1]][,5]),max(ouch[[2]][,5]),max(ouch[[3]][,5]),max(ouch[[4]][,5]),max(ouch[[5]][,5]),max(ouch[[6]][,5]),max(ouch[[7]][,5]),max(ouch[[8]][,5]),max(ouch[[9]][,5]),max(ouch[[10]][,5]),max(ouch[[11]][,5]),max(ouch[[12]][,5])),4)
names(ouch$MLE) <- names(ouch)[1:12]

#True Liklihood
tglik_cts <- round(max(na.omit(as.numeric(geiger_cts$MLE))),4)
gliktf_cts <- rep(NA,length(geiger_cts$mean_lik))
for (i in c(1:length(gliktf_cts))){if (is.na(geiger_cts$mean_lik[i])){gliktf_cts[i] <- "NA"}else{if (as.numeric(geiger_cts$mean_lik[i])==tglik_cts){gliktf_cts[i] <- "Y"}else{if (tglik_cts-.05 < as.numeric(geiger_cts$mean_lik[i]) && tglik_cts+.05 > as.numeric(geiger_cts$mean_lik[i])){gliktf_cts[i] <- "O"}else{gliktf_cts[i] <- "--"}}}}

for (i in c(1:length(well))){
  prptn <- rep(NA,length(geiger_cts[[i]][,1]))
    for (j in c(1:length(prptn))){
      if (is.na(geiger_cts[[i]][,1][j])){prptn[j] <- 0}else{if (round(as.numeric(geiger_cts[[i]][,1]),4)[j]==tglik_cts){prptn[j] <- 1}else{prptn[j] <- 0}}}
  geiger_cts[[i]] <- cbind(geiger_cts[[i]],prptn)
}

gprptn_cts <- c(round(tglik_cts,3),round(mean(geiger_cts[[1]][,10]),3),round(mean(geiger_cts[[2]][,10]),3),round(mean(geiger_cts[[3]][,10]),3),round(mean(geiger_cts[[4]][,10]),3),round(mean(geiger_cts[[5]][,10]),3),round(mean(geiger_cts[[6]][,10]),3),round(mean(geiger_cts[[7]][,10]),3),round(mean(geiger_cts[[8]][,10]),3),round(mean(geiger_cts[[9]][,10]),3),round(mean(geiger_cts[[10]][,10]),3),round(mean(geiger_cts[[11]][,10]),3),round(mean(geiger_cts[[12]][,10]),3))

tglik_disc <- round(max(na.omit(as.numeric(geiger_disc$MLE))),4)
gliktf_disc <- rep(NA,length(geiger_disc$mean_lik))
for (i in c(1:length(gliktf_disc))){if (is.na(geiger_disc$mean_lik[i])){gliktf_disc[i] <- "NA"}else{if (as.numeric(geiger_disc$mean_lik[i])==tglik_disc){gliktf_disc[i] <- "Y"}else{if (tglik_disc-.05 < as.numeric(geiger_disc$mean_lik[i]) && tglik_disc+.05 > as.numeric(geiger_disc$mean_lik[i])){gliktf_disc[i] <- "O"}else{gliktf_disc[i] <- "--"}}}}

for (i in c(1:length(well))){
  prptn <- rep(NA,length(geiger_disc[[i]][,1]))
    for (j in c(1:length(prptn))){
      if (is.na(geiger_disc[[i]][,7][j])){prptn[j] <- 0}else{if (round(as.numeric(geiger_disc[[i]][,7]),4)[j]==tglik_disc){prptn[j] <- 1}else{prptn[j] <- 0}}}
  geiger_disc[[i]] <- cbind(geiger_disc[[i]],prptn)
}

gprptn_disc <- c(round(tglik_disc,3),round(mean(geiger_disc[[1]][,10]),3),round(mean(geiger_disc[[2]][,10]),3),round(mean(geiger_disc[[3]][,10]),3),round(mean(geiger_disc[[4]][,10]),3),round(mean(geiger_disc[[5]][,10]),3),round(mean(geiger_disc[[6]][,10]),3),round(mean(geiger_disc[[7]][,10]),3),round(mean(geiger_disc[[8]][,10]),3),round(mean(geiger_disc[[9]][,10]),3),round(mean(geiger_disc[[10]][,10]),3),round(mean(geiger_disc[[11]][,10]),3),round(mean(geiger_disc[[12]][,10]),3))
            
talik_cts <- round(max(na.omit(as.numeric(ace_cts$MLE))),4)
aliktf_cts <- rep(NA,length(ace_cts$mean_lik))
for (i in c(1:length(aliktf_cts))){if (is.na(ace_cts$mean_lik[i])){aliktf_cts[i] <- "NA"}else{if (as.numeric(ace_cts$mean_lik[i])==talik_cts){aliktf_cts[i] <- "Y"}else{if (talik_cts-.05 < as.numeric(ace_cts$mean_lik[i]) && talik_cts+.05 > as.numeric(ace_cts$mean_lik[i])){aliktf_cts[i] <- "O"}else{aliktf_cts[i] <- "--"}}}}

for (i in c(1:length(well))){
  prptn <- rep(NA,length(ace_cts[[i]][,1]))
    for (j in c(1:length(prptn))){if (is.na(ace_cts[[i]][,5][j])){prptn[j] <- 0}else{if (round(ace_cts[[i]][,5],4)[j]==talik_cts){prptn[j] <- 1}else{prptn[j] <- 0}}}
  ace_cts[[i]] <- cbind(ace_cts[[i]],prptn)
}
aprptn_cts <- c(round(talik_cts,3),round(mean(ace_cts[[1]][,8]),3),round(mean(ace_cts[[2]][,8]),3),round(mean(ace_cts[[3]][,8]),3),round(mean(ace_cts[[4]][,8]),3),round(mean(ace_cts[[5]][,8]),3),round(mean(ace_cts[[6]][,8]),3),round(mean(ace_cts[[7]][,8]),3),round(mean(ace_cts[[8]][,8]),3),round(mean(ace_cts[[9]][,8]),3),round(mean(ace_cts[[10]][,8]),3),round(mean(ace_cts[[11]][,8]),3),round(mean(ace_cts[[12]][,8]),3))

talik_disc <- round(max(na.omit(as.numeric(ace_disc$MLE))),4)
aliktf_disc <- rep(NA,length(ace_disc$mean_lik))
for (i in c(1:length(aliktf_disc))){if (is.na(ace_disc$mean_lik[i])){aliktf_disc[i] <- "NA"}else{if (as.numeric(ace_disc$mean_lik[i])==talik_disc){aliktf_disc[i] <- "Y"}else{if (talik_disc-.05 < as.numeric(ace_disc$mean_lik[i]) && talik_disc+.05 > as.numeric(ace_disc$mean_lik[i])){aliktf_disc[i] <- "O"}else{aliktf_disc[i] <- "--"}}}}

for (i in c(1:length(well))){
  prptn <- rep(NA,length(ace_disc[[i]][,1]))
    for (j in c(1:length(prptn))){if (is.na(ace_disc[[i]][,5][j])){prptn[j] <- 0}else{if (round(ace_disc[[i]][,5],4)[j]==talik_disc){prptn[j] <- 1}else{prptn[j] <- 0}}}
  ace_disc[[i]] <- cbind(ace_disc[[i]],prptn)
}
aprptn_disc <- c(round(talik_disc,3),round(mean(ace_disc[[1]][,8]),3),round(mean(ace_disc[[2]][,8]),3),round(mean(ace_disc[[3]][,8]),3),round(mean(ace_disc[[4]][,8]),3),round(mean(ace_disc[[5]][,8]),3),round(mean(ace_disc[[6]][,8]),3),round(mean(ace_disc[[7]][,8]),3),round(mean(ace_disc[[8]][,8]),3),round(mean(ace_disc[[9]][,8]),3),round(mean(ace_disc[[10]][,8]),3),round(mean(ace_disc[[11]][,8]),3),round(mean(ace_disc[[12]][,8]),3))
            
tolik <- round(max(na.omit(as.numeric(ouch$MLE))),4)
oliktf <- rep(NA,length(ouch$mean_lik))
for (i in c(1:length(oliktf))){if (is.na(ouch$mean_lik[i])){oliktf[i] <- "NA"}else{if (as.numeric(ouch$mean_lik[i])==tolik){oliktf[i] <- "Y"}else{if (tolik-.05 < as.numeric(ouch$mean_lik[i]) && tolik+.05 > as.numeric(ouch$mean_lik[i])){oliktf[i] <- "O"}else{oliktf[i] <- "--"}}}}

for (i in c(1:length(well))){
  prptn <- rep(NA,length(ouch[[i]][,1]))
    for (j in c(1:length(prptn))){if (is.na(ouch[[i]][,5][j])){prptn[j] <- 0}else{if (round(ouch[[i]][,5],4)[j]==tolik){prptn[j] <- 1}else{prptn[j] <- 0}}}
  ouch[[i]] <- cbind(ouch[[i]],prptn)
}
oprptn <- c(round(tolik,3),round(mean(ouch[[1]][,10]),3),round(mean(ouch[[2]][,10]),3),round(mean(ouch[[3]][,10]),3),round(mean(ouch[[4]][,10]),3),round(mean(ouch[[5]][,10]),3),round(mean(ouch[[6]][,10]),3),round(mean(ouch[[7]][,10],3)),round(mean(ouch[[8]][,10]),3),round(mean(ouch[[9]][,10]),3),round(mean(ouch[[10]][,10]),3),round(mean(ouch[[11]][,10]),3),round(mean(ouch[[12]][,10]),3))

smoke <- matrix(c(gliktf_cts,aliktf_cts,oliktf,gliktf_disc,aliktf_disc),ncol=12,byrow=TRUE)
colnames(smoke) <- names(geiger_cts)[1:12]
rownames(smoke) <- c("geiger_cts","ace_cts","ouch","geiger_disc","ace_disc")
geospiza <- as.table(smoke)
geospiza

smoke1 <- matrix(c(gprptn_cts,aprptn_cts,oprptn,gprptn_disc,aprptn_disc),ncol=13,byrow=TRUE)
colnames(smoke1) <- c("Liklihood",names(ace_cts)[1:12])
rownames(smoke1) <- c("geiger_cts","ace_cts","ouch","geiger_disc","ace_disc")
geoprptn <- as.table(smoke1)
geoprptn

smoke2 <- matrix(c(p(geiger_cts),p(ace_cts),p(ouch),p(geiger_disc),p(ace_disc)),ncol=12,byrow=TRUE)
colnames(smoke2) <- names(geiger_cts)[1:12]
rownames(smoke2) <- c("geiger_cts","ace_cts","ouch","geiger_disc","ace_disc")
geo_lik <- as.table(smoke2)
geo_lik

#Plot geiger_cts
time.table <- vector("list", length(well))
for (i in c(1:length(well))){
  time.table[[i]] <- unlist(as.numeric(geiger_cts[[i]][,9]))*60}
names(time.table) <- well

proportion <- c(mean(geiger_cts[[1]][,10]),mean(geiger_cts[[2]][,10]),mean(geiger_cts[[3]][,10]),mean(geiger_cts[[4]][,10]),mean(geiger_cts[[5]][,10]),mean(geiger_cts[[6]][,10]),mean(geiger_cts[[7]][,10]),mean(geiger_cts[[8]][,10]),mean(geiger_cts[[9]][,10]),mean(geiger_cts[[10]][,10]),mean(geiger_cts[[11]][,10]),mean(geiger_cts[[12]][,10]))

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

png(file="geiger_geo_cts_propvstime.png")
plotCI(time.tab[,1],time.tab[,2],ui=UI,li=LI,xlim=c(min(MIN),max(MAX)),err="x",xlab="Time (Seconds)",ylab="Proportion",main="fitContinuous Geospiza",pch=NA,gap=3)
for (i in c(1:length(well)))
  {
  text(x=time.tab[,1][i],y=time.tab[,2][i],labels=well[i])
}
dev.off()

#Plot geiger_disc
time.table <- vector("list", length(well))
for (i in c(1:length(well))){
  time.table[[i]] <- unlist(as.numeric(geiger_disc[[i]][,9]))*60}
names(time.table) <- well

proportion <- c(mean(geiger_disc[[1]][,10]),mean(geiger_disc[[2]][,10]),mean(geiger_disc[[3]][,10]),mean(geiger_disc[[4]][,10]),mean(geiger_disc[[5]][,10]),mean(geiger_disc[[6]][,10]),mean(geiger_disc[[7]][,10]),mean(geiger_disc[[8]][,10]),mean(geiger_disc[[9]][,10]),mean(geiger_disc[[10]][,10]),mean(geiger_disc[[11]][,10]),mean(geiger_disc[[12]][,10]))

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

png(file="geiger_geo_disc_propvstime.png")
plotCI(time.tab[,1],time.tab[,2],ui=UI,li=LI,xlim=c(min(MIN),max(MAX)),err="x",xlab="Time (Seconds)",ylab="Proportion",main="fitDiscrete Geospiza",pch=NA,gap=3)
for (i in c(1:length(well)))
  {
  text(x=time.tab[,1][i],y=time.tab[,2][i],labels=well[i])
}
dev.off()

#Plot ace_cts
time.table <- vector("list", length(well))
for (i in c(1:length(well))){
  time.table[[i]] <- unlist(ace_cts[[i]][,7])*60}
names(time.table) <- well

proportion <- c(mean(ace_cts[[1]][,8]),mean(ace_cts[[2]][,8]),mean(ace_cts[[3]][,8]),mean(ace_cts[[4]][,8]),mean(ace_cts[[5]][,8]),mean(ace_cts[[6]][,8]),mean(ace_cts[[7]][,8]),mean(ace_cts[[8]][,8]),mean(ace_cts[[9]][,8]),mean(ace_cts[[10]][,8]),mean(ace_cts[[11]][,8]),mean(ace_cts[[12]][,8]))

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

png(file="ace_geo_cts_propvstime.png")
plotCI(time.tab[,1],time.tab[,2],ui=UI,li=LI,xlim=c(min(MIN),max(MAX)),err="x",xlab="Time (Seconds)",ylab="Proportion",main="Ace Cts Geospiza",pch=NA,gap=3)
for (i in c(1:length(well)))
  {
  text(x=time.tab[,1][i],y=time.tab[,2][i],labels=well[i])
}
dev.off()

#Plot ace_disc
time.table <- vector("list", length(well))
for (i in c(1:length(well))){
  time.table[[i]] <- unlist(ace_disc[[i]][,7])*60}
names(time.table) <- well

proportion <- c(mean(ace_disc[[1]][,8]),mean(ace_disc[[2]][,8]),mean(ace_disc[[3]][,8]),mean(ace_disc[[4]][,8]),mean(ace_disc[[5]][,8]),mean(ace_disc[[6]][,8]),mean(ace_disc[[7]][,8]),mean(ace_disc[[8]][,8]),mean(ace_disc[[9]][,8]),mean(ace_disc[[10]][,8]),mean(ace_disc[[11]][,8]),mean(ace_disc[[12]][,8]))

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

png(file="ace_geo_disc_propvstime.png")
plotCI(time.tab[,1],time.tab[,2],ui=UI,li=LI,xlim=c(min(MIN),max(MAX)),err="x",xlab="Time (Seconds)",ylab="Proportion",main="Ace Disc Geospiza",pch=NA,gap=3)
for (i in c(1:length(well)))
  {
  text(x=time.tab[,1][i],y=time.tab[,2][i],labels=well[i])
}
dev.off()

#Plot ouch
time.table <- vector("list", length(well))
for (i in c(1:length(well))){
  time.table[[i]] <- unlist(ouch[[i]][,9])*60}
names(time.table) <- well

proportion <- c(mean(ouch[[1]][,9]),mean(ouch[[2]][,9]),mean(ouch[[3]][,9]),mean(ouch[[4]][,9]),mean(ouch[[5]][,9]),mean(ouch[[6]][,9]),mean(ouch[[7]][,9]),mean(ouch[[8]][,9]),mean(ouch[[9]][,9]),mean(ouch[[10]][,9]),mean(ouch[[11]][,9]),mean(ouch[[12]][,9]))

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

png(file="ouch_geo_propvstime.png")
plotCI(time.tab[,1],time.tab[,2],ui=UI,li=LI,xlim=c(min(MIN),max(MAX)),err="x",xlab="Time (Seconds)",ylab="Proportion",main="Ouch Geospiza",pch=NA,gap=3)
for (i in c(1:length(well)))
  {
  text(x=time.tab[,1][i],y=time.tab[,2][i],labels=well[i])
}
dev.off()
