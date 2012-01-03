#Results for geospiza
load("/home/michels/Hallowed/Dropbox/repository/phylooptim/pkg/R/geiger/geigererror.RData")
geiger <- l2
l2 <- NULL

l <- NULL
geiger$mean_lik <- round(c(mean(geiger[[1]][,2]),mean(geiger[[2]][,2]),mean(geiger[[3]][,2]),mean(geiger[[4]][,2]),mean(geiger[[5]][,2]),mean(geiger[[6]][,2]),mean(geiger[[7]][,2]),mean(geiger[[8]][,2]),mean(geiger[[9]][,2]),mean(geiger[[10]][,2]),mean(geiger[[11]][,2]),mean(geiger[[12]][,2])),5)
names(geiger$mean_lik) <- names(geiger)[1:12]
                 
load("/home/michels/Hallowed/Dropbox/repository/phylooptim/pkg/R/ace/aceerror.RData")
ace <- l
l <- NULL
ace$mean_lik <- round(c(mean(ace[[1]][,5]),mean(ace[[2]][,5]),mean(ace[[3]][,5]),mean(ace[[4]][,5]),mean(ace[[5]][,5]),mean(ace[[6]][,5]),mean(ace[[7]][,5]),mean(ace[[8]][,5]),mean(ace[[9]][,5]),mean(ace[[10]][,5]),mean(ace[[11]][,5]),mean(ace[[12]][,5])),5)
names(ace$mean_lik) <- names(ace)[1:12]

load("/home/michels/Hallowed/Dropbox/repository/phylooptim/pkg/R/ouch/oucherror.RData")
ouch <- l
l <- NULL
ouch$mean_lik <- round(c(mean(ouch[[1]][,5]),mean(ouch[[2]][,5]),mean(ouch[[3]][,5]),mean(ouch[[4]][,5]),mean(ouch[[5]][,5]),mean(ouch[[6]][,5]),mean(ouch[[7]][,5]),mean(ouch[[8]][,5]),mean(ouch[[9]][,5]),mean(ouch[[10]][,5]),mean(ouch[[11]][,5]),mean(ouch[[12]][,5])),5)
names(ouch$mean_lik) <- names(ouch)[1:12]

#Find the mode
f <- function(x){
xt <- table(x)
t <- as.numeric(names(xt[xt == max(xt)]))  
return(t)}

#True Liklihood
tglik <- max(f(as.numeric(geiger$mean_lik)))
gliktf <- rep(NA,length(geiger$mean_lik))
for (i in c(1:length(gliktf))){if (as.numeric(geiger$mean_lik[i])==tglik){gliktf[i] <- "Y"}else{if (tglik-.05 < as.numeric(geiger$mean_lik[i]) && tglik+.05 > as.numeric(geiger$mean_lik[i])){gliktf[i] <- "O"}else{gliktf[i] <- "--"}}}

talik <- max(f(as.numeric(ace$mean_lik)))
aliktf <- rep(NA,length(ace$mean_lik))
for (i in c(1:length(aliktf))){if (as.numeric(ace$mean_lik[i])==talik){aliktf[i] <- "Y"}else{if (talik-.05 < as.numeric(ace$mean_lik[i]) && talik+.05 > as.numeric(ace$mean_lik[i])){aliktf[i] <- "O"}else{aliktf[i] <- "--"}}}

tolik <- max(f(as.numeric(ouch$mean_lik)))
oliktf <- rep(NA,length(ouch$mean_lik))
for (i in c(1:length(oliktf))){if (as.numeric(ouch$mean_lik[i])==tolik){oliktf[i] <- "Y"}else{if (tolik-.05 < as.numeric(ouch$mean_lik[i]) && tolik+.05 > as.numeric(ouch$mean_lik[i])){oliktf[i] <- "O"}else{oliktf[i] <- "--"}}}

smoke <- matrix(c(gliktf,aliktf,oliktf),ncol=12,byrow=TRUE)
colnames(smoke) <- names(geiger)[1:12]
rownames(smoke) <- c("geiger","ace","ouch")
geospiza <- as.table(smoke)
geospiza

#Plot geiger
time.table <- vector("list", length(well))
for (i in c(1:length(well))){
  time.table[[i]] <- unlist(geiger[[i]]$T)}
names(time.table) <- well

for (i in c(1:length(well))){
  prptn <- rep(NA,length(geiger[[i]]$T))
    for (j in c(1:length(prptn))){if (round(geiger[[i]][,2],5)[j]==tglik){prptn[j] <- 1}else{prptn[j] <- 0}}
  geiger[[i]] <- cbind(geiger[[i]],prptn)
}

proportion <- c(mean(geiger[[1]][,10]),mean(geiger[[2]][,10]),mean(geiger[[3]][,10]),mean(geiger[[4]][,10]),mean(geiger[[5]][,10]),mean(geiger[[6]][,10]),mean(geiger[[7]][,10]),mean(geiger[[8]][,10]),mean(geiger[[9]][,10]),mean(geiger[[10]][,10]),mean(geiger[[11]][,10]),mean(geiger[[12]][,10]))

time.tab <- matrix(NA,ncol=2,nrow=length(well))
for (i in c(1:length(well))){time.tab[i,] <- c(mean(time.table[[i]]),proportion[i])}

png(file="geiger_small_propvstime.png")
plot(time.tab[,1],time.tab[,2],pch=NA,xlab="Time",ylab="Proportion",main="fitContinuous")
for (i in c(1:length(well))){
  text(x=time.tab[,1][i],y=time.tab[,2][i],labels=well[i])
}
dev.off()

#Plot ace
time.table <- vector("list", length(well))
for (i in c(1:length(well))){
  time.table[[i]] <- unlist(ace[[i]]$T)}
names(time.table) <- well

for (i in c(1:length(well))){
  prptn <- rep(NA,length(ace[[i]]$T))
    for (j in c(1:length(prptn))){if (round(ace[[i]][,5],5)[j]==talik){prptn[j] <- 1}else{prptn[j] <- 0}}
  ace[[i]] <- cbind(ace[[i]],prptn)
}

proportion <- c(mean(ace[[1]][,7]),mean(ace[[2]][,7]),mean(ace[[3]][,7]),mean(ace[[4]][,7]),mean(ace[[5]][,7]),mean(ace[[6]][,7]),mean(ace[[7]][,7]),mean(ace[[8]][,7]),mean(ace[[9]][,7]),mean(ace[[10]][,7]),mean(ace[[11]][,7]),mean(ace[[12]][,7]))

time.tab <- matrix(NA,ncol=2,nrow=length(well))
for (i in c(1:length(well))){time.tab[i,] <- c(mean(time.table[[i]]),proportion[i])}

png(file="ace_small_propvstime.png")
plot(time.tab[,1],time.tab[,2],pch=NA,xlab="Time",ylab="Proportion",main="Ace")
for (i in c(1:length(well))){
  text(x=time.tab[,1][i],y=time.tab[,2][i],labels=well[i])
}
dev.off()

#Plot ouch
time.table <- vector("list", length(well))
for (i in c(1:length(well))){
  time.table[[i]] <- unlist(ouch[[i]]$time)}
names(time.table) <- well

for (i in c(1:length(well))){
  prptn <- rep(NA,length(ouch[[i]]$time))
    for (j in c(1:length(prptn))){if (round(ouch[[i]][,5],5)[j]==tolik){prptn[j] <- 1}else{prptn[j] <- 0}}
  ouch[[i]] <- cbind(ouch[[i]],prptn)
}

proportion <- c(mean(ouch[[1]][,9]),mean(ouch[[2]][,9]),mean(ouch[[3]][,9]),mean(ouch[[4]][,9]),mean(ouch[[5]][,9]),mean(ouch[[6]][,9]),mean(ouch[[7]][,9]),mean(ouch[[8]][,9]),mean(ouch[[9]][,9]),mean(ouch[[10]][,9]),mean(ouch[[11]][,9]),mean(ouch[[12]][,9]))

time.tab <- matrix(NA,ncol=2,nrow=length(well))
for (i in c(1:length(well))){time.tab[i,] <- c(mean(time.table[[i]]),proportion[i])}

png(file="ouch_small_propvstime.png")
plot(time.tab[,1],time.tab[,2],pch=NA,xlab="Time",ylab="Proportion",main="ouch")
for (i in c(1:length(well))){
  text(x=time.tab[,1][i],y=time.tab[,2][i],labels=well[i])
}
dev.off()


#Results for aquilegia data
load("/home/michels/Hallowed/Dropbox/repository/phylooptim/pkg/R/geiger/aquigeigererror.RData")
geiger <- l
l2 <- NULL
l <- NULL
geiger$mean_lik <- round(c(mean(geiger[[1]][,2]),mean(geiger[[2]][,2]),mean(geiger[[3]][,2]),mean(geiger[[4]][,2]),mean(geiger[[5]][,2]),mean(geiger[[6]][,2]),mean(geiger[[7]][,2]),mean(geiger[[8]][,2]),mean(geiger[[9]][,2]),mean(geiger[[10]][,2]),mean(geiger[[11]][,2]),mean(geiger[[12]][,2])),5)
names(geiger$mean_lik) <- names(geiger)[1:12]
                 
load("/home/michels/Hallowed/Dropbox/repository/phylooptim/pkg/R/ace/aquiaceerror.RData")
ace <- l
l <- NULL
ace$mean_lik <- round(c(mean(ace[[1]][,5]),mean(ace[[2]][,5]),mean(ace[[3]][,5]),mean(ace[[4]][,5]),mean(ace[[5]][,5]),mean(ace[[6]][,5]),mean(ace[[7]][,5]),mean(ace[[8]][,5]),mean(ace[[9]][,5]),mean(ace[[10]][,5]),mean(ace[[11]][,5]),mean(ace[[12]][,5])),5)
names(ace$mean_lik) <- names(ace)[1:12]

load("/home/michels/Hallowed/Dropbox/repository/phylooptim/pkg/R/ouch/aquioucherror.RData")
ouch <- l
l <- NULL
ouch$mean_lik <- round(c(mean(ouch[[1]][,5]),mean(ouch[[2]][,5]),mean(ouch[[3]][,5]),mean(ouch[[4]][,5]),mean(ouch[[5]][,5]),mean(ouch[[6]][,5]),mean(ouch[[7]][,5]),mean(ouch[[8]][,5]),mean(ouch[[9]][,5]),mean(ouch[[10]][,5]),mean(ouch[[11]][,5]),mean(ouch[[12]][,5])),5)
names(ouch$mean_lik) <- names(ouch)[1:12]

#Find the mode
f <- function(x){
xt <- table(x)
t <- as.numeric(names(xt[xt == max(xt)]))  
return(t)}

#True Liklihood
tglik <- max(f(as.numeric(geiger$mean_lik)))
gliktf <- rep(NA,length(geiger$mean_lik))
for (i in c(1:length(gliktf))){if (as.numeric(geiger$mean_lik[i])==tglik){gliktf[i] <- "Y"}else{if (tglik-.05 < as.numeric(geiger$mean_lik[i]) && tglik+.05 > as.numeric(geiger$mean_lik[i])){gliktf[i] <- "O"}else{gliktf[i] <- "--"}}}

talik <- max(f(as.numeric(ace$mean_lik)))
aliktf <- rep(NA,length(ace$mean_lik))
for (i in c(1:length(aliktf))){if (as.numeric(ace$mean_lik[i])==talik){aliktf[i] <- "Y"}else{if (talik-.05 < as.numeric(ace$mean_lik[i]) && talik+.05 > as.numeric(ace$mean_lik[i])){aliktf[i] <- "O"}else{aliktf[i] <- "--"}}}

tolik <- max(f(as.numeric(ouch$mean_lik)))
oliktf <- rep(NA,length(ouch$mean_lik))
for (i in c(1:length(oliktf))){if (as.numeric(ouch$mean_lik[i])==tolik){oliktf[i] <- "Y"}else{if (tolik-.05 < as.numeric(ouch$mean_lik[i]) && tolik+.05 > as.numeric(ouch$mean_lik[i])){oliktf[i] <- "O"}else{oliktf[i] <- "--"}}}

smoke <- matrix(c(gliktf,aliktf,oliktf),ncol=12,byrow=TRUE)
colnames(smoke) <- names(geiger)[1:12]
rownames(smoke) <- c("geiger","ace","ouch")
aquilegia <- as.table(smoke)
aquilegia

#Plot geiger
time.table <- vector("list", length(well))
for (i in c(1:length(well))){
  time.table[[i]] <- unlist(geiger[[i]]$time)}
names(time.table) <- well

for (i in c(1:length(well))){
  prptn <- rep(NA,length(geiger[[i]]$time))
    for (j in c(1:length(prptn))){if (round(geiger[[i]][,2],5)[j]==tglik){prptn[j] <- 1}else{prptn[j] <- 0}}
  geiger[[i]] <- cbind(geiger[[i]],prptn)
}

proportion <- c(mean(geiger[[1]][,10]),mean(geiger[[2]][,10]),mean(geiger[[3]][,10]),mean(geiger[[4]][,10]),mean(geiger[[5]][,10]),mean(geiger[[6]][,10]),mean(geiger[[7]][,10]),mean(geiger[[8]][,10]),mean(geiger[[9]][,10]),mean(geiger[[10]][,10]),mean(geiger[[11]][,10]),mean(geiger[[12]][,10]))

time.tab <- matrix(NA,ncol=2,nrow=length(well))
for (i in c(1:length(well))){time.tab[i,] <- c(mean(time.table[[i]]),proportion[i])}

png(file="geiger_med_propvstime.png")
plot(time.tab[,1],time.tab[,2],pch=NA,xlab="Time",ylab="Proportion",main="fitContinuous")
for (i in c(1:length(well))){
  text(x=time.tab[,1][i],y=time.tab[,2][i],labels=well[i])
}
dev.off()

#Plot ace
time.table <- vector("list", length(well))
for (i in c(1:length(well))){
  time.table[[i]] <- unlist(ace[[i]]$time)}
names(time.table) <- well

for (i in c(1:length(well))){
  prptn <- rep(NA,length(ace[[i]]$time))
    for (j in c(1:length(prptn))){if (round(ace[[i]][,5],5)[j]==talik){prptn[j] <- 1}else{prptn[j] <- 0}}
  ace[[i]] <- cbind(ace[[i]],prptn)
}

proportion <- c(mean(ace[[1]][,7]),mean(ace[[2]][,7]),mean(ace[[3]][,7]),mean(ace[[4]][,7]),mean(ace[[5]][,7]),mean(ace[[6]][,7]),mean(ace[[7]][,7]),mean(ace[[8]][,7]),mean(ace[[9]][,7]),mean(ace[[10]][,7]),mean(ace[[11]][,7]),mean(ace[[12]][,7]))

time.tab <- matrix(NA,ncol=2,nrow=length(well))
for (i in c(1:length(well))){time.tab[i,] <- c(mean(time.table[[i]]),proportion[i])}

png(file="ace_med_propvstime.png")
plot(time.tab[,1],time.tab[,2],pch=NA,xlab="Time",ylab="Proportion",main="Ace")
for (i in c(1:length(well))){
  text(x=time.tab[,1][i],y=time.tab[,2][i],labels=well[i])
}
dev.off()

#Plot ouch
time.table <- vector("list", length(well))
for (i in c(1:length(well))){
  time.table[[i]] <- unlist(ouch[[i]]$time)}
names(time.table) <- well

for (i in c(1:length(well))){
  prptn <- rep(NA,length(ouch[[i]]$time))
    for (j in c(1:length(prptn))){if (round(ouch[[i]][,5],5)[j]==tolik){prptn[j] <- 1}else{prptn[j] <- 0}}
  ouch[[i]] <- cbind(ouch[[i]],prptn)
}

proportion <- c(mean(ouch[[1]][,9]),mean(ouch[[2]][,9]),mean(ouch[[3]][,9]),mean(ouch[[4]][,9]),mean(ouch[[5]][,9]),mean(ouch[[6]][,9]),mean(ouch[[7]][,9]),mean(ouch[[8]][,9]),mean(ouch[[9]][,9]),mean(ouch[[10]][,9]),mean(ouch[[11]][,9]),mean(ouch[[12]][,9]))

time.tab <- matrix(NA,ncol=2,nrow=length(well))
for (i in c(1:length(well))){time.tab[i,] <- c(mean(time.table[[i]]),proportion[i])}

png(file="ouch_med_propvstime.png")
plot(time.tab[,1],time.tab[,2],pch=NA,xlab="Time",ylab="Proportion",main="ouch")
for (i in c(1:length(well))){
  text(x=time.tab[,1][i],y=time.tab[,2][i],labels=well[i])
}
dev.off()
