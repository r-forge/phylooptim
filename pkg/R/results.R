rm(list = ls())

require(gplots)

well <- c("spg", "Rcgmin", "Rvmmin", "bobyqa","L-BFGS-B","nlminb","ucminf","Nelder-Mead","nlm","CG","BFGS","newuoa")

diff <- function(x, j, y, prefix="")
{
    v <- as.numeric(x$MLE[j])-x$'Overall MLE'
    txt <- format(c(v, 0.123456789), digits=y)
    q <- paste(prefix, txt, sep="")[1:12]
    return(q)
}

MLE <- function(x, y, prefix="")
{
    r <- x$MLE
    txt <- format(c(r, 0.123456789), digits=y)
    txt <- paste(prefix, txt, sep="")[1:12]
#   if(missing(cex.cor)) cex.cor <- 0.5/strwidth(txt)
#   text(0.5, 0.5, txt, cex = cex.cor )
    return(txt)
}

sd2 <- function(x, j)
{
    if (is.na(x$sd_lik[j])){q <- "NA"}else{
    q <- x$sd_lik[j]}
    return(q)
}

max1 <- function(x, y, prefix="")
{
    t <- na.omit(as.numeric(x[,5]))
    q <- round(max(t),y)
    txt <- format(c(q, 0.123456789), digits=y)
    txt <- paste(prefix, txt, sep="")[1:12]
    return(txt)
}

table1 <- function(x, j, y, prefix="")
{
    if (is.na(x$mean_lik[j])){q <- "NA"}else{
    if (round(as.numeric(x$mean_lik[j]),y)==x$'Overall MLE'){q <- "Y"}else{
    v <- x$mean_lik[j]-x$'Overall MLE'
    txt <- format(c(v, 0.123456789), digits=y)
    q <- paste(prefix, txt, sep="")[1:12]}}
    return(q)
}

table2 <- function(x, j, y, prefix="")
{
    if (is.na(x$median_lik[j])){q <- "NA"}else{
    if (round(as.numeric(x$median_lik[j]),y)==x$'Overall MLE'){q <- "Y"}else{
    v <- x$median_lik[j]-x$'Overall MLE'
    txt <- format(c(v, 0.123456789), digits=y)
    q <- paste(prefix, txt, sep="")[1:12]}}
    return(q)
}

tablefinal <- function(y,j)
{
table <- matrix(c(mean_likelihood[j,],median_likelihood[j,],proportion[j,],MLE_table[j,],c("",sd_table[j,])),ncol=13,byrow=TRUE)
colnames(table) <- colnames(proportion)
rownames(table) <- c("Mean MLE diff from true MLE ","Median MLE diff from true MLE ","Proportion","MLE","SD of MLE")
table <- as.table(table)
names(table) <- y$name1
return(table)
}

plot1 <- function(yy,k)
{
time.table <- vector("list", length(well))
for (i in c(1:length(well))){
  time.table[[i]] <- unlist(as.numeric(yy[[i]][,6]))*60}
names(time.table) <- well

time.tab <- matrix(NA,ncol=2,nrow=length(well))
MIN <- rep(NA,length(well))
MAX <- rep(NA,length(well))
UI <- rep(NA,length(well))
LI <- rep(NA,length(well))
for (i in c(1:length(well))){
  time.tab[i,] <- c(mean(time.table[[i]]),as.numeric(proportion[k,])[i+1])
  MIN[i] <- min(time.table[[i]])
  MAX[i] <- max(time.table[[i]])
  UI[i] <- mean(time.table[[i]])+qnorm(.975)*sd(time.table[[i]])
  LI[i] <- mean(time.table[[i]])+qnorm(.025)*sd(time.table[[i]])
}

png(file=paste(yy$name2,"_propvstime.png",sep=""))
plotCI(time.tab[,1],time.tab[,2],ui=UI,li=LI,xlim=c(min(MIN),max(MAX)),err="x",xlab="Time (Seconds)",ylab="Proportion",main=paste(yy$name1),pch=NA,gap=3)
for (i in c(1:length(well)))
  {
  text(x=time.tab[,1][i],y=time.tab[,2][i],labels=well[i])
}
dev.off()

q <- tablefinal(yy,k)

return(q)
}

meanprptn <- function(x,j,y)
{
    t <- as.numeric(x[[j]][,8])
    q <- round(mean(t),y)
    return(q)
}

#Discrete Geospiza Data
load("/home/michels/Hallowed/repository/phylooptim/pkg/R/ace/geoacedisc.RData")
ace_geo <- l
l <- NULL

#Discrete Aquilegia Data
load("/home/michels/Hallowed/repository/phylooptim/pkg/R/ace/aquiacedisc.RData")
ace_aqui <- l
l <- NULL

#Discrete Monocot Data
load("/home/michels/Hallowed/repository/phylooptim/pkg/R/ace/monoacedisc.RData")
ace_mono <- l
l <- NULL

#Discrete Mammal Data
load("/home/michels/Hallowed/repository/phylooptim/pkg/R/ace/mamacedisc.RData")
ace_mam <- l
l <- NULL

#Result tables
didu <- list()
didu[[1]] <- ace_geo
didu[[2]] <- ace_aqui
didu[[3]] <- ace_mono
didu[[4]] <- ace_mam

mean_lik <- vector("list",length(didu))
for (i in c(1:length(didu))){mean_lik[[i]] <- c(didu[[i]]$'Overall MLE',rep(NA,length(well)))}
names(mean_lik) <- c(didu[[1]]$name3,didu[[2]]$name3,didu[[3]]$name3,didu[[4]]$name3)
median_lik <- vector("list",length(didu))
for (i in c(1:length(didu))){median_lik[[i]] <- c(didu[[i]]$'Overall MLE',rep(NA,length(well)))}
names(median_lik) <- c(didu[[1]]$name3,didu[[2]]$name3,didu[[3]]$name3,didu[[4]]$name3)
sd_lik <- vector("list",length(didu))
for (i in c(1:length(didu))){sd_lik[[i]] <- c("",rep(NA,length(well)))}
names(sd_lik) <- c(didu[[1]]$name3,didu[[2]]$name3,didu[[3]]$name3,didu[[4]]$name3)
prop <- vector("list",length(didu))
for (i in c(1:length(didu))){prop[[i]] <- c(didu[[i]]$'Overall MLE',rep(NA,length(well)))}
names(prop) <- c(didu[[1]]$name3,didu[[2]]$name3,didu[[3]]$name3,didu[[4]]$name3)
MLEt <- vector("list",length(didu))
for (i in c(1:length(didu))){MLEt[[i]] <- c(didu[[i]]$'Overall MLE',rep(NA,length(well)))}
names(MLEt) <- c(didu[[1]]$name3,didu[[2]]$name3,didu[[3]]$name3,didu[[4]]$name3)


for (a in c(1:length(didu))){
   for (i in c(1:length(well))){
   prop[[a]][i+1] <- meanprptn(didu[[a]],i,4)
   }
      for (i in c(1:length(well))){
   mean_lik[[a]][i+1] <- table1(didu[[a]],i,4)
   }
   for (i in c(1:length(well))){
   median_lik[[a]][i+1] <- table2(didu[[a]],i,4)
   }
   for (i in c(1:length(well))){
   sd_lik[[a]][i+1] <- sd2(didu[[a]],i)
   }
   for (i in c(1:length(well))){
   MLEt[[a]][i+1] <- diff(didu[[a]],i,4)
   }
}

smoke <- matrix(c(mean_lik$geo,mean_lik$aqui,mean_lik$mono,mean_lik$mam),ncol=13,byrow=TRUE)
colnames(smoke) <- c("-2*log Lklhood",well)
rownames(smoke) <- c("Geospiza","Aquilegia","Monocot","Mammal")
mean_likelihood <- as.table(smoke)

smoke4 <- matrix(c(median_lik$geo,median_lik$aqui,median_lik$mono,median_lik$mam),ncol=13,byrow=TRUE)
colnames(smoke4) <- c("-2*log Lklhood",well)
rownames(smoke4) <- c("Geospiza","Aquilegia","Monocot","Mammal")
median_likelihood <- as.table(smoke4)

smoke3 <- matrix(c(sd_lik$geo[2:13],sd_lik$aqui[2:13],sd_lik$mono[2:13],sd_lik$mam[2:13]),ncol=12,byrow=TRUE)
colnames(smoke3) <- well
rownames(smoke3) <- c("Geospiza","Aquilegia","Monocot","Mammal")
sd_table <- as.table(smoke3)

smoke1 <- matrix(c(prop$geo,prop$aqui,prop$mono,prop$mam),ncol=13,byrow=TRUE)
colnames(smoke1) <- c("-2*log Lklhood",well)
rownames(smoke1) <- c("Geospiza","Aquilegia","Monocot","Mammal")
proportion <- as.table(smoke1)

smoke2 <- matrix(c(MLEt$geo,MLEt$aqui,MLEt$mono,MLEt$mam),ncol=13,byrow=TRUE)
colnames(smoke2) <- c("-2*log Lklhood",well)
rownames(smoke2) <- c("Geospiza","Aquilegia","Monocot","Mammal")
MLE_table <- as.table(smoke2)

o.table <- list()
for (j in c(1:length(didu))){
  o.table[[j]] <- plot1(didu[[j]],j)
}
names(o.table) <- c(didu[[1]]$name1,didu[[2]]$name1,didu[[3]]$name1,didu[[4]]$name1)

o.table$'Median -2*log Lklhood difference from true -2*log Lklhood' <- median_likelihood
o.table$'Mean -2*log Lklhood difference from true -2*log Lklhood' <- mean_likelihood
o.table$'Proportion Correct Lklhood' <- proportion
o.table$'Difference from true -2*log Lklhood' <- MLE_table
o.table$'SD of each -2*log Lklhood' <- sd_table
o.table

rm(list=ls()[-1*c(which(ls()=="o.table"),which(ls()=="ace_geo"),which(ls()=="ace_aqui"),which(ls()=="ace_mono"),which(ls()=="ace_mam"))])

save.image("/home/michels/Hallowed/repository/phylooptim/pkg/R/acediscresults.RData")
