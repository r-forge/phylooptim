rm(list = ls())

load("/home/michels/Hallowed/repository/phylooptim/pkg/R/ouch/mamouchcts.RData")
l[[1]]
#What colummn is likelihood in?
lc <- 5

mean1 <- function(x, y)
{
    t <- na.omit(as.numeric(x[,lc]))
    q <- round(mean(t),y)
    return(q)
}

median1 <- function(x, y)
{
    t <- na.omit(as.numeric(x[,lc]))
    q <- round(median(t),y)
    return(q)
}

sd1 <- function(x, y)
{
    t <- na.omit(as.numeric(x[,lc]))
    q <- round(sd(t),y)
    return(q)
}

max1 <- function(x, y, prefix="")
{
    t <- na.omit(as.numeric(x[,lc]))  
    q <- round(max(t),y)
    txt <- format(c(q, 0.123456789), digits=y)
    txt <- paste(prefix, txt, sep="")[1:12]
    return(txt)
}

CI1 <- function(x,j,y)
{
    q <- round(as.numeric(x$mean_lik[j])+c(-1,1)*qnorm(.975)*as.numeric(x$sd_lik[j]),y)
    return(q)
}

oMLE <- function(x,y)
{
    t <- na.omit(as.numeric(x$MLE))
    q <- round(max(t),y)
    return(q)
}

prptn1 <- function(x,y)
{
for (i in c(1:length(well))){
  prptn <- rep(NA,length(x[[i]][,lc]))
    for (j in c(1:length(prptn))){
      if (is.na(x[[i]][,lc][j])){prptn[j] <- 0}else{
        if (round(as.numeric(x[[i]][,lc]),y)[j]==oMLE(x,y)){prptn[j] <- 1}else{
          prptn[j] <- 0}}}
  x[[i]] <- cbind(x[[i]],prptn)}
  return(x)
}

well <- c("spg", "Rcgmin", "Rvmmin", "bobyqa","L-BFGS-B","nlminb","ucminf","Nelder-Mead","nlm","CG","BFGS","newuoa")

l$name1 <- "Hansen Cts Mammal"
l$name2 <- "ouch_mam_disc"
l$name3 <- "mam"
l$mean_lik <- rep(NA,length(well))
l$median_lik <- rep(NA,length(well))
l$sd_lik <- rep(NA,length(well))
l$MLE <- rep(NA,length(well))
l$CI_lik <- list()
for (i in c(1:length(well))){
l$mean_lik[i] <- mean1(l[[i]],4)
l$median_lik[i] <- median1(l[[i]],4)
l$sd_lik[i] <- sd1(l[[i]],4)
l$MLE[i] <- max1(l[[i]],4)
l$CI_lik[[i]] <- CI1(l,i,4)
}
names(l$mean_lik) <- names(l)[1:12]
names(l$median_lik) <- names(l)[1:12]
names(l$sd_lik) <- names(l)[1:12]
names(l$MLE) <- names(l)[1:12]
names(l$CI_lik) <- names(l)[1:12]

l <- prptn1(l,4)

l$'Overall MLE' <- oMLE(l,4)

rm(list=ls()[-1*c(which(ls()=="a.trait"),which(ls()=="l"),which(ls()=="kk"),which(ls()=="dv"))])

save.image("/home/michels/Hallowed/repository/phylooptim/pkg/R/ouch/mamouchcts.RData")
