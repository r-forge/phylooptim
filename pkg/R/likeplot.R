rm(list=ls())
require(grid)

prop <- list()
load("/home/michels/Hallowed/repository/phylooptim/pkg/R/acediscresults.RData")
prop[[1]] <- o.table[[8]]
rm(list=ls()[-1*which(ls()=="prop")])
load("/home/michels/Hallowed/repository/phylooptim/pkg/R/geigerctsresults.RData")
prop[[2]] <- o.table[[8]]
rm(list=ls()[-1*which(ls()=="prop")])
load("/home/michels/Hallowed/repository/phylooptim/pkg/R/geigerdiscresults.RData")
prop[[3]] <- o.table[[8]]
rm(list=ls()[-1*which(ls()=="prop")])
load("/home/michels/Hallowed/repository/phylooptim/pkg/R/ouchctsresults.RData")
prop[[4]] <- o.table[[8]]
rm(list=ls()[-1*which(ls()=="prop")])
names(prop) <- c("acedisc","geigercts","geigerdisc","ouchcts")

well <- c("spg", "Rcgmin", "Rvmmin", "bobyqa","L-BFGS-B","nlminb","ucminf","Nelder-Mead","nlm","CG","BFGS","newuoa")

propnew <- matrix(NA,ncol=12,nrow=16)
for (k in c(1:4)){
  for (i in c(1:12)){
    for (j in c(1:4)){
if (k==1){propnew[j,i] <- abs(as.numeric(prop[[k]][j,i+1]))}
if (k==2){propnew[4+j,i] <- abs(as.numeric(prop[[k]][j,i+1]))}
if (k==3){propnew[8+j,i] <- abs(as.numeric(prop[[k]][j,i+1]))}
if (k==4){propnew[12+j,i] <- abs(as.numeric(prop[[k]][j,i+1]))}
    }
  }
}

#Identify large values
propnew[4,3] <- 1
propnew[3,3] <- 1
propnew[3,5] <- 1
propnew[4,5] <- 1
propnew[11,2] <- 1
propnew[11,3] <- 1
propnew[11,5] <- 1
propnew[12,5] <- 1
propnew[11,9] <- 1
propnew[11,10] <- 1
propnew[11,11] <- 1
propnew[12,9] <- 1
propnew[12,10] <- 1
propnew[12,11] <- 1

#propnew[,] <- 1
for(i in c(1:16)){
  for (j in c(1:12)){
    if (propnew[i,j]>0.25 & propnew[i,j]<0.75){propnew[i,j] <- 0.50}
  }
}
for(i in c(1:16)){
  for (j in c(1:12)){
    if (propnew[i,j]>0 & propnew[i,j]<0.05){propnew[i,j] <- 0.20}
  }
}

png("likehoriplot.png")
plot(x=c(-1,12),y=c(0,17.75),type="n",axes=FALSE,xlab="",ylab="",main="Closeness of MLE")
for (y in c(1:4)) {
  for (x in c(1:12))  {
     rect(xleft=x-1,xright=x,ybottom=17-y,ytop=17-(y+1),col=gray(1-propnew[y,x]))
   }
}

for (y in c(5:8)) {
  for (x in c(1:12))  {
     rect(xleft=x-1,xright=x,ybottom=17-y-0.25,ytop=17-(y+1)-0.25,col=gray(1-propnew[y,x]))
   }
}

for (y in c(9:12)) {
  for (x in c(1:12))  {
     rect(xleft=x-1,xright=x,ybottom=17-y-0.5,ytop=17-(y+1)-0.5,col=gray(1-propnew[y,x]))
   }
}

for (y in c(13:16)) {
  for (x in c(1:12))  {
     rect(xleft=x-1,xright=x,ybottom=17-y-0.75,ytop=17-(y+1)-0.75,col=gray(1-propnew[y,x]))
   }
}

text(.5, 16.5, "spg",cex = .8,srt=45)
text(1.5, 16.8, "Rcgmin",cex = .8,srt=45)
text(2.5, 16.8, "Rvmmin",cex = .8,srt=45)
text(3.5, 16.8, "bobyqa",cex = .8,srt=45)
text(4.7, 17.05, "L-BFGS-B",cex = .8,srt=45)
text(5.5, 16.8, "nlminb",cex = .8,srt=45)
text(6.5, 16.8, "ucminf",cex = .8,srt=45)
text(7.8, 17.2, "Nelder-Mead",cex = .8,srt=45)
text(8.5, 16.5, "nlm",cex = .8,srt=45)
text(9.5, 16.5, "CG",cex = .8,srt=45)
text(10.5, 16.7, "BFGS",cex = .8,srt=45)
text(11.5, 16.8, "newuoa",cex = .8,srt=45)

text(-.5, 15.5, "Geo",cex = .8)
text(-.5, 14.5, "Aqui",cex = .8)
text(-.5, 13.5, "Mono",cex = .8)
text(-.5, 12.5, "Mam",cex = .8)

text(-.5, 11.5-.25, "Geo",cex = .8)
text(-.5, 10.5-.25, "Aqui",cex = .8)
text(-.5, 9.5-.25, "Mono",cex = .8)
text(-.5, 8.5-.25, "Mam",cex = .8)

text(-.5, 7.5-.5, "Geo",cex = .8)
text(-.5, 6.5-.5, "Aqui",cex = .8)
text(-.5, 5.5-.5, "Mono",cex = .8)
text(-.5, 4.5-.5, "Mam",cex = .8)

text(-.5, 3.5-.75, "Geo",cex = .8)
text(-.5, 2.5-.75, "Aqui",cex = .8)
text(-.5, 1.5-.75, "Mono",cex = .8)
text(-.5, 0.5-.75, "Mam",cex = .8)

text(-1.25, 14, "Ace Discrete",cex = .8,srt=90)
text(-1.25, 10-.25, "fitContinuous",cex = .8,srt=90)
text(-1.25, 6-.5, "fitDiscrete",cex = .8,srt=90)
text(-1.25, 2-.75, "Hansen",cex = .8,srt=90)

dev.off()

#This is the way for the paper

propnew <- matrix(NA,ncol=16,nrow=12)
for (k in c(1:4)){
  for (i in c(1:12)){
    for (j in c(1:4)){
if (k==1){propnew[i,j] <- abs(as.numeric(prop[[k]][j,i+1]))}
if (k==2){propnew[i,4+j] <- abs(as.numeric(prop[[k]][j,i+1]))}
if (k==3){propnew[i,8+j] <- abs(as.numeric(prop[[k]][j,i+1]))}
if (k==4){propnew[i,12+j] <- abs(as.numeric(prop[[k]][j,i+1]))}
    }
  }
}

row.names(propnew) <- well
propnew <- cbind(propnew,c(5,6,9,7,10,2,1,3,8,11,12,4))

propnew <- propnew[order(propnew[,17]) , ][,c(1:16)]

for(i in c(1:12)){
  for (j in c(1:16)){
    if (propnew[i,j]<2 & propnew[i,j]>0){propnew[i,j] <- 0.50}
  }
}
for(i in c(1:12)){
  for (j in c(1:16)){
    if (propnew[i,j]>=2){propnew[i,j] <- 1}
  }
}

png("likevertplot.png")
plot(x=c(-1,17),y=c(-1,12),type="n",axes=FALSE,xlab="",ylab="")
for (x in c(1:4)) {
  for (y in c(1:12))  {
     rect(xleft=x,xright=x+1,ybottom=12-(y+1),ytop=12-y,col=gray(propnew[y,x]))  
   }
}

for (x in c(5:8)) {
  for (y in c(1:12))  {
     rect(xleft=x+0.25,xright=x+1+0.25,ybottom=12-(y+1),ytop=12-y,col=gray(propnew[y,x]))  
   }
}

for (x in c(9:12)) {
  for (y in c(1:12))  {
     rect(xleft=x+0.5,xright=x+1+0.5,ybottom=12-(y+1),ytop=12-y,col=gray(propnew[y,x]))    
   }
}

for (x in c(13:16)) {
  for (y in c(1:12))  {
     rect(xleft=x+0.75,xright=x+1+0.75,ybottom=12-(y+1),ytop=12-y,col=gray(propnew[y,x]))    
   }
}

text(0.5, 10.5, "ucminf",cex = .8)
text(0, 9.5, "nlminb",cex = .8)
text(0, 8.5, "Nelder-Mead",cex = .8)
text(0.1, 7.5, "bobyqa",cex = .8)
text(-0.25, 6.5, "L-BFGS-B",cex = .8)
text(0.15, 5.5, "nlminb",cex = .8)
text(0.15, 4.5, "ucminf",cex = .8)
text(-0.45, 3.5, "Nelder-Mead",cex = .8)
text(0.5, 2.5, "nlm",cex = .8)
text(0.5, 1.5, "CG",cex = .8)
text(0.2, 0.5, "BFGS",cex = .8)
text(0.1, -0.5, "newuoa",cex = .8)

text(1.5, 11.5, "Geo",cex = .8,srt=90)
text(2.5, 11.55, "Aqui",cex = .8,srt=90)
text(3.5, 11.6, "Mono",cex = .8,srt=90)
text(4.5, 11.55, "Mam",cex = .8,srt=90)

text(5.5+.25, 11.5, "Geo",cex = .8,srt=90)
text(6.5+.25, 11.55, "Aqui",cex = .8,srt=90)
text(7.5+.25, 11.6, "Mono",cex = .8,srt=90)
text(8.5+.25, 11.55, "Mam",cex = .8,srt=90)

text(9.5+.5, 11.5, "Geo",cex = .8,srt=90)
text(10.5+.5, 11.55, "Aqui",cex = .8,srt=90)
text(11.5+.5, 11.6, "Mono",cex = .8,srt=90)
text(12.5+.5, 11.55, "Mam",cex = .8,srt=90)

text(13.5+.75, 11.5, "Geo",cex = .8,srt=90)
text(14.5+.75, 11.55, "Aqui",cex = .8,srt=90)
text(15.5+.75, 11.6, "Mono",cex = .8,srt=90)
text(16.5+.75, 11.55, "Mam",cex = .8,srt=90)

text(3, 12.35, "Ace Discrete",cex = .8)
text(7+.25, 12.35, "fitContinuous",cex = .8)
text(11+.5, 12.35, "fitDiscrete",cex = .8)
text(15+.75, 12.35, "Hansen",cex = .8)

dev.off()
