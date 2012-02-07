#Required to run this part.
require(ouch)
require(optimx)
require(geiger)
load("ouchbackup.RData")

#Default inputs, if you want more iterations, change it, then rerun this piece of code.
#N <- detectCores()                       #Number of cores
N <- 22
it <- 4  	                         #Number of iterations
z <- length(well)                        #Number of optimizers
MIN=0.01                                 #min upper value
MAX=10                                   #max upper value
strt <- seq(MIN,MAX,length=it)           #Start values

a <- ceiling(it / N)
h <- 1
m <- 1
v <- vector("list",N)
for (b in c(1:N)){v[[b]]<-matrix(NA,nrow=a,ncol=1)}
for (d in c(1:a)){
  for (j in c(1:N)){
    if (m > it*z) {break}else{
    if (h <= 12){
      if (j < N ){v[[j]][d] <- c(strt[m]);h <- h+1;m <- m+1}else{v[[j]][d] <- c(strt[m]);h<-h+1;m <- m+1;next}
                 }else{
      h <- 1
      if (j < N ){v[[j]][d] <- c(strt[m]);h <- h+1;m <- m+1}else{v[[j]][d] <- c(strt[m]);h<-h+1;m <- m+1;next}}}
  }
}

for (b in c(1:N)){v[[b]] <- na.omit(v[[b]])}

#Mammal Data
tree<-read.nexus("mammalChar1.nex")
tree$edge.length[tree$edge.length<1e-5]=1e-5
trait <- read.csv("mammallogChar1.csv")
trait1 <- read.csv("mammallogChar1.csv", row.names=1)
nc <- name.check(tree,trait1)
tree <- drop.tip(tree,nc$Tree.not.data)

require(multicore)
b.time <-proc.time()
jobs <- lapply(v, function(x) parallel(f(x),silent=TRUE))
results <- collect(jobs,wait=TRUE)
total.time <- as.numeric(proc.time()[3]-b.time[3])/(60)

l <- vector("list",length(well))
names(l) <- well
for (i in c(1:length(well))){l[[i]] <- matrix(NA,ncol=8,nrow=N);names(l[[i]]) <- c("I","lb","ub","T","L","S","A","time")}

for (i in c(1:N)){for (j in c(1:length(well))){l[[j]][i,] <- as.numeric(results[[i]][[j]][1,])}}

save.image("/home/michels/repository/phylooptim/pkg/R/ouch/mamoucherror.RData")

#save.image("/Users/michels/phylooptim/pkg/R/geiger/mamoucherror.RData")

#save.image("/Users/omearalabguest/phylooptim/pkg/R/geiger/mamoucherror.RData")
