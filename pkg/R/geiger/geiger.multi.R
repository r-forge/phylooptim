require(geiger)
require(optimx)
require(ape)
require(multicore)
load("geigerbackup.RData")

#Run this
M <- 22
it <- 50  	                         #Number of iterations
z <- length(well)                        #Number of optimizers
MIN=2                                    #min upper value
MAX=500                                  #max upper value
ub <- seq(MIN,MAX,length=it)             #Random beta start values
N <- min(it*z,M)

bb <- unlist(lapply(ub, function(x) rep(x,z)))
a <- ceiling(it*z / N)
h <- 1
m <- 1
v <- vector("list",N)
for (b in c(1:N)){v[[b]]<-matrix(NA,nrow=a,ncol=2)}
for (d in c(1:a)){
  for (j in c(1:N)){
    if (m > it*z) {break}else{
    if (h <= 12){
      if (j < N ){v[[j]][d,] <- c(well[h],bb[m]);h <- h+1;m <- m+1}else{v[[j]][d,] <- c(well[h],bb[m]);h<-h+1;m <- m+1;next}
                 }else{
      h <- 1
      if (j < N ){v[[j]][d,] <- c(well[h],bb[m]);h <- h+1;m <- m+1}else{v[[j]][d,] <- c(well[h],bb[m]);h<-h+1;m <- m+1;next}}}
  }
}

for (b in c(1:N)){v[[b]] <- na.omit(v[[b]])}

#td <- treedata(read.tree("BJO.Monocot.tre"),read.delim("BJO.monocot_GS")[,3],sort=T)
tree <- read.nexus("mammalChar1.nex")
tree$edge.length[tree$edge.length<1e-5]=1e-5
td <- treedata(tree,read.csv("mammallogChar1.csv", row.names=1),sort=T)
ntax=length(td$phy$tip.label)
chdata <- read.csv("mammallogChar1.csv")[,2]
tree<- td$phy# Tree
n<- length(chdata)

#----- MINIMIZE NEGATIVE LOG LIKELIHOOD

start=log(c(1.5, 0.1))

dplot<-TRUE
foostore=list(x1=NULL,x2=NULL)
  print(start)

b.time <-proc.time()
jobs <- lapply(v, function(x) parallel(f(x),silent=TRUE))
results <- collect(jobs,wait=TRUE)
total.time <- as.numeric(proc.time()[3]-b.time[3])/(60)
#what?

lt <- vector("list",length(well))

h <- 1
  for (i in c(1:N)){
    for (j in c(1:length(v[[i]][,1]))){
lt[[h]]<-c(unlist(v[[i]][j]),round(as.numeric(results[[i]][[j]]$lnl),5),round(as.numeric(results[[i]][[j]]$beta.start),5),round(as.numeric(exp(results[[i]][[j]]$beta.bnd[1])),5),round(exp(as.numeric(results[[i]][[j]]$beta.bnd[2])),5),round(results[[i]][[j]]$beta,5),round(as.numeric(results[[i]][[j]]$delta.start),5),round(exp(as.numeric(results[[i]][[j]]$delta.bnd[1])),5),round(exp(as.numeric(results[[i]][[j]]$delta.bnd[2])),5),round(results[[i]][[j]]$delta,5),round(results[[i]][[j]]$time,5));h <- h+1}}

ff <- data.frame(matrix(NA,ncol=length(lt[[1]]),nrow=length(lt)))
colnames(ff)<-c("name","L","bLB","bUB","b","dLB","dUB","d","time")
for (i in c(1:length(lt))){
  for (j in c(1:length(lt[[i]]))){
    ff[i,j] <- lt[[i]][j]}}
ff <- ff[order(ff$name,ff$dUB) , ]

ll <- vector("list",length(well))
names(ll) <- c("BFGS","bobyqa","CG","L-BFGS-B","Nelder-Mead","newuoa","nlm","nlminb","Rcgmin","Rvmmin","spg","ucminf")

aa <- seq(1,it*z,by=it)
bb <- seq(aa[2]-1,it*z,by=it)
for (i in c(1:length(aa))){ll[[i]] <- ff[aa[i]:bb[i],2:9]}

l <- vector("list",length(well))
names(l) <- well
l[[1]] <- ll[[11]]
l[[2]] <- ll[[9]]
l[[3]] <- ll[[10]]
l[[4]] <- ll[[2]]
l[[5]] <- ll[[4]]
l[[6]] <- ll[[8]]
l[[7]] <- ll[[12]]
l[[8]] <- ll[[5]]
l[[9]] <- ll[[7]]
l[[10]] <- ll[[3]]
l[[11]] <- ll[[1]]
l[[12]] <- ll[[6]]

save.image("/Users/michels/phylooptim/pkg/R/geiger/mamgeigererror.RData")

save.image("/Users/omearalabguest/phylooptim/pkg/R/geiger/mamgeigererror.RData")
