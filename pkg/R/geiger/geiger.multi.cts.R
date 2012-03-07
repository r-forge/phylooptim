#Make sure to have geigerbackupcts.RData in the correct directory
rm(list = ls())
require(geiger)
require(optimx)
require(ape)
require(multicore)
source("geiger.mine.cts.R")
source("f.geiger.cts.R")

wellt <- c("spg", "Rcgmin", "Rvmmin", "bobyqa","L-BFGS-B",1,1,1,1,1,1,1)
well <- c("spg", "Rcgmin", "Rvmmin", "bobyqa","L-BFGS-B","nlminb","ucminf","Nelder-Mead","nlm","CG","BFGS","newuoa")

#The following will detect the number of cores you have.
if(.Platform$OS.type == "windows") {
detectCores <- function(all.tests = FALSE)
    .Call("ncpus", FALSE, package = "parallel")
} else {
detectCores <- function(all.tests = FALSE)
{
    systems <-
        list(darwin  = "/usr/sbin/sysctl -n hw.ncpu 2>/dev/null",
             freebsd = "/sbin/sysctl -n hw.ncpu 2>/dev/null",
             linux   = "grep processor /proc/cpuinfo 2>/dev/null | wc -l",
             irix    = c("hinv | grep Processors | sed 's: .*::'",
                         "hinv | grep '^Processor '| wc -l"),
             solaris = "/usr/sbin/psrinfo -p")
    for (i in seq(systems))
        if(all.tests ||
           length(grep(paste("^", names(systems)[i], sep=''), R.version$os)))
            for (cmd in systems[i]) {
                a <- gsub("^ +","", system(cmd, TRUE)[1])
                if (length(grep("^[1-9]", a))) return(as.integer(a))
            }
    NA_integer_
}
}

M <- detectCores()                      #Use if you want all cores used (Total # of cores)
#M <- 2                                   #Use if you want to choose the # of cores to use
it <- 50  	                         #Number of iterations (at least 2)
z <- length(well)                        #Number of optimizers
MIN=2                                    #min upper value
MAX=500                                  #max upper value
ub <- seq(MIN,MAX,length=it)             #Upper bound values
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

##Geospiza Data
data(geospiza)
name <- c("geo")
start=log(c(0.1, 1.5))
tree <- geospiza$geospiza.tree
a.trait <- geospiza$geospiza.data

##Which column of data do you want?
kk <- 1

##Aquilegia Data
#name <- c("aqui")
#start=log(c(.1, 1.5))
#tree <- read.tree("Aquilegia.new.tre")
#a.trait <- read.delim("Aquilegia.traits",row.names=1)

##Which column of data do you want?
#kk <- 2

##Monocot Data
#name <- c("mono")
#start=log(c(0.1, 1.5))
#tree <- read.tree("BJO.Monocot.tre")
#a.trait <- read.delim("BJO.monocot_GS",row.names=1)

##Which column of data do you want?
#kk <- 2

##Mammal Data
#name <- c("mam")
#start=log(c(1.5, 0.1))
#tree<-read.nexus("mammalChar1.nex")
#a.trait <- read.csv("mammallogChar1.csv", row.names=1)

##Which column of data do you want?
#kk <- 1

if (dim(a.trait)[2] == 1){a.trait <- data.frame(a.trait,well=rep(1,length(a.trait)))}
nc <- name.check(tree,a.trait)
if (nc[1]=="OK"){nc$Tree.not.data <- NULL}
tree <- drop.tip(tree,nc$Tree.not.data)
tree$edge.length[tree$edge.length<1e-5]=1e-5

a.trait <- a.trait[order(rownames(a.trait)),]
d <- data.frame(name=tree$tip.label,num=seq(1,length(tree$tip.label),by=1))
dd <- d[order(d[,1]) , ]
a.trait$new <- dd[,2]
a.trait <- a.trait[order(a.trait$new),]

td <- treedata(tree,a.trait[,kk],sort=T)

chdata<- as.vector(td$data) # TIP data
ntax=length(td$phy$tip.label)
tree<- td$phy# Tree
n <- length(chdata)

dplot<-TRUE
foostore=list(x1=NULL,x2=NULL)
  print(start)

b.time <-proc.time()
jobs <- lapply(v, function(x) parallel(f(x),silent=TRUE))
results <- collect(jobs,wait=TRUE)
total.time <- as.numeric(proc.time()[3]-b.time[3])/(60)

lt <- vector("list",length(well))

h <- 1
  for (i in c(1:N)){
    for (j in c(1:length(v[[i]][,1]))){
lt[[h]]<-c(unlist(v[[i]][j]),round(as.numeric(results[[i]][[j]]$lnl),5),round(as.numeric(results[[i]][[j]]$beta.start),5),round(as.numeric(exp(results[[i]][[j]]$beta.bnd[1])),5),round(exp(as.numeric(results[[i]][[j]]$beta.bnd[2])),5),round(results[[i]][[j]]$beta,5),round(as.numeric(results[[i]][[j]]$delta.start),5),round(exp(as.numeric(results[[i]][[j]]$delta.bnd[1])),5),round(exp(as.numeric(results[[i]][[j]]$delta.bnd[2])),5),round(results[[i]][[j]]$delta,5),as.numeric(results[[i]][[j]]$conv),round(results[[i]][[j]]$time,5));h <- h+1}}

ff <- data.frame(matrix(NA,ncol=length(lt[[1]]),nrow=length(lt)))
colnames(ff)<-c("name","L","bLB","bUB","b","dLB","dUB","d","conv","time")
for (i in c(1:length(lt))){
  for (j in c(1:length(lt[[i]]))){
    ff[i,j] <- lt[[i]][j]}}
ff <- ff[order(ff$name,ff$dUB) , ]

ll <- vector("list",length(well))
names(ll) <- c("BFGS","bobyqa","CG","L-BFGS-B","Nelder-Mead","newuoa","nlm","nlminb","Rcgmin","Rvmmin","spg","ucminf")

aa <- seq(1,it*z,by=it)
bb <- seq(aa[2]-1,it*z,by=it)
for (i in c(1:length(aa))){ll[[i]] <- ff[aa[i]:bb[i],2:10]}

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

#Saves the results in the directory
save.image(paste(getwd(),"/",name,"geigercts.RData",sep=""))
