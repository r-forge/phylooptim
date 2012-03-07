#Make sure to have acebackupcts.RData in the correct directory
rm(list = ls())
require(ouch)
require(optimx)
require(geiger)
require(multicore)
source("ace.mine.cts.R")
source("f.ace.cts.R")

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
#M <- 2                                  #Use if you want to choose the # of cores to use
it <- 10  	                         #Number of iterations (at least 2)
z <- length(well)                        #Number of optimizers
MIN=0.01                                 #min upper value
MAX=10                                   #max upper value
strt <- seq(MIN,MAX,length=it)           #Start values
N <- min(M,it)

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

oli <- rep(NA,N)
sum.oli <- rep(NA,N)
for (b in c(1:N)){v[[b]] <- na.omit(v[[b]])}
for (b in c(1:N)){oli[b] <- length(v[[b]])}
for (b in c(1:N)){sum.oli[b] <- sum(oli[b:1])}

#DATA TO BE USED (Default is geospiza)
##Geospiza data
data(geospiza)
name <- c("geo")
tree <- geospiza$geospiza.tree
a.trait <- geospiza$geospiza.data

##Which column of data do you want?
kk <- 1

##Aquilegia Data
#name <- c("aqui")
#tree <- read.tree("Aquilegia.new.tre")
#a.trait <- read.delim("Aquilegia.traits",row.names=1)

##Which column of data do you want?
#kk <- 1

##Monocot Data
#name <- c("mono")
#tree <- read.tree("BJO.Monocot.tre")
#a.trait <- read.delim("BJO.monocot_GS",row.names=1)

##Which column of data do you want?
#kk <- 3

##Mammal Data
#name <- c("mam")
#tree<-read.nexus("mammalChar4.nex")
#a.trait <- read.csv("mammalChar4.csv", row.names=1)

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

dv <- treedata(tree,a.trait[,kk],sort=T)

b.time <-proc.time()
jobs <- lapply(v, function(x) parallel(f(x),silent=TRUE))
results <- collect(jobs,wait=TRUE)
total.time <- as.numeric(proc.time()[3]-b.time[3])/(60)

l <- vector("list",length(well))
names(l) <- well
for (i in c(1:length(well))){l[[i]] <- matrix(NA,ncol=7,nrow=it);colnames(l[[i]]) <- c("I","lb","ub","P","L","time","conv")}

for (j in c(1:N)){
  for (i in c(1:z)){
    for (d in c(1:oli[j])){
      if (j == 1 ){l[[i]][d,] <- as.numeric(results[[j]][[i]][d,])}else{
                   l[[i]][d+sum.oli[j-1],] <- as.numeric(results[[j]][[i]][d,])}
    }
  }
}

for (i in c(1:z)){l[[i]] <- l[[i]][order(l[[i]][,1]) , ]}

#Saves the results in the directory
save.image(paste(getwd(),"/",name,"acects.RData",sep=""))
