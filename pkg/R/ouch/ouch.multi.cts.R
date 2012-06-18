#Make sure to have ouchbackupcts.RData in the correct directory,
rm(list = ls())
require(ouch)
rm(hansen)
require(optimx)
require(geiger)
require(multicore)
#load("ouchbackupcts.RData")
source("ouch.error2.R")
source("f.ouch.cts.R")
source("ouch.wrap2.x.R")

well <- c("spg", "Rcgmin", "Rvmmin", "bobyqa","L-BFGS-B","nlminb","ucminf","Nelder-Mead","nlm","CG","BFGS","newuoa")
wellt <- c("spg", "Rcgmin", "Rvmmin", "bobyqa","L-BFGS-B",1,1,1,1,1,1,1)

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

M <- detectCores()                       #Use if you want all cores used (Total # of cores)
#M <- 2                                  #Use if you want to choose the # of cores to use
#it <- 2  	                         #Number of iterations (at least 2)
iter <- 500                              #Number of max itn's (in function)
z <- length(well)                        #Number of optimizers
#MIN=0.01                                #min upper value
#MAX=10                                  #max upper value
#strt <- seq(MIN,MAX,length=it)          #Start values
#strt <- c(0.4857143,4.7671429)
strt <- list()
strt[[1]] <- c(0.5,1.5,2.5,3.5,4.5)
strt[[2]] <- c(0.5,1.5,2.5,3.5,4.5)
names(strt) <- c("sqrt.alpha","sigma")
it <- length(strt[[1]])

N <- min(M,it)

a <- ceiling(it / N)
m <- 1
v <- vector("list",N)

for (b in c(1:N)){v[[b]]<-vector("list",length(strt))}
for (b in c(1:N)){for (t in c(1:N)){v[[b]][[t]] <- matrix(NA,nrow=a,ncol=length(strt[[2]]))}}
for (d in c(1:a)){
  for (j in c(1:N)){
    for (k in c(1:length(strt))){
       if (k==1) {if (j < N ){v[[k]][[j]][d,] <- strt[[2]]}else{v[[k]][[j]][d,] <- rep(strt[[1]][m],length(strt[[2]]));m <- m+1}}
       if (k==2) {if (j < N ){v[[k]][[j]][d,] <- strt[[2]]}else{v[[k]][[j]][d,] <- rep(strt[[1]][m],length(strt[[2]]));m <- m+1}}
    }
  }
}

for (i in c(1:N)){
  for (k in c(1:length(strt[[2]]))){
    for (l in c(1:a)){
      if (is.na(v[[i]][[2]][l,k])){v[[i]][[1]][l,k] <- NA}
    }
  }
}

for (b in c(1:N)){for (i in c(1:length(strt))){v[[b]][[i]] <- na.omit(v[[b]][[i]])}}

#DATA TO BE USED (Default is geospiza)
##Geospiza data
#data(geospiza)
#name <- c("geo")
#tree <- geospiza$geospiza.tree
#a.trait <- geospiza$geospiza.data

##Which column of data do you want?
#kk <- 1

##Aquilegia data
name <- c("aqui")
tree <- read.tree("Aquilegia.new.tre")
a.trait <- read.delim("Aquilegia.traits",row.names=1)

##Which column of data do you want?
kk <- 2

##Monocot Data
#name <- c("mono")
#tree <- read.tree("BJO.Monocot.tre")
#a.trait <- read.delim("BJO.monocot_GS",row.names=1)

##Which column of data do you want?
#kk <- 2

##Mammal Data
#name <- c("mam")
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

dv <- treedata(tree,a.trait[,kk],sort=T)
tree <- dv$phy
trait <- data.frame(taxa=rownames(dv$data),as.vector(dv$data))

b.time <-proc.time()
jobs <- lapply(v, function(x) parallel(f(x),silent=TRUE))
results <- collect(jobs,wait=TRUE)
total.time <- as.numeric(proc.time()[3]-b.time[3])/(60)

l <- vector("list",length(well))
for (i in c(1:length(well))){l[[i]] <- vector("list",2)}
for (i in c(1:length(well))){l[[i]][[2]] <- vector("list",length(v))}
names(l) <- well
for (i in c(1:length(well))){l[[i]][[1]] <- matrix(NA,ncol=9,nrow=length(strt[[1]])*length(strt[[2]]));colnames(l[[i]][[1]]) <- c("sa.I","sig.I","lb","ub","L","SA","SIG","conv","time")}

oli <- rep(NA,N)
sum.oli <- rep(NA,N)
for (i in c(1:N)){oli[i] <- length(results[[i]][[1]][[1]][,1])}
for (b in c(1:N)){sum.oli[b] <- sum(oli[b:1])}

for (j in c(1:N)){
  for (i in c(1:z)){
    for (d in c(1:oli[j])){
      if (j == 1 ){l[[i]][[1]][d,] <- as.numeric(results[[j]][[i]][[1]][d,])}else{
                   l[[i]][[1]][d+sum.oli[j-1],] <- as.numeric(results[[j]][[i]][[1]][d,])}
    }
  }
}

for (j in c(1:N)){
  for (i in c(1:z)){
    l[[i]][[2]][[j]] <- results[[j]][[i]][[2]]
  }
}
for (i in c(1:z)){l[[i]][[1]] <- l[[i]][[1]][order(l[[i]][[1]][,1]) , ]}

ll <- vector("list",length(well))
names(ll) <- well
for (j in c(1:length(well))){for (i in c(1:N)){ll[[j]] <- c(ll[[j]],l[[j]][[2]][[i]])}}

for (i in c(1:length(well))){a <- order(as.numeric(names(ll[[i]])));ll[[i]] <- ll[[i]][a]}

for (i in c(1:z)){l[[i]][[2]] <- ll[[i]]}

ll <- NULL

#Saves the results in the directory
save.image(paste(getwd(),"/",name,"ouchctspath.RData",sep=""))
