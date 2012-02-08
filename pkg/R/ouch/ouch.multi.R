#Make sure to have ouchbackup.RData in the correct directory, 
require(ouch)
require(optimx)
require(geiger)
require(multicore)
load("ouchbackup.RData")

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

#M <- detectCores()                      #Use if you want all cores used (Total # of cores)
M <- 2                                   #Use if you want to choose the # of cores to use
it <- 4  	                         #Number of iterations (at least 2)
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

for (b in c(1:N)){v[[b]] <- na.omit(v[[b]])}

#DATA TO BE USED (Default is geospiza)

##Geospiza data
data(geospiza)
sapply(geospiza,class)
name <- c("geo")
nc <- with(geospiza,name.check(geospiza.tree,geospiza.data))
tree <- with(geospiza,drop.tip(geospiza.tree,nc$Tree.not.data))
trait <- data.frame(taxa=rownames(geospiza$geospiza.data),geospiza$geospiza.data)

##Aquilegia data
#source("ouch.wrap.x.R")
#name <- c("aqui")
#tree<-read.tree("Aquilegia.new.tre")
#trait<-read.delim("Aquilegia.traits")

##Monocot Data
#source("ouch.wrap.x.R")
#name <- c("mono")
#tree<-read.tree("BJO.Monocot.tre")
#trait<-read.delim("BJO.monocot_GS")

##Mammal Data
#name <- c("mam")
#tree<-read.nexus("mammalChar1.nex")
#tree$edge.length[tree$edge.length<1e-5]=1e-5
#trait <- read.csv("mammallogChar1.csv")
#trait1 <- read.csv("mammallogChar1.csv", row.names=1)
#nc <- name.check(tree,trait1)
#tree <- drop.tip(tree,nc$Tree.not.data)

b.time <-proc.time()
jobs <- lapply(v, function(x) parallel(f(x),silent=TRUE))
results <- collect(jobs,wait=TRUE)
total.time <- as.numeric(proc.time()[3]-b.time[3])/(60)

l <- vector("list",length(well))
names(l) <- well
for (i in c(1:length(well))){l[[i]] <- matrix(NA,ncol=9,nrow=N);colnames(l[[i]]) <- c("I","lb","ub","T","L","S","A","conv","time")}

for (i in c(1:N)){for (j in c(1:length(well))){l[[j]][i,] <- as.numeric(results[[i]][[j]][1,])}}

#Saves the results in the directory
save.image(paste(getwd(),"/",name,"oucherror.RData",sep=""))
