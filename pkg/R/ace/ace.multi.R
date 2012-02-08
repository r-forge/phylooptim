#Make sure to have acebackup.RData in the correct directory
require(ouch)
require(optimx)
require(geiger)
require(multicore)
load("acebackup.RData")

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
name <- c("geo")
tree <- geospiza$geospiza.tree
a <- data.frame(geospiza$geospiza.data,T=as.factor(geospiza$geospiza.data[,1]>4.2))

##Which column of data do you want?
jj <- 6

##Aquilegia Data
#name <- c("aqui")
#tree <- read.tree("Aquilegia.new.tre")
#a <- read.delim("Aquilegia.traits",row.names=1)

##Which column of data do you want?
#jj <- 1

##Monocot Data
#name <- c("mono")
#tree <- read.tree("BJO.Monocot.tre")
#a <- read.delim("BJO.monocot_GS",row.names=1)

##Which column of data do you want?
#jj <- 3

##Mammal Data
#name <- c("mam")
#tree<-read.nexus("mammalChar4.nex")
#a <- read.csv("mammalChar4.csv", row.names=1)

##Which column of data do you want?
#jj <- 1

trait <- data.frame(T=a[,jj])
rownames(trait) <- rownames(a)
name <- c("geo")
tree$edge.length[tree$edge.length<1e-5]=1e-5
nc <- name.check(tree,a)
if (nc[1]=="OK"){nc$Tree.not.data <- NULL}
tree <- drop.tip(tree,nc$Tree.not.data)
dv <- treedata(tree,trait,sort=T)

b.time <-proc.time()
jobs <- lapply(v, function(x) parallel(f(x),silent=TRUE))
results <- collect(jobs,wait=TRUE)
total.time <- as.numeric(proc.time()[3]-b.time[3])/(60)

l <- vector("list",length(well))
names(l) <- well
for (i in c(1:length(well))){l[[i]] <- matrix(NA,ncol=7,nrow=N);colnames(l[[i]]) <- c("I","lb","ub","P","L","conv","time")}

for (i in c(1:N)){for (j in c(1:length(well))){l[[j]][i,] <- as.numeric(results[[i]][[j]][1,])}}

#Saves the results in the directory
save.image(paste(getwd(),"/",name,"aceerror.RData",sep=""))
