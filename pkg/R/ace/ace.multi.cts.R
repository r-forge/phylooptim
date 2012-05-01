#Make sure to have acebackupcts.RData in the correct directory
rm(list = ls())
require(ouch)
require(optimx)
require(geiger)
source("ace.mine.cts.R")

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
#kk <- 2

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
tree$edge.length[tree$edge.length<.1]=.1

a.trait <- a.trait[order(rownames(a.trait)),]
d <- data.frame(name=tree$tip.label,num=seq(1,length(tree$tip.label),by=1))
dd <- d[order(d[,1]) , ]
a.trait$new <- dd[,2]
a.trait <- a.trait[order(a.trait$new),]

dv <- treedata(tree,a.trait[,kk],sort=T)

b.time <-proc.time()
l <- ace(as.vector(dv$data),dv$phy,type='continuous',CI=FALSE)
total.time <- as.numeric(proc.time()[3]-b.time[3])/(60)

#Saves the results in the directory
save.image(paste(getwd(),"/",name,"acects.RData",sep=""))
