rm(list=ls())

require(ouch)
require(geiger)

well <- c("spg", "Rcgmin", "Rvmmin", "bobyqa","L-BFGS-B","nlminb","ucminf","Nelder-Mead","nlm","CG","BFGS","newuoa")
#actual order of well
#well[c(7,6,8,12,1,2,4,9,3,5,10,11)]
#load("/home/michels/Hallowed/repository/phylooptim/pkg/R/ouch/geoouchctspath2.RData")

#load("/home/michels/Hallowed/repository/phylooptim/pkg/R/ouch/aquiouchctspath.RData")
load("aquiouchctspath.RData")

##Geospiza data
#data(geospiza)
#name <- c("geo")
#tree <- geospiza$geospiza.tree
#a.trait <- geospiza$geospiza.data

##Which column of data do you want?
#kk <- 1

##Aquilegia Data
#tree <- read.tree("/home/michels/Hallowed/repository/phylooptim/pkg/R/ouch/Aquilegia.new.tre")
tree <- read.tree("Aquilegia.new.tre")
#a.trait <- read.delim("/home/michels/Hallowed/repository/phylooptim/pkg/R/ouch/Aquilegia.traits",row.names=1)
a.trait <- read.delim("Aquilegia.traits",row.names=1)
##Which column of data do you want?
kk <- 2

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

ot <- ape2ouch(tree)
otd <- as(ot,"data.frame")
trait <- data.frame(phenoTrait=trait[,2], labels=trait[,1])
otd <- merge(otd,trait, by="labels", all=TRUE)
rownames(otd) <- otd$nodes
ot <- with(otd,ouchtree(nodes=nodes,ancestors=ancestors,times=times,labels=labels))
	
otd$regimes <- as.factor("global")

tree=ot
data=otd["phenoTrait"]
regimes=otd["regimes"]
sqrt.alpha=1
sigma=1

#maxit=10000
#fit = TRUE
#method = c("Nelder-Mead","subplex","BFGS","L-BFGS-B")
#hessian = FALSE

source("glssoln.R")
dyn.load("weight-matrix.so")
dyn.load("covar-matrix.so")
## note that, on input, alpha and sigma are full symmetric matrices
ou.lik.fn <- function (tree, alpha, sigma, beta, dat) {
  n <- length(dat)
  ev <- eigen(alpha,symmetric=TRUE)
  w <- .Call("ouch_weights",object=tree,lambda=ev$values,S=ev$vectors,beta=beta)
  v <- .Call("ouch_covar",object=tree,lambda=ev$values,S=ev$vectors,sigma.sq=sigma)
  gsol <- try(
              glssoln(w,dat,v),
              silent=FALSE
              )
  if (inherits(gsol,'try-error')) { # return Inf deviance (so that optimizer can keep trying)
    e <- rep(NA,n)
    theta <- rep(NA,ncol(w))
    dev <- Inf
  } else {                              # return finite deviance
    e <- gsol$residuals
    theta <- gsol$coeff
    q <- e%*%solve(v,e)
    det.v <- determinant(v,logarithm=TRUE)
    if (det.v$sign!=1)
      stop("ou.lik.fn error: non-positive determinant",call.=FALSE)
    dev <- n*log(2*pi)+as.numeric(det.v$modulus)+q[1,1]
  }
  list(
       deviance=dev,
       coeff=theta,
       weight=w,
       vcov=v,
       resids=e
       )
}

sym.par <- function (x) {
  nchar <- floor(sqrt(2*length(x)))
  if (nchar*(nchar+1)!=2*length(x)) {
    stop("a symmetric matrix is parameterized by a triangular number of parameters",call.=FALSE)
  }
  y <- matrix(0,nchar,nchar)
  y[lower.tri(y,diag=TRUE)] <- x
  y%*%t(y)
}

sym.unpar <- function (x) {
  y <- t(chol(x))
  y[lower.tri(y,diag=TRUE)]
}

sets.of.regimes <- function (object, regimes) {
  lapply(regimes,function(x)sort(unique(x)))
}

regime.spec <- function (object, regimes) {
  nterm <- object@nterm
  nchar <- length(regimes)
  reg <- sets.of.regimes(object,regimes)
  nreg <- sapply(reg,length)
  beta <- vector(mode='list',length=nterm)
  for (i in seq_len(nterm)) {
    p <- object@lineages[[object@term[i]]]
    np <- length(p)
    beta[[i]] <- vector(mode='list',length=nchar)
    for (n in seq_len(nchar)) {
      beta[[i]][[n]] <- matrix(data=NA,nrow=np,ncol=nreg[n])
      for (ell in seq_len(nreg[n])) {
        beta[[i]][[n]][,ell] <- ifelse(regimes[[n]][p]==reg[[n]][ell],1,0)
      }
    }
  }
  beta
}

## solve the matrix equation
##   A . X + X . A = B
## for X, where we have assumed A = A'.
sym.solve <- function (a, b) {
  n <- nrow(a)
  d <- array(data=0,dim=c(n,n,n,n))
  for (k in seq_len(n)) {
    d[k,,k,] <- d[k,,k,] + a
    d[,k,,k] <- d[,k,,k] + a
  }
  dim(b) <- n*n
  dim(d) <- c(n*n,n*n)
  x <- solve(d,b)
  dim(x) <- c(n,n)
  x
}

hansen.deviate <- function (n = 1, object) {
  ev <- eigen(sym.par(object@sqrt.alpha),symmetric=TRUE)
  w <- .Call(ouch_weights,object=object,lambda=ev$values,S=ev$vectors,beta=object@beta)
  v <- .Call(ouch_covar,object=object,lambda=ev$values,S=ev$vectors,sigma.sq=sym.par(object@sigma))
  X <- array(
             data=NA,
             dim=c(object@nnodes,object@nchar,n),
             dimnames=list(
               object@nodes,
               names(object@data),
               paste('rep',seq(n),sep='.')
               )
             )

  theta <- do.call(c,object@theta)

  X[object@term,,] <- array(
                            data=rmvnorm(
                              n=n,
                              mean=as.numeric(w%*%theta),
                              var=v
                              ),
                            dim=c(object@nterm,object@nchar,n)
                            )
  apply(X,3,as.data.frame)
}

setMethod(
          'simulate',
          'hansentree',
          function (object, nsim = 1, seed = NULL, ...) {
            if (!is.null(seed)) {
              if (!exists('.Random.seed',envir=.GlobalEnv)) runif(1)
              save.seed <- get('.Random.seed',envir=.GlobalEnv)
              set.seed(seed)
            }
            X <- hansen.deviate(n=nsim,object)
            if (!is.null(seed)) {
              assign('.Random.seed',save.seed,envir=.GlobalEnv)
            }
            X
          }
          )

setMethod(
          'update',
          'hansentree',
          function (object, data, regimes, sqrt.alpha, sigma, ...) {
            if (missing(sqrt.alpha)) sqrt.alpha <- object@sqrt.alpha
            if (missing(sigma)) sigma <- object@sigma
            hansen(
                   data=data,
                   tree=object,
                   regimes=regimes,
                   sqrt.alpha=sqrt.alpha,
                   sigma=sigma,
                   ...
                   )
          }
          )


setMethod(
          "bootstrap",
          "hansentree",
          function (object, nboot = 200, seed = NULL, ...) {
            simdata <- simulate(object,nsim=nboot,seed=seed)
            results <- vector(mode='list',length=nboot)
            toshow <- c("alpha","sigma.squared","optima","loglik","aic","aic.c","sic","dof")
            for (b in seq_len(nboot)) {
              results[[b]] <- summary(update(object,data=simdata[[b]],...))
            }
            as.data.frame(t(sapply(results,function(x)unlist(x[toshow]))))
          }
          )
		
  if (!is(tree,'ouchtree'))
    stop(sQuote("tree")," must be an object of class ",sQuote("ouchtree"))

  if (missing(data)) {
    if (is(tree,"hansentree")) {
      data <- tree@data
    } else {
      stop(sQuote("data")," must be specified")
    }
  }
  if (is.data.frame(data)) {
    nm <- rownames(data)
    data <- lapply(as.list(data),function(x){names(x)<-nm;x})
  }
  if (is.numeric(data)) {
    nm <- deparse(substitute(data))[1]
    data <- list(data)
    names(data) <- nm
  }
  if (is.list(data)) {
    if (
        any(sapply(data,class)!='numeric') ||
        any(sapply(data,length)!=tree@nnodes)
        )
      stop(sQuote("data")," vector(s) must be numeric, with one entry per node of the tree")
    if (any(sapply(data,function(x)(is.null(names(x)))||(!setequal(names(x),tree@nodes)))))
      stop(sQuote("data"), " vector names (or data-frame row names) must match node names of ", sQuote("tree"))
    for (xx in data) {
      no.dats <- which(is.na(xx[tree@nodes[tree@term]]))
      if (length(no.dats)>0)
        stop("missing data on terminal node(s): ",paste(tree@nodes[tree@term[no.dats]],collapse=','))
    }
  } else
  stop(sQuote("data")," must be either a single numeric data set or a list of numeric data sets")

  nchar <- length(data)
  if (is.null(names(data))) names(data) <- paste('char',seq_len(nchar),sep='')
  
  if (any(sapply(data,function(x)(is.null(names(x)))||(!setequal(names(x),tree@nodes)))))
    stop("each data set must have names corresponding to the node names")
  data <- lapply(data,function(x)x[tree@nodes])
  dat <- do.call(c,lapply(data,function(y)y[tree@term]))
  
  nsymargs <- nchar*(nchar+1)/2
  nalpha <- length(sqrt.alpha)
  nsigma <- length(sigma)

  if (nalpha!=nsymargs)
    stop("the length of ",sQuote("sqrt.alpha")," must be a triangular number")

  if (nsigma!=nsymargs)
    stop("the length of ",sQuote("sigma")," must be a triangular number")

  if (missing(regimes)) {
    if (is(tree,"hansentree")) {
      regimes <- tree@regimes
      beta <- tree@beta
    } else {
      stop(sQuote("regimes")," must be specified")
    }
  }
  if (is.data.frame(regimes)) {
    nm <- rownames(regimes)
    regimes <- lapply(as.list(regimes),function(x){names(x)<-nm;x})
  }
  if (is.list(regimes)) {
    if (any(sapply(regimes,length)!=tree@nnodes))
      stop("each element in ",sQuote("regimes")," must be a vector with one entry per node of the tree")
  } else {
    if (length(regimes)!=tree@nnodes)
      stop("there must be one entry in ",sQuote("regimes")," per node of the tree")
    nm <- deparse(substitute(regimes))[1]
    regimes <- list(regimes)
    names(regimes) <- nm
  }
  
  if (any(!sapply(regimes,is.factor)))
    stop(sQuote("regimes")," must be of class ",sQuote("factor")," or a list of ",sQuote("factor")," objects")

  if (length(regimes)==1)
    regimes <- rep(regimes,nchar)

  if (length(regimes) != nchar)
    stop("you must supply a regime-specification vector for each character")

  if (any(sapply(regimes,function(x)(is.null(names(x)))||(!setequal(names(x),tree@nodes)))))
    stop("each regime specification must have names corresponding to the node names")
  regimes <- lapply(regimes,function(x)x[tree@nodes])
  beta <- regime.spec(tree,regimes)

  optim.diagn <- vector(mode='list',length=0)

##For clean contour plot
##find min & max of sqrt.alpha & sigma
#min.sa <- rep(NA,length(well));max.sa <- rep(NA,length(well));min.sig <- rep(NA,length(well));max.sig <- rep(NA,length(well))

#for (i in c(1:length(well))){min.sa[i] <- min(sa[[i]]);max.sa[i] <- max(sa[[i]]);min.sig[i] <- min(sig[[i]]);max.sig[i] <- max(sig[[i]])}

#sqrt.alpha
#x1range <- seq(min(na.omit(min.sa)),max(na.omit(max.sa)),length=100)

#sigma
#x2range <- seq(min(na.omit(min.sig)),floor(max(na.omit(max.sig)))+1,length=100)

#The response
#  res<-t(sapply(x1range,function(x,y=x2range) sapply(y,function(z) ou.lik.fn(
#                               tree=tree,
#                               alpha=x^2,
#                               sigma=z^2,
#                               beta=beta,
#                               dat=dat
#                               )$deviance)))

#For gray contour plot
#Exclude nlm & BFGS since predicts negative values.
#min.sa.list <- list()
#max.sa.list <- list()
#min.sig.list <- list()
#max.sig.list <- list()

#for (j in c(1:length(well))){
#  min.sa.list[[j]] <- rep(NA,length(c(1:length(l[[j]][[2]]))[seq(1,length(l[[j]][[2]]),by=1)]))
#  max.sa.list[[j]] <- rep(NA,length(c(1:length(l[[j]][[2]]))[seq(1,length(l[[j]][[2]]),by=1)]))
#  min.sig.list[[j]] <- rep(NA,length(c(1:length(l[[j]][[2]]))[seq(1,length(l[[j]][[2]]),by=1)]))
#  max.sig.list[[j]] <- rep(NA,length(c(1:length(l[[j]][[2]]))[seq(1,length(l[[j]][[2]]),by=1)]))
#}

#for (j in c(1:length(well))){
#   for (i in c(1:length(l[[j]][[2]]))[seq(1,length(l[[j]][[2]]),by=1)]){
#     min.sa.list[[j]][i] <- min(na.omit(l[[j]][[2]][[i]][[1]]$sqrt.alpha))
#     max.sa.list[[j]][i] <- max(na.omit(l[[j]][[2]][[i]][[1]]$sqrt.alpha))
#     min.sig.list[[j]][i] <- min(na.omit(l[[j]][[2]][[i]][[1]]$sigma))
#     max.sig.list[[j]][i] <- max(na.omit(l[[j]][[2]][[i]][[1]]$sigma))
#   }
# }

#min.sa <- rep(NA,length(well))
#max.sa <- rep(NA,length(well))
#min.sig <- rep(NA,length(well))
#max.sig <- rep(NA,length(well))

#for (j in c(1:length(well))){
#  min.sa[j] <- min(min.sa.list[[j]])
#  max.sa[j] <- max(max.sa.list[[j]])
#  min.sig[j] <- min(min.sig.list[[j]])
#  max.sig[j] <- max(max.sig.list[[j]])
#}


for (i in c(1:length(well))){
  for (j in c(1:length(l[[i]][[2]]))){
    l[[i]][[2]][[j]][[1]]$sqrt.alpha <- abs(l[[i]][[2]][[j]][[1]]$sqrt.alpha)
    l[[i]][[2]][[j]][[1]]$sigma <- abs(l[[i]][[2]][[j]][[1]]$sigma)
  }
}

#c(2,5,16,17,18,20,22,23,24,25)
#for (i in c(1:length(well))){
#  i <- 7
#  print(abs(l[[i]][[1]][,c(1,2,6,7)]))
#}
#for (i in c(1:length(well))){
#  print(abs(l[[i]][[1]][c(2,18,22,25),c(1,2,6,7)]))
#}

#sqrt.alpha
#
#x1range <- seq(min(min.sa),max(max.sa),length=100)
x1range1 <- seq(0,2,length=100)
x1range <- seq(2.05,8.5,length=100)
x1range <- c(x1range1,x1range)
#sigma
#x2range <- seq(min(na.omit(min.sig)),floor(max(na.omit(max.sig)))+1,length=100)
x2range1 <- seq(0,2,length=100)
x2range <- seq(2.05,8.5,length=100)
x2range <- c(x2range1,x2range)

#The response
  res<-t(sapply(x1range,function(x,y=x2range) sapply(y,function(z) ou.lik.fn(
                               tree=tree,
                               alpha=x^2,
                               sigma=z^2,
                               beta=beta,
                               dat=dat
                               )$deviance/-2)))

#Default use filled.contour to get side legend

#Filled contour without legend.
source("Filled.contour2.R")

#Needed for a filled contour plot, in grayscale
      numcol<-50
#     nlevels=50
#     levels=pretty(range(res,na.rm=TRUE),nlevels,finite=TRUE)
#     k<-pretty(trunc(range(res,na.rm=T)),nlevels)
      wr.pal<-colorRampPalette(c("white","grey"))
      wr <- wr.pal(numcol)
#     wr<-wr.pal(round(numcol*sum(k>=0)/length(k)))
      bw.pal<-colorRampPalette(c("grey","black"))
#     bw<-bw.pal(round(numcol*sum(k<=0)/length(k)))
      bw <- bw.pal(numcol)
      cols<-c(wr,bw)
#     rcols<-cols[round((1:(length(levels)))*length(cols)/length(levels))]

#x.ticks <- c(sa[[1]][1],sa[[12]][length(sa[[12]])],sa[[1]][length(sa[[1]])],sa[[10]][length(sa[[10]])],sa[[5]][length(sa[[5]])],sa[[11]][length(sa[[11]])],sa[[9]][length(sa[[9]])])
#x.labels <- c(round(sa[[1]][1],2),round(sa[[12]][length(sa[[12]])],2),round(sa[[1]][length(sa[[1]])],2),round(sa[[10]][length(sa[[10]])],2),round(sa[[5]][length(sa[[5]])],2),round(sa[[11]][length(sa[[11]])],2),round(sa[[9]][length(sa[[9]])],2))
#y.ticks <- c(c(sig[[1]][1],sig[[12]][length(sig[[12]])],sig[[1]][length(sig[[1]])],sig[[10]][length(sig[[10]])],sig[[5]][length(sig[[5]])]),sig[[11]][length(sig[[11]])],sig[[9]][length(sig[[9]])],-4)
#y.labels <- c("","","","","","","","")
x.labels <- c(0:8)
y.labels <- c(0:8)

#lev <- c(seq(-30-33*10,-40,by=10),seq(-30,-22,length.out=17))
lev <- c(seq(-30-83*.5,-30.5,by=.5),seq(-30,-22,length.out=17))
#filled.contour2(x1range,x2range,res,levels=lev,col=cols,xlab="sqrt.alpha",ylab="sigma",main="Aquilegia log likelihood surface")
#png(file="ouchcontourtotalgrey.png")
#filled.contour2(x1range,x2range,res,levels=levels1,col=rcols1,xlab="sqrt.alpha",ylab="sigma",main="Aquilegia log likelihood surface")
#filled.contour(x1range,x2range,res,levels=levels,col=rcols)
for (j in c(1:length(well))){
png(file=paste("ouchcontour",well[j],"black.png",sep=""))
filled.contour2(x1range,x2range,res,levels=lev,col=cols,axes=FALSE,
      		plot.axes={
                  axis(side = 1, at = x.labels,labels=x.labels)
                  axis(side = 2, at = y.labels,labels=y.labels)
            plot.title=title(xlab="sqrt.alpha", ylab="sigma",main=paste(well[j]),cex.main=2)
                      }
               )
#First start value
for (i in c(1:length(l[[j]][[2]]))[c(2,18,22,25)]){
  points(l[[j]][[2]][[i]][[1]]$sqrt.alpha[1],l[[j]][[2]][[i]][[1]]$sigma[1],col="black",pch=1,cex=1.5,lwd=2)
}

#Path
for (i in c(1:length(l[[j]][[2]]))[c(2,18,22,25)]){
  points(na.omit(l[[j]][[2]][[i]][[1]]$sqrt.alpha),na.omit(l[[j]][[2]][[i]][[1]]$sigma),col="black",type="l",lwd=2)
}

#End start value
for (i in c(1:length(l[[j]][[2]]))[c(2,18,22,25)]){
    points(l[[j]][[2]][[i]][[1]]$sqrt.alpha[length(l[[j]][[2]][[i]][[1]]$sqrt.alpha)],l[[j]][[2]][[i]][[1]]$sigma[length(l[[j]][[2]][[i]][[1]]$sigma)],col="white",pch=20,cex=1.25,lwd=2)
  points(l[[j]][[2]][[i]][[1]]$sqrt.alpha[length(l[[j]][[2]][[i]][[1]]$sqrt.alpha)],l[[j]][[2]][[i]][[1]]$sigma[length(l[[j]][[2]][[i]][[1]]$sigma)],col="black",pch=1,cex=1.5,lwd=2)
}

dev.off()
}

#mtext(at=2.2, line=c(-1.1,-9,-10,-11,-12,-13,-19.5,-25), text=c(round(sig[[1]][1],2),round(sig[[12]][length(sig[[12]])],2),round(sig[[1]][length(sig[[1]])],2),round(sig[[10]][length(sig[[10]])],2),round(sig[[5]][length(sig[[5]])],2),round(sig[[11]][length(sig[[11]])],2),round(sig[[9]][length(sig[[9]])],2),-4),cex = 1)
#lev <- c(seq(-25000,-15,by=5),-13,-11.5,-9.5,-7.5,-6,-4,-2,-0.5,3.5,5,7,9)
#contour(x1range,x2range,res,levels=lev,add=TRUE)

#  arrows(sa[[1]][s], sig[[1]][s], sa[[1]][s+1],sig[[1]][s+1],length=.1,col=1,type = "l",lwd=2,lty=1)
#  arrows(sa[[2]][s], sig[[2]][s], sa[[2]][s+1],sig[[2]][s+1],length=.1,col="green2",type = "l",lwd=2,lty=2)
#  arrows(na.omit(sa[[3]])[sy], na.omit(sig[[3]])[sy], na.omit(sa[[3]])[sy+1],na.omit(sig[[3]])[sy+1],length=.1,col="yellow",type = "l",lwd=2,lty=3)
#  arrows(sa[[4]][s], sig[[4]][s], sa[[4]][s+1],sig[[4]][s+1],length=.1,col="pink3",type = "l",lwd=2,lty=4)
#  arrows(sa[[5]][s], sig[[5]][s], sa[[5]][s+1],sig[[5]][s+1],length=.1,col="yellow3",type = "l",lwd=2,lty=5)
#  arrows(sa[[6]][s], sig[[6]][s], sa[[6]][s+1],sig[[6]][s+1],length=.1,col="red3",type = "l",lwd=2,lty=6)
#  arrows(sa[[7]][s], sig[[7]][s], sa[[7]][s+1],sig[[7]][s+1],length=.1,col="green4",type = "l",lwd=2,lty=7)
#  arrows(sa[[8]][s], sig[[8]][s], sa[[8]][s+1],sig[[8]][s+1],length=.1,col="blue4",type = "l",lwd=2,lty=8)
#  arrows(sa[[9]][s], sig[[9]][s], sa[[9]][s+1],sig[[9]][s+1],length=.1,col="red4",type = "l",lwd=2,lty=9)
#  arrows(sa[[10]][s], sig[[10]][s], sa[[10]][s+1],sig[[10]][s+1],length=.1,col=6,type = "l",lwd=2,lty=10)
#  arrows(sa[[11]][s], sig[[11]][s], sa[[11]][s+1],sig[[11]][s+1],length=.1,col="blue",type = "l",lwd=2,lty=11)
#  arrows(sa[[12]][s], sig[[12]][s], sa[[12]][s+1],sig[[12]][s+1],length=.1,col="yellow4",type = "l",lwd=2,lty=12)

#Start-points, same for every optimizer
#  points(sa[[1]][1],sig[[1]][1],col="green",pch=19,cex=2)
                        
#End-points                     
#  points(sa[[1]][length(sa[[1]])],sig[[1]][length(sig[[1]])],col=1,pch=1,cex=3,lwd=2)
#  points(sa[[2]][length(sa[[2]])],sig[[2]][length(sig[[2]])],col="green2",pch=2,cex=3,lwd=2)
#  points(na.omit(sa[[3]])[length(na.omit(sa[[3]]))],na.omit(sig[[3]])[length(na.omit(sig[[3]]))],col="yellow",pch=3,cex=3,lwd=2)
#  points(sa[[4]][length(sa[[4]])],sig[[4]][length(sig[[4]])],col="pink3",pch=4,cex=3,lwd=2)
#  points(sa[[5]][length(sa[[5]])],sig[[5]][length(sig[[5]])],col="yellow3",pch=5,cex=3,lwd=2)
#  points(sa[[6]][length(sa[[6]])],sig[[6]][length(sig[[6]])],col="red3",pch=6,cex=3,lwd=2)
#  points(sa[[7]][length(sa[[7]])],sig[[7]][length(sig[[7]])],col="green4",pch=7,cex=3,lwd=2)
#  points(sa[[8]][length(sa[[8]])],sig[[8]][length(sig[[8]])],col="blue4",pch=8,cex=3,lwd=2)
#  points(sa[[9]][length(sa[[9]])],sig[[9]][length(sig[[9]])],col="red4",pch=9,cex=3,lwd=2)
#  points(sa[[10]][length(sa[[10]])],sig[[10]][length(sig[[10]])],col=6,pch=10,cex=3,lwd=2)
#  points(sa[[11]][length(sa[[11]])],sig[[11]][length(sig[[11]])],col="blue",pch=11,cex=3,lwd=2)
#  points(sa[[12]][length(sa[[12]])],sig[[12]][length(sig[[12]])],col="yellow4",pch=12,cex=3,lwd=2)

#True End-Point
#  points(sa[[12]][length(sa[[12]])],sig[[12]][length(sig[[12]])],col="red",pch=19,cex=2)

#Legend
#legend(2.75,-1, c("start value","true end value"),
#     col = c("green","red"),
#     text.col =c(1) ,
#     bty="n",pch=c(19),
#     merge = TRUE)

#legend(2.68,-1.75, well[1:5],
#     col = c(1,"green2","yellow","pink3","yellow3"),
#     text.col =c(1,"green2","yellow","pink3","yellow3") ,
#     lty = c(1:5),bty="n",pch=c(1:5),
#     merge = TRUE)

#legend(5.5,-1, well[6:12],
#     col = c("red3","green4","blue4","red4",6,"blue","yellow4"),
#     text.col =c("red3","green4","blue4","red4",6,"blue","yellow4") ,
#     lty = c(6:12),bty="n",pch=c(6:12),
#     merge = TRUE)

#dev.off()
