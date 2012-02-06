require(ouch)

#source("ouch.wrap.new.R")

source("ouch.wrap.new.x.R")

rm(hansen)

setClass(
         'hansentree',
         contains='ouchtree',
         representation=representation(
           call='call',
           nchar='integer',
           optim.diagn='list',
           hessian='matrix',
           data='list',
           regimes='list',
           beta='list',
           theta='list',
           sigma='numeric',
           sqrt.alpha='numeric',
           loglik='numeric'
           )
         )

hansen <- function (data, tree, regimes, sqrt.alpha, sigma,
                    fit = TRUE,
                    method = c("Nelder-Mead","subplex","BFGS","L-BFGS-B"),
                    hessian = FALSE,
                    ...) {

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

  if (fit) { ## maximize the likelihood

    method <- match.arg(method)

    if (method=='subplex') {
      opt <- subplex(
                     par=c(sqrt.alpha,sigma),
                     fn = function (par) {
                       ou.lik.fn(
                                 tree=tree,
                                 alpha=sym.par(par[seq(nalpha)]),
                                 sigma=sym.par(par[nalpha+seq(nsigma)]),
                                 beta=beta,
                                 dat=dat
                                 )$deviance
                     },
                     hessian=hessian,
                     control=list(...)
                     )
      if (opt$convergence!=0) {
        message("unsuccessful convergence, code ",opt$convergence,", see documentation for `subplex'")
        if (!is.null(opt$message))
          message("`subplex' message:",opt$message)
        warning("unsuccessful convergence")
      }
    } else {

wellt <- c("spg", "Rcgmin", "Rvmmin", "bobyqa","L-BFGS-B",1,1,1,1,1,1,1)
well <- c("spg", "Rcgmin", "Rvmmin", "bobyqa","L-BFGS-B","nlminb","ucminf","Nelder-Mead","nlm","CG","BFGS","newuoa")
out <- vector("list", length(well))
for (i in c(1:length(out))){out[[i]] <- vector("list", 2)}
names(out)<-well

  for (i in c(1:length(out))){
    if (well[i]==wellt[i]){
        begin.time <-proc.time()
	out[[i]][[1]] <- tryCatch(optimx(
                      fn = function (par) {
                     ou.lik.fn(
                               tree=tree,
                               alpha=sym.par(par[seq(nalpha)]),
                               sigma=sym.par(par[nalpha+seq(nsigma)]),
                               beta=beta,
                               dat=dat
                               )$deviance
                   }, 
                      p=c(sqrt.alpha,sigma),
                      lower = rep(0.001,(nalpha+nsigma)), upper = rep(1e50,(nalpha+nsigma)),
                      method=well[i]
                      ), error=function(x){
                                           x$par[seq(nalpha)] <- NA
                                           x$par[nalpha+seq(nsigma)] <- NA
                                           return(x)})
        out[[i]][[2]] <- as.numeric(proc.time()[3]-begin.time[3])/(60)}else{
        begin.time <-proc.time()
	out[[i]][[1]] <- tryCatch(optimx(
                      fn = function (par) {
                     ou.lik.fn(
                               tree=tree,
                               alpha=sym.par(par[seq(nalpha)]),
                               sigma=sym.par(par[nalpha+seq(nsigma)]),
                               beta=beta,
                               dat=dat
                               )$deviance
                   }, 
                      p=c(sqrt.alpha,sigma),
                      method=well[i]
                      ), error=function(x){
                                           x$par[seq(nalpha)] <- NA
                                           x$par[nalpha+seq(nsigma)] <- NA
                                           return(x)})
        out[[i]][[2]] <- as.numeric(proc.time()[3]-begin.time[3])/(60)}
}

      	opt <- optim(
                   par=c(sqrt.alpha,sigma),
                   fn = function (par) {
                     ou.lik.fn(
                               tree=tree,
                               alpha=sym.par(par[seq(nalpha)]),
                               sigma=sym.par(par[nalpha+seq(nsigma)]),
                               beta=beta,
                               dat=dat
                               )$deviance
                   },
                   gr=NULL,
                   hessian=hessian,
                   method=method,
                   control=list(...)
                   )

      if (opt$convergence!=0) {
        message("unsuccessful convergence, code ",opt$convergence,", see documentation for `optim'")
        if (!is.null(opt$message))
          message("`optim' message:",opt$message)
        warning("unsuccessful convergence")
      }
    }
  }
well <- c("spg", "Rcgmin", "Rvmmin", "bobyqa","L-BFGS-B","nlminb","ucminf","Nelder-Mead","nlm","CG","BFGS","newuoa")
  obj.list <- vector("list", length(out))
  for (i in c(1:length(obj.list))){obj.list[[i]] <- vector("list", 2)}
  names(obj.list) <- well
  for (i in c(1:length(obj.list))){

    sqrt.alpha <- unlist(out[[i]][[1]]$par)[seq(nalpha)]
#   sqrt.alpha <- opt$par[seq(nalpha)]
    sigma <- unlist(out[[i]][[1]]$par)[nalpha+seq(nsigma)]
#   sigma <- opt$par[nalpha+seq(nsigma)]
    time <- unlist(out[[i]][[2]])
    optim.diagn <- list(convergence=opt$convergence,message=opt$message)

  if (hessian) {
    hs <- opt$hessian
#     se <- sqrt(diag(solve(0.5*hs)))
#     se.alpha <- se[seq(nalpha)]
#     se.sigma <- se[nalpha+seq(nsigma)]
  } else {
    hs <- matrix(NA,0,0)
#     se.alpha <- rep(NA,nalpha)
#     se.sigma <- rep(NA,nalpha)
  }

  sol <- ou.lik.fn(
                   tree=tree,
                   alpha=sym.par(sqrt.alpha),
                   sigma=sym.par(sigma),
                   beta=beta,
                   dat=dat
                   )
  theta.x <- sol$coeff
  reg <- sets.of.regimes(tree,regimes)
  theta <- vector('list',nchar)
  names(theta) <- names(data)
  count <- 1
  for (n in seq_len(nchar)) {
    theta[[n]] <- theta.x[seq(from=count,length=length(reg[[n]]),by=1)]
    names(theta[[n]]) <- as.character(reg[[n]])
    count <- count+length(reg[[n]])
  }
  
  obj.list[[i]][[1]] <- new(
      'hansentree',
      as(tree,'ouchtree'),
      call=match.call(),
      nchar=nchar,
      optim.diagn=optim.diagn,
      hessian=hs,
      data=as.list(data),
      regimes=as.list(regimes),
      beta=beta,
      theta=theta,
      sigma=sigma,
      sqrt.alpha=sqrt.alpha,
      loglik=-0.5*sol$deviance
      )
  obj.list[[i]][[2]] <- time
  }
obj.list
}

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

wellt <- c("spg", "Rcgmin", "Rvmmin", "bobyqa","L-BFGS-B",1,1,1,1,1,1,1)
well <- c("spg", "Rcgmin", "Rvmmin", "bobyqa","L-BFGS-B","nlminb","ucminf","Nelder-Mead","nlm","CG","BFGS","newuoa")

#detectCores <- function(all.tests = FALSE)
#{
#    systems <-
#        list(darwin  = "/usr/sbin/sysctl -n hw.ncpu 2>/dev/null",
#             freebsd = "/sbin/sysctl -n hw.ncpu 2>/dev/null",
#             linux   = "grep processor /proc/cpuinfo 2>/dev/null | wc -l",
#             irix    = c("hinv | grep Processors | sed 's: .*::'",
#                         "hinv | grep '^Processor '| wc -l"),
#             solaris = "/usr/sbin/psrinfo -v | grep 'Status of.*processor' | wc -l"
#             solaris = "/usr/sbin/psrinfo -p")
#    for (i in seq(systems))
#        if(all.tests ||
#           length(grep(paste("^", names(systems)[i], sep=''), R.version$os)))
#            for (cmd in systems[i]) {
#                a <- gsub("^ +","", system(cmd, TRUE)[1]) 
#                if (length(grep("^[1-9]", a))) return(as.integer(a))
#            }
#    NA_integer_
#}
#Input the following, the function will automatically find the number of cores but if you want to use a specific number of cores enter yourself.
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
require(optimx)
require(geiger)
tree<-read.nexus("mammalChar1.nex")
tree$edge.length[tree$edge.length<1e-5]=1e-5
trait <- read.csv("mammallogChar1.csv")
trait1 <- read.csv("mammallogChar1.csv", row.names=1)
nc <- name.check(tree,trait1)
tree <- drop.tip(tree,nc$Tree.not.data)
#ouch.wrap(tree,trait,model=c("ou1"))

begin.time <-proc.time()
l <- ouch.wrap.new(tree,trait,model=c("ou1"))
total.time <- as.numeric(proc.time()[3]-begin.time[3])/(60)
#what?

f <- function(x)
{ 
yy <- x

m <- length(yy)
jj <- yy

optimx <- vector("list", length(well))
names(optimx)<-well

for (i in c(1:length(optimx))){
      optimx[[i]]<-vector("list",m)
      names(optimx[[i]]) <- jj
}

lt <- vector("list",length(optimx))
names(lt) <- well
time <- vector("list",length(optimx))
for (i in c(1:length(optimx))){time[[i]] <- rep(NA,m)}

k<-matrix(jj,1,m)
k<-rbind(k,apply(k,2,function(x){ouch.wrap.x(tree,trait,model=c("ou1"),x)}))

for (i in c(1:length(optimx))){for (j in c(1:m)){optimx[[i]][[j]] <- k[2,][[j]][[i]][[1]];time[[i]][j] <- k[2,][[j]][[i]][[2]]}}

for(j in c(1:length(optimx))){
  if (well[j]==wellt[j]){
    lt[[j]]<-as.data.frame(t(rbind(unlist(k[1,]),as.data.frame(matrix(unlist(lapply(optimx[[j]],function(x) return(c(0.001,1e50,as.numeric(x@theta),x@loglik,x@sigma,x@sqrt.alpha)))),6,ncol(k))),time[[j]])))
    colnames(lt[[j]])<-c("I","lb","ub","T","L","S","A","time")}else{
    lt[[j]]<-as.data.frame(t(rbind(unlist(k[1,]),as.data.frame(matrix(unlist(lapply(optimx[[j]],function(x) return(c(NA,NA,as.numeric(x@theta),x@loglik,x@sigma,x@sqrt.alpha)))),6,ncol(k))),time[[j]])))
    colnames(lt[[j]])<-c("I","lb","ub","T","L","S","A","time")}
}

return(lt)
}

#Required to run this part.
require(ouch)
require(optimx)
require(geiger)
load("ouchbackup.RData")

#Default inputs, if you want more iterations, change it, then rerun this piece of code.
M <- 22
it <- 50  	                         #Number of iterations
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

require(multicore)
b.time <-proc.time()
jobs <- lapply(v, function(x) parallel(f(x),silent=TRUE))
results <- collect(jobs,wait=TRUE)
total.time <- as.numeric(proc.time()[3]-b.time[3])/(60)
#what?

#save.image("/home/michels/repository/phylooptim/pkg/R/ouch/geooucherror.RData")
#save.image("/home/michels/repository/phylooptim/pkg/R/ouch/aquioucherror.RData")
#save.image("/Users/michels/phylooptim/pkg/R/ouch/monooucherror.RData")
save.image("/Users/michels/phylooptim/pkg/R/ouch/mamoucherror.RData")
#save.image("/home/michels/repository/phylooptim/pkg/R/ouch/monooucherror.RData")
