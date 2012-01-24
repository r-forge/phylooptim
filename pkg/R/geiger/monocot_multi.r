
#setwd("/media/GM/Dropbox/repository/phylooptim/pkg/R/geiger")
#/home/michels/Hallowed/Dropbox/repository/phylooptim/pkg/R/geiger
setwd("/Users/michels/phylooptim/pkg/R/geiger")

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
it <- 50  	                         #Number of iterations
z <- length(well)                        #Number of optimizers
MIN=2                                    #min upper value
MAX=500                                  #max upper value
ub <- seq(MIN,MAX,length=it)             #Random beta start values

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

require(geiger)
require(optimx)
require(ape)

meserr=0 #Measurement error
userstart=NULL #Starting position
badLnL=1000 #Value returned when the function to optimize is undefined
numcol<-50 #number of colors to use
precision=100 #Number of points per variable to compute
plot3d=FALSE #Surface plot. Note, doesn't work if part of the surface to be plotted is undefined

td <- treedata(read.tree("BJO.Monocot.tre"),read.delim("BJO.monocot_GS")[,3],sort=T)
ntax=length(td$phy$tip.label)
chdata<- read.delim("BJO.monocot_GS")[,3] # TIP data
tree<- td$phy# Tree
n<- length(chdata)

#----- MINIMIZE NEGATIVE LOG LIKELIHOOD

start=log(c(0.1, 1.5))

k<-3

phylogMean<-function(phyvcv, data) 
{
	o<-rep(1, length(data))
	ci<-solve(phyvcv)
	
	m1<-solve(t(o) %*% ci %*% o)
	m2<-t(o) %*% ci %*% data
	
	return(m1 %*% m2)
	
}

foo<-function(x) {
	#exp(x[1])=beta, exp(x[2])=delta
	if(dplot){
		foostore$x1<<-c(foostore$x1,x[1])
		foostore$x2<<-c(foostore$x2,x[2])
	}
	t<-deltaTree(tree, delta=exp(x[2]))
	vcv<-vcv.phylo(t)
	
	
	vv<-exp(x[1])*vcv
	
	
	diag(vv)<-diag(vv)+meserr^2
	determinantVCV<-det(vv)
	if (determinantVCV==0){ #old delta had bounds, so couldn't get very low values and so didn't get singular matrices. Now that can happen, so have to guard against it
		warning("Possibly singular variance-covariance matrix, so giving this particular parameter combination a very bad likelihood score (rather than crashing)")
		return(badLnL)
	}
	mu<-phylogMean(vv, chdata)
	mu<-rep(mu, n)
	-dmvnorm(chdata, mu, vv, log=T)
}

dplot<-TRUE
foostore=list(x1=NULL,x2=NULL)
  print(start)

f <- function(x)
{ 
y <- x[,1]
cow <- as.numeric(x[,2])

obj.list <- vector("list", length(y))
names(obj.list)<-y

po <- seq(1,12,by=1)
w <- vector("list",length(y))
for (b in c(1:length(y))){w[[b]]<-rep(NA,length(well))}
for (i in c(1:length(y))){
  for (j in c(1:length(well))){
    if (y[i]==well[j]){w[[i]][j]<-po[j]}else{w[[i]][j]<-0}}}
t <- rep(NA,length(y))
for (i in c(1:length(t))){t[i]<-sum(w[[i]])}

for (i in c(1:length(y))){

bounds=list(delta=c(0.0003,cow[i])) #Optimization bounds
  
#---- DEFAULT BOUNDS
bounds.default <- matrix(c(0.00000001, 20, 0.0000001,1, 0.000001, 1, 0.00001, 2.999999, 0.0000001, 50, -3, 0, 0.0000000001, 100, -100, 100), nrow=8, ncol=2, byrow=TRUE)
rownames(bounds.default) <- c("beta", "lambda", "kappa", "delta", "alpha", "a", "nv", "mu");
colnames(bounds.default) <- c("min", "max")

#---- USER DEFINED PARAMETER BOUNDS
if (is.null(bounds)) {
	bounds <- bounds.default       # USE DEFAULTS
}else{
	if (class(bounds)!="list"){
		stop("Please specify user defined parameter bounds as a list()")
	}else{
		specified   <- !c(is.null(bounds$beta), is.null(bounds$lambda), 
				is.null(bounds$kappa), is.null(bounds$delta),  is.null(bounds$alpha), is.null(bounds$a),
				is.null(bounds$nv), is.null(bounds$mu)
		)
		bounds.user <- matrix(c(bounds$beta, bounds$lambda, bounds$kappa, bounds$delta, bounds$alpha, bounds$a, bounds$nv, bounds$mu), 
				nrow=sum(specified), ncol=2, byrow=TRUE
		)
		rownames(bounds.user) <- c("beta", "lambda", "kappa", "delta", "alpha", "a",  "nv", "mu")[specified]
		colnames(bounds.user) <- c("min", "max")
		
		#----  SET FINAL SEARCH BOUNDS
		bounds <- bounds.default
		
		bounds[specified,] <- bounds.user     # Final Bounds
		print(bounds)
	} # END if list
}  # END user bound if loop>
#--------------------------------
#---   APPEND MODEL SETTINGS  ---
#--------------------------------
bounds <- data.frame(t(bounds))

lower=log(bounds[1,c("beta","delta")])
upper=log(bounds[2,c("beta","delta")])
  
    if (t[i]==1|2|3|4|5){
                begin.time <-proc.time()
		o <- tryCatch(optimx(foo, p=start, lower=unlist(lower), upper=unlist(upper),
				     method=y[i]), error=function(x){
                                                                     x$fvalues$fvalues <- NA
                                                                     x$par$par[1] <- NA
                                                                     x$par$par[2] <- NA
                                                                     return(x)})
		lbb <- lower$beta;ubb <- upper$beta;lbd <- lower$delta;ubd <-upper$delta 
                time <- as.numeric(proc.time()[3]-begin.time[3])/(60)}else{
                begin.time <-proc.time()
		o <- tryCatch(optimx(foo, p=start,method=y[i]), 
                                     error=function(x){
                                                       x$fvalues$fvalues <- NA
                                                       x$par$par[1] <- NA
                                                       x$par$par[2] <- NA
                                                       return(x)})
		lbb <- NA;ubb <- NA;lbd <- NA;ubd <- NA
                time <- as.numeric(proc.time()[3]-begin.time[3])/(60)}
          
		results<-list(lnl=-o$fvalues$fvalues, beta=exp(o$par$par[1]), delta=exp(o$par$par[2]),beta.bnd=c(lbb,ubb),delta.bnds=c(lbd,ubd),time=time)

	results$aic<-2*k-2*results$lnl
	results$aicc<-2*k*(n-1)/(n-k-2)-2*results$lnl
	results$k<-k
	obj.list[[i]]<-results}
return(obj.list)
}

require(multicore)
begin.time <-proc.time()
jobs <- lapply(v, function(x) parallel(f(x),silent=TRUE))
results <- collect(jobs,wait=TRUE)
total.time <- as.numeric(proc.time()[3]-begin.time[3])/(60)
#what?

lt <- vector("list",length(well))

h <- 1
  for (i in c(1:N)){
    for (j in c(1:length(v[[i]][,1]))){
lt[[h]]<-c(unlist(v[[i]][j]),round(as.numeric(results[[i]][[j]]$lnl),5),round(as.numeric(results[[i]][[j]]$beta.start),5),round(as.numeric(results[[i]][[j]]$beta.bnd[1]),5),round(as.numeric(results[[i]][[j]]$beta.bnd[2]),5),round(results[[i]][[j]]$beta,5),round(as.numeric(results[[i]][[j]]$delta.start),5),round(as.numeric(results[[i]][[j]]$delta.bnd[1]),5),round(as.numeric(results[[i]][[j]]$delta.bnd[2]),5),round(results[[i]][[j]]$delta,5),round(results[[i]][[j]]$time,5));h <- h+1}}

ff <- data.frame(matrix(NA,ncol=length(lt[[1]]),nrow=length(lt)))
colnames(ff)<-c("name","L","bLB","bUB","b","dLB","dUB","d","time")
for (i in c(1:length(lt))){
  for (j in c(1:length(lt[[i]]))){
    ff[i,j] <- lt[[i]][j]}}
ff <- ff[order(ff$name,ff$dUB) , ]

optx <- vector("list",length(well))
names(optx) <- c("BFGS","bobyqa","CG","L-BFGS-B","Nelder-Mead","newuoa","nlm","nlminb","Rcgmin","Rvmmin","spg","ucminf")

aa <- seq(1,it*z,by=it)
bb <- seq(aa[2]-1,it*z,by=it)
for (i in c(1:length(aa))){optx[[i]] <- ff[aa[i]:bb[i],2:9]}

save.image("/Users/michels/phylooptim/pkg/R/geiger/monogeigererror.RData")

#po <- seq(1,12,by=1)
#w <- vector("list",N)
#for (a in c(1:N)){w[[a]] <- vector("list",length(v[[a]]))}
#for (a in c(1:N)){for (b in c(1:length(v[[a]]))){w[[a]][[b]]<-rep(NA,length(well))}}
#for (k in c(1:N)){
#  for (i in c(1:length(v[[k]]))){
#    for (j in c(1:length(well))){
#      if (v[[k]][i]==well[j]){w[[k]][[i]][j]<-po[j]}else{w[[k]][[i]][j]<-0}}}}
#t <- vector("list",N)
#for (a in c(1:N)){t[[a]] <- rep(NA,length(v[[a]]))}
#for (a in c(1:N)){for (i in c(1:length(t[[a]]))){t[[a]][i]<-sum(w[[a]][[i]])}}
