recomputebkg=T
if(recomputebkg || is.na(res)){
  rm(list=ls())
  recomputebkg=TRUE
}

library(rgl)
library(geiger)
library(RColorBrewer)
library(optimx)
data(geospiza)
attach(geospiza)

#Bounds for upper bound?
lu <- 400            #lower bound
uu <- 900            #upper bound

#Bounds for lower bound?
ll <- 0            #lower bound
ul <- 20            #upper bound

#Start value, must be larger then the upper bound for the lower bound
s <- 1.5             #start value

meserr=0 #Measurement error
userstart=NULL #Starting position
badLnL=1000 #Value returned when the function to optimize is undefined
numcol<-50 #number of colors to use
precision=1 #Number of points per variable to compute
plot3d=FALSE #Surface plot. Note, doesn't work if part of the surface to be plotted is undefined

td<-treedata(geospiza.tree,geospiza.data[,1],sort=T)
ntax=length(td$phy$tip.label)

#Setting up the range upon which we want to plot the function
x1range<-seq(-15,10,length.out=precision) #beta
x2range<-seq(-5,8,length.out=precision) #delta

phylogMean<-function(phyvcv, data) 
{
	o<-rep(1, length(data))
	ci<-solve(phyvcv)
	
	m1<-solve(t(o) %*% ci %*% o)
	m2<-t(o) %*% ci %*% data
	
	return(m1 %*% m2)
	
}

dplot <- FALSE

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

chdata<- geospiza.data[,1]# TIP data
tree<- td$phy# Tree
n<- length(chdata)

#So this is making the matrix to input into the contour plot, exactly how is it being formed?  What is the idea behind it?

if(recomputebkg){
# The contour function requires a matrix as its input
# here I transform the z vector into the proper matrix
 res<-t(sapply(x1range,function(x,y=x2range) sapply(y,function(z) foo(c(x,z)) )))

#As the likelihood across the range changes dramatically, I log transform the likelihoods
#However, because some of the z values are negative, I have to process them
if(any(res[!is.na(res)]<0)) {
	res[!is.na(res) & res>=0 & res<1]<-0
	res[!is.na(res) & res<0]<--1*log(-1*res[!is.na(res) & res<0])
}
res[!is.na(res) & res>0]<-log(res[!is.na(res) & res>0])
}

beta.start<-var(chdata)/max(branching.times(tree))
#beta.start <- b[j]
start=log(c(beta.start, s))

if(!is.null(userstart)) {
	start<-log(userstart)
}

u <- seq(lu,uu,length.out=precision)
l <- seq(ll,ul,length.out=precision)
#t <- rep(NA,length(k)/5)
#  for (i in k)
#  {
#t[i/5]<-k[i/5]*5
#  }
o.list <- vector("list", length(l))
for (j in c(1:(length(l))))
  {
  o.list[[j]]<-vector("list", length(u))
  }

begin.time <- proc.time()
for (j in c(1:length(l)))
  {
    for (i in c(1:length(u)))
      {

bounds=list(delta=c(l[j],u[i])) #Optimization bounds

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
#----- MINIMIZE NEGATIVE LOG LIKELIHOOD

#How was this decided?

o_all<-optimx(foo, p=start, lower=unlist(lower), upper=unlist(upper), control=list(all.methods=TRUE, trace=1))
o.list[[j]][[i]]<-o_all
}
  }
total.time <- as.numeric(proc.time()[3]-begin.time[3])/(60*60)

bobyqa <- matrix(NA,nrow=length(l),ncol=length(u))
bfgs <- matrix(NA,nrow=length(l),ncol=length(u))
nlminb <- matrix(NA,nrow=length(l),ncol=length(u))
rcgmin <- matrix(NA,nrow=length(l),ncol=length(u))
rvmmin <- matrix(NA,nrow=length(l),ncol=length(u))
spg <- matrix(NA,nrow=length(l),ncol=length(u))
bobyqac <- matrix(NA,nrow=length(l)*length(u),ncol=2)
colnames(bobyqac) <- c("lower","upper")
bfgsc <- matrix(NA,nrow=length(l)*length(u),ncol=2)
colnames(bfgsc) <- c("lower","upper")
nlminbc <- matrix(NA,nrow=length(l)*length(u),ncol=2)
colnames(nlminbc) <- c("lower","upper")
rcgminc <- matrix(NA,nrow=length(l)*length(u),ncol=2)
colnames(rcgminc) <- c("lower","upper")
rvmminc <- matrix(NA,nrow=length(l)*length(u),ncol=2)
colnames(rvmminc) <- c("lower","upper")
spgc <- matrix(NA,nrow=length(l)*length(u),ncol=2) 
colnames(spgc) <- c("lower","upper")
bobyqaw <- matrix(NA,nrow=length(l)*length(u),ncol=2)
colnames(bobyqaw) <- c("lower","upper")
bfgsw <- matrix(NA,nrow=length(l)*length(u),ncol=2)
colnames(bfgsw) <- c("lower","upper")
nlminbw <- matrix(NA,nrow=length(l)*length(u),ncol=2)
colnames(nlminbw) <- c("lower","upper")
rcgminw <- matrix(NA,nrow=length(l)*length(u),ncol=2)
colnames(rcgminw) <- c("lower","upper")
rvmminw <- matrix(NA,nrow=length(l)*length(u),ncol=2)
colnames(rvmminw) <- c("lower","upper")
spgw <- matrix(NA,nrow=length(l)*length(u),ncol=2) 
colnames(spgw) <- c("lower","upper")

m<-1

for (k in c(1:length(l))){
  for (j in c(1:length(u))){
    
p<-matrix(NA,nrow=6,ncol=2)
name<-rep(NA,6)
for (i in c(1:6)){if (o.list[[k]][[j]][,1][[i]][2]<3.1){p[i,]<-c(1,o.list[[k]][[j]][,1][[i]][2])}else{p[i,]<-c(0,o.list[[k]][[j]][,1][[i]][2])}
name[i]<-o.list[[k]][[j]]$method[[i]]
                }
data<-cbind(p,name)
data<-data[sort.list(data[,3]), ]

if (data[,1][1]=="1"){bobyqac[m,] <- c(l[k],u[j])}else{bobyqaw[m,] <- c(l[k],u[j])}
bobyqa[k,j] <- data[,2][1]
if (data[,1][2]=="1"){bfgsc[m,] <- c(l[k],u[j])}else{bfgsw[m,] <- c(l[k],u[j])}
bfgs[k,j] <- data[,2][1]
if (data[,1][3]=="1"){nlminbc[m,] <- c(l[k],u[j])}else{nlminbw[m,] <- c(l[k],u[j])}
nlminb[k,j] <- data[,2][1]
if (data[,1][4]=="1"){rcgminc[m,] <- c(l[k],u[j])}else{rcgminw[m,] <- c(l[k],u[j])}
rcgmin[k,j] <- data[,2][1]
if (data[,1][5]=="1"){rvmminc[m,] <- c(l[k],u[j])}else{rvmminw[m,] <- c(l[k],u[j])}
rvmmin[k,j] <- data[,2][1]
if (data[,1][6]=="1"){spgc[m,] <- c(l[k],u[j])}else{spgw[m,] <- c(l[k],u[j])}
spg[k,j] <- data[,2][1]
m<-1+m
     }
     }

bobyqac <- na.omit(bobyqac)
bfgsc <- na.omit(bfgsc)
nlminbc <- na.omit(nlminbc)
rcgminc <- na.omit(rcgminc)
rvmminc <- na.omit(rvmminc)
spgc <- na.omit(spgc)
bobyqaw <- na.omit(bobyqaw)
bfgsw <- na.omit(bfgsw)
nlminbw <- na.omit(nlminbw)
rcgminw <- na.omit(rcgminw)
rvmminw <- na.omit(rvmminw)
spgw <- na.omit(spgw)

#png(file="bobyqa.png")
#plot(bobyqac[,2],bobyqac[,1],type="p",
#xlab='upper',ylab='lower',main="bobyqa",cex=.5)
#points(bobyqaw[,2],bobyqaw[,1],pch=2,col="red")
#dev.off()
#
#png(file="bobyqa_contour.png")
#filled.contour(u, l, bobyqa, color = terrain.colors,
#    plot.title = title(main = "bobyqa critial values w/ changing delta bounds",
#    xlab = "upper", ylab = "lower"),
#)
#dev.off()
#
#png(file="L-BFGS-B.png")
#plot(bfgsc[,2],bfgsc[,1],type="p",
#xlab='upper',ylab='lower',main="L-BFGS-B",cex=.5)
#points(bfgsw[,2],bfgsw[,1],pch=2,col="red")
#dev.off()
#
#png(file="bfgs_contour.png")
#filled.contour(u, l, bobyqa, color = terrain.colors,
#    plot.title = title(main = "bfgs critial values w/ changing delta bounds",
#    xlab = "upper", ylab = "lower"),
#)
#dev.off()
#
#png(file="nlminb.png")
#plot(nlminbc[,2],nlminbc[,1],type="p",
#xlab='upper',ylab='lower',main="nlminb",cex=.5)
#points(nlminbw[,2],nlminbw[,1],pch=2,col="red")
#dev.off()
#
#png(file="nlminb_contour.png")
#filled.contour(u, l, bobyqa, color = terrain.colors,
#    plot.title = title(main = "nlminb critial values w/ changing delta bounds",
#    xlab = "upper", ylab = "lower"),
#)
#dev.off()
#
#png(file="rcgmin.png")
#plot(rcgminc[,2],rcgminc[,1],type="p",
#xlab='upper',ylab='lower',main="rcgmin",cex=.5)
#points(rcgminw[,2],rcgminw[,1],pch=2,col="red")
#dev.off()
#
#png(file="rcgmin_contour.png")
#filled.contour(u, l, bobyqa, color = terrain.colors,
#    plot.title = title(main = "rcgmin critial values w/ changing delta bounds",
#    xlab = "upper", ylab = "lower"),
#)
#dev.off()
#
#png(file="rvmmin.png")
#plot(rvmminc[,2],rvmminc[,1],type="p",
#xlab='upper',ylab='lower',main="rvmmin",cex=.5)
#points(rvmminw[,2],rvmminw[,1],pch=2,col="red")
#dev.off()
#
#png(file="rvmmin_contour.png")
#filled.contour(u, l, bobyqa, color = terrain.colors,
#    plot.title = title(main = "rvmmin critial values w/ changing delta bounds",
#    xlab = "upper", ylab = "lower"),
#)
#dev.off()
#
#png(file="spg.png")
#plot(spgc[,2],spgc[,1],type="p",
#xlab='upper',ylab='lower',main="spg",cex=.5)
#points(spgw[,2],spgw[,1],pch=2,col="red")
#dev.off()
#
#png(file="spg_contour.png")
#filled.contour(u, l, bobyqa, color = terrain.colors,
#    plot.title = title(main = "spg critial values w/ changing delta bounds",
#    xlab = "upper", ylab = "lower"),
#)
#dev.off()
