
begin.time <- proc.time()

detectCores <- function(all.tests = FALSE)
{
    systems <-
        list(darwin  = "/usr/sbin/sysctl -n hw.ncpu 2>/dev/null",
             freebsd = "/sbin/sysctl -n hw.ncpu 2>/dev/null",
             linux   = "grep processor /proc/cpuinfo 2>/dev/null | wc -l",
             irix    = c("hinv | grep Processors | sed 's: .*::'",
                         "hinv | grep '^Processor '| wc -l"),
#             solaris = "/usr/sbin/psrinfo -v | grep 'Status of.*processor' | wc -l"
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
N <- detectCores()
#Precision?
p <- 50
#Lower lower bound
llb <- .0002
#upper lower bound
ulb <- .001
#lower upper bound
lub <- 400
#upper upper bound
uub <- 900
#start value
s <- .375825

    U <- seq(lub,uub,length.out=p)
    L <- seq(llb,ulb,length.out=p)
    MM <- (uub-lub)/N
    M <- c(lub,rep(NA,N))
    for (i in c(2:(N+1))){M[i]<-M[i-1]+MM}
    LOW<-c(lub,rep(NA,N*2-2),uub)
    j<-2
    for (m in c(2:N)){
    for (i in c(1:p)){if ((U[i]< M[m]) && (M[m] < U[i+1])){LOW[j]<-U[i];LOW[j+1]<-U[i+1]}}
    j<-j+2}

    mm <- p/N
    nn <- floor(mm)+1
    ll <- floor(mm)
    if ((mm-ll)<0){rr <- (mm-ll)*N; ww <- c(rep(hh,rr),rep(ll,N-rr))} else {ww<-rep(mm,N)}
    j<-1
    k <- vector("list",N)
    for(i in c(1:N))
    {
      k[[i]]<-c(llb,ulb,LOW[j],LOW[j+1],ww[i],p,s);j<-j+2
    }

library(geiger)
library(optimx)
data(geospiza)
attach(geospiza)

meserr=0 #Measurement error
userstart=NULL #Starting position
badLnL=1000 #Value returned when the function to optimize is undefined
precision=10 #Number of points per variable to compute

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

f <- function(x)
{
ll <- x[1] # - lower lower bound
ul <- x[2] # - upper lower bound

lu <- x[3] # - lower upper bound
uu <- x[4] # - upper upper bound

w <- x[5] # - precision for width

p <- x[6] # - precision
s <- x[7] # - start
beta.start<-var(chdata)/max(branching.times(tree))
#beta.start <- b[j]
start=log(c(beta.start, s))

if(!is.null(userstart)) {
	start<-log(userstart)
}
    
u <- seq(lu,uu,length.out=w)
l <- seq(ll,ul,length.out=p)
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

u#How was this decided?

o_all<-optimx(foo, p=start, lower=unlist(lower), upper=unlist(upper), control=list(all.methods=TRUE, trace=1))
names(o_all)<-c(l[j],u[i])
o.list[[j]][[i]]<-o_all
}
  }
return(o.list)
}

require(multicore)
jobs <- lapply(k, function(x) parallel(f(x),silent=TRUE))
results <- collect(jobs,wait=TRUE)

#Time in hours
total.time <- as.numeric(proc.time()[3]-begin.time[3])/(60*60)

bobyqa <- matrix(NA,nrow=p,ncol=p)
bfgs <- matrix(NA,p,ncol=p)
nlminb <- matrix(NA,nrow=p,ncol=p)
rcgmin <- matrix(NA,nrow=p,ncol=p)
rvmmin <- matrix(NA,nrow=p,ncol=p)
spg <- matrix(NA,nrow=p,ncol=p)
bobyqac <- matrix(NA,nrow=p*p,ncol=2)
colnames(bobyqac) <- c("lower","upper")
bfgsc <- matrix(NA,nrow=p*p,ncol=2)
colnames(bfgsc) <- c("lower","upper")
nlminbc <- matrix(NA,nrow=p*p,ncol=2)
colnames(nlminbc) <- c("lower","upper")
rcgminc <- matrix(NA,nrow=p*p,ncol=2)
colnames(rcgminc) <- c("lower","upper")
rvmminc <- matrix(NA,nrow=p*p,ncol=2)
colnames(rvmminc) <- c("lower","upper")
spgc <- matrix(NA,nrow=p*p,ncol=2) 
colnames(spgc) <- c("lower","upper")
bobyqaw <- matrix(NA,nrow=p*p,ncol=2)
colnames(bobyqaw) <- c("lower","upper")
bfgsw <- matrix(NA,nrow=p*p,ncol=2)
colnames(bfgsw) <- c("lower","upper")
nlminbw <- matrix(NA,nrow=p*p,ncol=2)
colnames(nlminbw) <- c("lower","upper")
rcgminw <- matrix(NA,nrow=p*p,ncol=2)
colnames(rcgminw) <- c("lower","upper")
rvmminw <- matrix(NA,nrow=p*p,ncol=2)
colnames(rvmminw) <- c("lower","upper")
spgw <- matrix(NA,nrow=p*p,ncol=2) 
colnames(spgw) <- c("lower","upper")

m<-1

for (h in c(1:N)){
  for (i in c(1:p)){
    for (j in c(1:ww[h])){
    
P<-matrix(NA,nrow=6,ncol=2)
name<-rep(NA,6)
      for (k in c(1:6))
        {
if (results[[h]][[i]][[j]][,1][[k]][2]<3.1){P[k,]<-c(1,results[[h]][[i]][[j]][,1][[k]][2])}else{P[k,]<-c(0,results[[h]][[i]][[j]][,1][[k]][2])}
name[k]<-results[[h]][[i]][[j]][,3][[k]]
                }
data<-cbind(P,name)
data<-data[sort.list(data[,3]), ]

if (h==1){
if (data[,1][1]=="1"){bobyqac[m,] <- c(as.numeric(names(results[[h]][[i]][[j]])[1]),as.numeric(names(results[[h]][[i]][[j]])[2]))}else{bobyqaw[m,] <- c(as.numeric(names(results[[h]][[i]][[j]])[1]),as.numeric(names(results[[h]][[i]][[j]])[2]))}
bobyqa[i,j] <- data[,2][1]
if (data[,1][2]=="1"){bfgsc[m,] <- c(as.numeric(names(results[[h]][[i]][[j]])[1]),as.numeric(names(results[[h]][[i]][[j]])[2]))}else{bfgsw[m,] <- c(as.numeric(names(results[[h]][[i]][[j]])[1]),as.numeric(names(results[[h]][[i]][[j]])[2]))}
bfgs[i,j] <- data[,2][2]
if (data[,1][3]=="1"){nlminbc[m,] <- c(as.numeric(names(results[[h]][[i]][[j]])[1]),as.numeric(names(results[[h]][[i]][[j]])[2]))}else{nlminbw[m,] <- c(as.numeric(names(results[[h]][[i]][[j]])[1]),as.numeric(names(results[[h]][[i]][[j]])[2]))}
nlminb[i,j] <- data[,2][3]
if (data[,1][4]=="1"){rcgminc[m,] <- c(as.numeric(names(results[[h]][[i]][[j]])[1]),as.numeric(names(results[[h]][[i]][[j]])[2]))}else{rcgminw[m,] <- c(as.numeric(names(results[[h]][[i]][[j]])[1]),as.numeric(names(results[[h]][[i]][[j]])[2]))}
rcgmin[i,j] <- data[,2][4]
if (data[,1][5]=="1"){rvmminc[m,] <- c(as.numeric(names(results[[h]][[i]][[j]])[1]),as.numeric(names(results[[h]][[i]][[j]])[2]))}else{rvmminw[m,] <- c(as.numeric(names(results[[h]][[i]][[j]])[1]),as.numeric(names(results[[h]][[i]][[j]])[2]))}
rvmmin[i,j] <- data[,2][5]
if (data[,1][6]=="1"){spgc[m,] <- c(as.numeric(names(results[[h]][[i]][[j]])[1]),as.numeric(names(results[[h]][[i]][[j]])[2]))}else{spgw[m,] <- c(as.numeric(names(results[[h]][[i]][[j]])[1]),as.numeric(names(results[[h]][[i]][[j]])[2]))}
spg[i,j] <- data[,2][6]}else{b<-rep(0,N+1);for (a in c(1:2)){b[a+1]<-ww[a]+b[a]};
if (data[,1][1]=="1"){bobyqac[m,] <- c(as.numeric(names(results[[h]][[i]][[j]])[1]),as.numeric(names(results[[h]][[i]][[j]])[2]))}else{bobyqaw[m,] <- c(as.numeric(names(results[[h]][[i]][[j]])[1]),as.numeric(names(results[[h]][[i]][[j]])[2]))}
bobyqa[i,j+b[h]] <- data[,2][1]
if (data[,1][2]=="1"){bfgsc[m,] <- c(as.numeric(names(results[[h]][[i]][[j]])[1]),as.numeric(names(results[[h]][[i]][[j]])[2]))}else{bfgsw[m,] <- c(as.numeric(names(results[[h]][[i]][[j]])[1]),as.numeric(names(results[[h]][[i]][[j]])[2]))}
bfgs[i,j+b[h]] <- data[,2][2]
if (data[,1][3]=="1"){nlminbc[m,] <- c(as.numeric(names(results[[h]][[i]][[j]])[1]),as.numeric(names(results[[h]][[i]][[j]])[2]))}else{nlminbw[m,] <- c(as.numeric(names(results[[h]][[i]][[j]])[1]),as.numeric(names(results[[h]][[i]][[j]])[2]))}
nlminb[i,j+b[h]] <- data[,2][3]
if (data[,1][4]=="1"){rcgminc[m,] <- c(as.numeric(names(results[[h]][[i]][[j]])[1]),as.numeric(names(results[[h]][[i]][[j]])[2]))}else{rcgminw[m,] <- c(as.numeric(names(results[[h]][[i]][[j]])[1]),as.numeric(names(results[[h]][[i]][[j]])[2]))}
rcgmin[i,j+b[h]] <- data[,2][4]
if (data[,1][5]=="1"){rvmminc[m,] <- c(as.numeric(names(results[[h]][[i]][[j]])[1]),as.numeric(names(results[[h]][[i]][[j]])[2]))}else{rvmminw[m,] <- c(as.numeric(names(results[[h]][[i]][[j]])[1]),as.numeric(names(results[[h]][[i]][[j]])[2]))}
rvmmin[i,j+b[h]] <- data[,2][5]
if (data[,1][6]=="1"){spgc[m,] <- c(as.numeric(names(results[[h]][[i]][[j]])[1]),as.numeric(names(results[[h]][[i]][[j]])[2]))}else{spgw[m,] <- c(as.numeric(names(results[[h]][[i]][[j]])[1]),as.numeric(names(results[[h]][[i]][[j]])[2]))}
spg[i,j+b[h]] <- data[,2][6]}
  
m<-1+m
     }
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

png(file="bobyqa.png")
plot(bobyqac[,2],bobyqac[,1],type="p",xlim=c(min(U),max(U)),ylim=c(min(L),max(L)),
xlab='upper',ylab='lower',main="bobyqa",cex=.5)
points(bobyqaw[,2],bobyqaw[,1],pch=2,col="red")
dev.off()

png(file="bobyqa_contour.png")
filled.contour(U, L, t(bobyqa), color = terrain.colors,
    plot.title = title(main = "bobyqa critial values w/ changing delta bounds",
    xlab = "upper", ylab = "lower"),
)
dev.off()

png(file="L-BFGS-B.png")
plot(bfgsc[,2],bfgsc[,1],type="p",xlim=c(min(U),max(U)),ylim=c(min(L),max(L)),
xlab='upper',ylab='lower',main="L-BFGS-B",cex=.5)
points(bfgsw[,2],bfgsw[,1],pch=2,col="red")
dev.off()

png(file="bfgs_contour.png")
filled.contour(U, L, t(bfgs), color = terrain.colors,
    plot.title = title(main = "bfgs critial values w/ changing delta bounds",
    xlab = "upper", ylab = "lower"),
)
dev.off()

png(file="nlminb.png")
plot(nlminbc[,2],nlminbc[,1],type="p",xlim=c(min(U),max(U)),ylim=c(min(L),max(L)),
xlab='upper',ylab='lower',main="nlminb",cex=.5)
points(nlminbw[,2],nlminbw[,1],pch=2,col="red")
dev.off()

png(file="nlminb_contour.png")
filled.contour(U, L, t(nlminb), color = terrain.colors,
    plot.title = title(main = "nlminb critial values w/ changing delta bounds",
    xlab = "upper", ylab = "lower"),
)
dev.off()

png(file="rcgmin.png")
plot(rcgminc[,2],rcgminc[,1],type="p",xlim=c(min(U),max(U)),ylim=c(min(L),max(L)),
xlab='upper',ylab='lower',main="rcgmin",cex=.5)
points(rcgminw[,2],rcgminw[,1],pch=2,col="red")
dev.off()

png(file="rcgmin_contour.png")
filled.contour(U, L, t(rcgmin), color = terrain.colors,
    plot.title = title(main = "rcgmin critial values w/ changing delta bounds",
    xlab = "upper", ylab = "lower"),
)
dev.off()

png(file="rvmmin.png")
plot(rvmminc[,2],rvmminc[,1],type="p",xlim=c(min(U),max(U)),ylim=c(min(L),max(L)),
xlab='upper',ylab='lower',main="rvmmin",cex=.5)
points(rvmminw[,2],rvmminw[,1],pch=2,col="red")
dev.off()

png(file="rvmmin_contour.png")
filled.contour(U, L, t(rvmmin), color = terrain.colors,
    plot.title = title(main = "rvmmin critial values w/ changing delta bounds",
    xlab = "upper", ylab = "lower"),
)
dev.off()

png(file="spg.png")
plot(spgc[,2],spgc[,1],type="p",xlim=c(min(U),max(U)),ylim=c(min(L),max(L)),
xlab='upper',ylab='lower',main="spg",cex=.5)
points(spgw[,2],spgw[,1],pch=2,col="red")
dev.off()

png(file="spg_contour.png")
filled.contour(U, L, t(spg), color = terrain.colors,
    plot.title = title(main = "spg critial values w/ changing delta bounds",
    xlab = "upper", ylab = "lower"),
)
dev.off()
