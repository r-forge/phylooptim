
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
p <- 4
#Lower lower bound
llb <- .0002
#upper lower bound
ulb <- .001
#lower upper bound
lub <- 400
#upper upper bound
uub <- 900
#lower start value
ls <- .00025
#upper start value
us <- .0011

#This section splits up the number of points evenly between each core
    U <- seq(lub,uub,length.out=p)
    L <- seq(llb,ulb,length.out=p)
    S <- seq(ls,us,length.out=p)
    MM <- (uub-lub)/N
    M <- c(lub,rep(NA,N))
    for (i in c(2:(N+1))){M[i]<-M[i-1]+MM}
    LOW<-c(lub,rep(NA,N*2-2),uub)

    j<-2
    for (m in c(2:N)){
    for (i in c(1:p)){if ((U[i]<= M[m]) && (M[m] <= U[i+1])){LOW[j]<-U[i];LOW[j+1]<-U[i+1]}}
    j<-j+2}

    buk <- vector("list", N)
    a<-1
    for(j in seq(1,length(LOW),by=2)){
    k<-.5*j+.5
    i<-a
    while(U[i] <= LOW[j+1] && U[i] >= LOW[j]){
      buk[[k]][i]<-U[i];
      i<-i+1;a<-i;if (a>p){break}}  
    }
    for (i in c(1:N)){buk[[i]]<-na.omit(buk[[i]])}
    ww <- rep(NA,N)
    for (i in c(1:N)){ww[i]<-length(buk[[i]])}

    j<-1
    k <- vector("list",N)
    for(i in c(1:N))
    {
      k[[i]]<-c(llb,ulb,LOW[j],LOW[j+1],ww[i],p,ls,us);j<-j+2
    }

library(geiger)
library(optimx)
data(geospiza)
attach(geospiza)

meserr=0 #Measurement error
userstart=NULL #Starting position
badLnL=1000 #Value returned when the function to optimize is undefined
#precision=10 #Number of points per variable to compute

td<-treedata(geospiza.tree,geospiza.data[,1],sort=T)
ntax=length(td$phy$tip.label)

#Setting up the range upon which we want to plot the function
x1range<-seq(-15,10,length.out=p) #beta
x2range<-seq(-5,8,length.out=p) #delta

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

ls <- x[7] # - lower start value
us <- x[8] # - upper start value

u <- seq(lu,uu,length.out=w)
l <- seq(ll,ul,length.out=p)
s <- seq(ls,us,length.out=p)
#t <- rep(NA,length(k)/5)
#  for (i in k)
#  {
#t[i/5]<-k[i/5]*5
#  }

o.list <- vector("list", length(s))
for (k in c(1:(length(s))))
  {
  o.list[[k]]<-vector("list", length(l))
  for (j in c(1:length(l)))
  {
    o.list[[k]][[j]]<-vector("list", length(u))
  }
}

for (k in c(1:length(s)))
  {
    for (j in c(1:length(l)))
      {
        for (i in c(1:length(u)))
          {

if (s[k]>l[j]){
beta.start<-var(chdata)/max(branching.times(tree))
#beta.start <- b[j]
start=log(c(beta.start, s[k]))

if(!is.null(userstart)) {
	start<-log(userstart)
}
        
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
names(o_all)<-c(l[j],u[i],s[k])
o.list[[k]][[j]][[i]]<-o_all
}else{d<-matrix(NA,ncol=10,nrow=6);d[,3][[1]][1]<-"A";
names(d)<-c(l[j],u[i],s[k]);o.list[[k]][[j]][[i]]<-d}
    }
  }
}
return(o.list)
}

#Think about how the start value must be BIGGER than the lower bound, so use an if then statement in the function to deal with that.

require(multicore)
jobs <- lapply(k, function(x) parallel(f(x),silent=TRUE))
results <- collect(jobs,wait=TRUE)

#for (h in c(1:N)){
#  for (v in c(1:length(S))){
#    for (i in c(1:p)){
#      for (j in c(1:ww[h])){
#if(abs(results[[h]][[v]][[i]][[j]][,1][[1]][1])<1){next}else{d<-matrix(NA,ncol=10,nrow=6);d[,1][[1]][1]<-"A";
#names(d)<-c(L[j],U[i],S[v]);results[[h]][[v]][[i]][[j]]<-d}
#}}}}

#results<-res

#Time in hours
total.time <- as.numeric(proc.time()[3]-begin.time[3])/(60*60)

bobyqac <- vector("list", length(S))
names(bobyqac) <- S
bfgsc <- vector("list", length(S))
names(bfgsc) <- S
nlminbc <- vector("list", length(S))
names(nlminbc) <- S
rcgminc <- vector("list", length(S))
names(rcgminc) <- S
rvmminc <- vector("list", length(S))
names(rvmminc) <- S
spgc <- vector("list", length(S))
names(spgc) <- S
bobyqaw <- vector("list", length(S))
names(bobyqaw) <- S
bfgsw <- vector("list", length(S))
names(bfgsw) <- S
nlminbw <- vector("list", length(S))
names(nlminbw) <- S
rcgminw <- vector("list", length(S))
names(rcgminw) <- S
rvmminw <- vector("list", length(S))
names(rvmminw) <- S
spgw <- vector("list", length(S))
names(spgw) <- S
bobyqa <- vector("list", length(S))
names(bobyqa) <- S
bfgs <- vector("list", length(S))
names(bfgs) <- S
nlminb <- vector("list", length(S))
names(nlminb) <- S
rcgmin <- vector("list", length(S))
names(rcgmin) <- S
rvmmin <- vector("list", length(S))
names(rvmmin) <- S
spg <- vector("list", length(S))
names(spg) <- S

for (v in c(1:length(S)))
{
bobyqa[[v]] <- matrix(NA,nrow=p,ncol=p)
bfgs[[v]] <- matrix(NA,p,ncol=p)
nlminb[[v]] <- matrix(NA,nrow=p,ncol=p)
rcgmin[[v]] <- matrix(NA,nrow=p,ncol=p)
rvmmin[[v]] <- matrix(NA,nrow=p,ncol=p)
spg[[v]] <- matrix(NA,nrow=p,ncol=p)
bobyqac[[v]] <- matrix(NA,nrow=p*p,ncol=2)
colnames(bobyqac[[v]]) <- c("lower","upper")
bfgsc[[v]] <- matrix(NA,nrow=p*p,ncol=2)
colnames(bfgsc[[v]]) <- c("lower","upper")
nlminbc[[v]] <- matrix(NA,nrow=p*p,ncol=2)
colnames(nlminbc[[v]]) <- c("lower","upper")
rcgminc[[v]] <- matrix(NA,nrow=p*p,ncol=2)
colnames(rcgminc[[v]]) <- c("lower","upper")
rvmminc[[v]] <- matrix(NA,nrow=p*p,ncol=2)
colnames(rvmminc[[v]]) <- c("lower","upper")
spgc[[v]] <- matrix(NA,nrow=p*p,ncol=2) 
colnames(spgc[[v]]) <- c("lower","upper")
bobyqaw[[v]] <- matrix(NA,nrow=p*p,ncol=2)
colnames(bobyqaw[[v]]) <- c("lower","upper")
bfgsw[[v]] <- matrix(NA,nrow=p*p,ncol=2)
colnames(bfgsw[[v]]) <- c("lower","upper")
nlminbw[[v]] <- matrix(NA,nrow=p*p,ncol=2)
colnames(nlminbw[[v]]) <- c("lower","upper")
rcgminw[[v]] <- matrix(NA,nrow=p*p,ncol=2)
colnames(rcgminw[[v]]) <- c("lower","upper")
rvmminw[[v]] <- matrix(NA,nrow=p*p,ncol=2)
colnames(rvmminw[[v]]) <- c("lower","upper")
spgw[[v]] <- matrix(NA,nrow=p*p,ncol=2) 
colnames(spgw[[v]]) <- c("lower","upper")
}

b<-rep(0,N+1);for (a in c(1:N)){b[a+1]<-ww[a]+b[a]};

#m<-1
#while (m <= p*p){
v<-1
while (v<=length(S)){
m<-1
  for (h in c(1:N)){
    for (i in c(1:p)){
      for (j in c(1:ww[h])){
if (results[[h]][[v]][[i]][[j]][,3][[1]][1]=="A"){
bobyqa[[v]][i,j+b[h]]<-NA
bfgs[[v]][i,j+b[h]]<-NA
nlminb[[v]][i,j+b[h]]<-NA
rcgmin[[v]][i,j+b[h]]<-NA
rvmmin[[v]][i,j+b[h]]<-NA
spg[[v]][i,j+b[h]]<-NA
m<-m+1
}else{

P<-matrix(NA,nrow=6,ncol=4)
name<-rep(NA,6)
      for (k in c(1:6))
        {
if (is.na(results[[h]][[v]][[i]][[j]][,1][[k]][2])==TRUE){P[k,]<-c(0,results[[h]][[v]][[i]][[j]][,1][[k]][2],as.numeric(names(results[[h]][[v]][[i]][[j]])[1]),as.numeric(names(results[[h]][[v]][[i]][[j]])[2]))}else{
if (results[[h]][[v]][[i]][[j]][,1][[k]][2]<3.1){P[k,]<-c(1,results[[h]][[v]][[i]][[j]][,1][[k]][2],as.numeric(names(results[[h]][[v]][[i]][[j]])[1]),as.numeric(names(results[[h]][[v]][[i]][[j]])[2]))}else{P[k,]<-c(0,results[[h]][[v]][[i]][[j]][,1][[k]][2],as.numeric(names(results[[h]][[v]][[i]][[j]])[1]),as.numeric(names(results[[h]][[v]][[i]][[j]])[2]))}}
name[k]<-results[[h]][[v]][[i]][[j]][,3][[k]]
                }

data<-cbind(P,name)
data<-data[sort.list(data[,5]), ]

#one v-value, one i, one j
if (m>p*p){next}else{

if (data[,1][1]=="1"){bobyqac[[v]][m,] <- c(as.numeric(data[,3][1]),as.numeric(data[,4][1]))}else{bobyqaw[[v]][m,] <- c(as.numeric(data[,3][1]),as.numeric(data[,4][1]))}
bobyqa[[v]][i,j+b[h]] <- as.numeric(data[,2][1])

if (data[,1][2]=="1"){bfgsc[[v]][m,] <- c(as.numeric(data[,3][2]),as.numeric(data[,4][2]))}else{bfgsw[[v]][m,] <- c(as.numeric(data[,3][2]),as.numeric(data[,4][2]))}
bfgs[[v]][i,j+b[h]] <- as.numeric(data[,2][2])

if (data[,1][3]=="1"){nlminbc[[v]][m,] <- c(as.numeric(data[,3][3]),as.numeric(data[,4][3]))}else{nlminbw[[v]][m,] <- c(as.numeric(data[,3][3]),as.numeric(data[,4][3]))}
nlminb[[v]][i,j+b[h]] <- as.numeric(data[,2][3])

if (data[,1][4]=="1"){rcgminc[[v]][m,] <- c(as.numeric(data[,3][4]),as.numeric(data[,4][4]))}else{rcgminw[[v]][m,] <- c(as.numeric(data[,3][4]),as.numeric(data[,4][4]))}
rcgmin[[v]][i,j+b[h]] <- as.numeric(data[,2][4])

if (data[,1][5]=="1"){rvmminc[[v]][m,] <- c(as.numeric(data[,3][5]),as.numeric(data[,4][5]))}else{rvmminw[[v]][m,] <- c(as.numeric(data[,3][5]),as.numeric(data[,4][5]))}
rvmmin[[v]][i,j+b[h]] <- as.numeric(data[,2][5])

if (data[,1][6]=="1"){spgc[[v]][m,] <- c(as.numeric(data[,3][6]),as.numeric(data[,4][6]))}else{spgw[[v]][m,] <- c(as.numeric(data[,3][6]),as.numeric(data[,4][6]))}
spg[[v]][i,j+b[h]] <- as.numeric(data[,2][6])

m<-1+m}

}}}};v<-v+1}

for (v in c(1:length(S))){
bobyqac[[v]] <- na.omit(bobyqac[[v]])
bfgsc[[v]] <- na.omit(bfgsc[[v]])
nlminbc[[v]] <- na.omit(nlminbc[[v]])
rcgminc[[v]] <- na.omit(rcgminc[[v]])
rvmminc[[v]] <- na.omit(rvmminc[[v]])
spgc[[v]] <- na.omit(spgc[[v]])
bobyqaw[[v]] <- na.omit(bobyqaw[[v]])
bfgsw[[v]] <- na.omit(bfgsw[[v]])
nlminbw[[v]] <- na.omit(nlminbw[[v]])
rcgminw[[v]] <- na.omit(rcgminw[[v]])
rvmminw[[v]] <- na.omit(rvmminw[[v]])
spgw[[v]] <- na.omit(spgw[[v]])
}

#Single Plots
for (v in c(1:length(S))){
png(file=paste("bobyqa start",S[v],".png"))
plot(bobyqac[[v]][,2],bobyqac[[v]][,1],type="p",xlim=c(min(U),max(U)),ylim=c(min(L),max(L)),
xlab='upper delta bound',ylab='lower delta bounds',main=paste("bobyqa:  start: ",S[v],", for changing delta bounds"),cex=.5)
points(bobyqaw[[v]][,2],bobyqaw[[v]][,1],pch=2,col="red")
dev.off()

png(file=paste("bobyqa countour start",S[v],".png"))
filled.contour(U, L, t(bobyqa[[v]]), color = terrain.colors,
    plot.title = title(main = paste("bobyqa:  start: ",S[v],", for changing delta bounds"),
    xlab = "upper delta bound", ylab = "lower delta bound"),
)
dev.off()

png(file=paste("L-BFGS-B start",S[v],".png"))
plot(bfgsc[[v]][,2],bfgsc[[v]][,1],type="p",xlim=c(min(U),max(U)),ylim=c(min(L),max(L)),
xlab='upper delta bound',ylab='lower delta bounds',main=paste("L-BFGS-B:  start: ",S[v],", for changing delta bounds"),cex=.5)
points(bfgsw[[v]][,2],bfgsw[[v]][,1],pch=2,col="red")
dev.off()

png(file=paste("L-BFGS-B countour start",S[v],".png"))
filled.contour(U, L, t(bfgs[[v]]), color = terrain.colors,
    plot.title = title(main = paste("L-BFGS-B:  start: ",S[v],", for changing delta bounds"),
    xlab = "upper delta bound", ylab = "lower delta bound"),
)
dev.off()

png(file=paste("nlminb start",S[v],".png"))
plot(nlminbc[[v]][,2],nlminbc[[v]][,1],type="p",xlim=c(min(U),max(U)),ylim=c(min(L),max(L)),
xlab='upper delta bound',ylab='lower delta bounds',main=paste("nlminb:  start: ",S[v],", for changing delta bounds"),cex=.5)
points(nlminbw[[v]][,2],nlminbw[[v]][,1],pch=2,col="red")
dev.off()

png(file=paste("nlminb countour start",S[v],".png"))
filled.contour(U, L, t(nlminb[[v]]), color = terrain.colors,
    plot.title = title(main = paste("nlminb:  start: ",S[v],", for changing delta bounds"),
    xlab = "upper delta bound", ylab = "lower delta bound"),
)
dev.off()
 
png(file=paste("rcgmin start",S[v],".png"))
plot(rcgminc[[v]][,2],rcgminc[[v]][,1],type="p",xlim=c(min(U),max(U)),ylim=c(min(L),max(L)),
xlab='upper delta bound',ylab='lower delta bounds',main=paste("rcgmin:  start: ",S[v],", for changing delta bounds"),cex=.5)
points(rcgminw[[v]][,2],rcgminw[[v]][,1],pch=2,col="red")
dev.off()

png(file=paste("rcgmin countour start",S[v],".png"))
filled.contour(U, L, t(rcgmin[[v]]), color = terrain.colors,
    plot.title = title(main = paste("rcgmin:  start: ",S[v],", for changing delta bounds"),
    xlab = "upper delta bound", ylab = "lower delta bound"),
)
dev.off()

png(file=paste("rvmmin start",S[v],".png"))
plot(rvmminc[[v]][,2],rvmminc[[v]][,1],type="p",xlim=c(min(U),max(U)),ylim=c(min(L),max(L)),
xlab='upper delta bound',ylab='lower delta bounds',main=paste("rvmmin:  start: ",S[v],", for changing delta bounds"),cex=.5)
points(rvmminw[[v]][,2],rvmminw[[v]][,1],pch=2,col="red")
dev.off()

png(file=paste("rvmmin countour start",S[v],".png"))
filled.contour(U, L, t(rvmmin[[v]]), color = terrain.colors,
    plot.title = title(main = paste("rvmmin:  start: ",S[v],", for changing delta bounds"),
    xlab = "upper delta bound", ylab = "lower delta bound"),
)
dev.off()

png(file=paste("spg start",S[v],".png"))
plot(spgc[[v]][,2],spgc[[v]][,1],type="p",xlim=c(min(U),max(U)),ylim=c(min(L),max(L)),
xlab='upper delta bound',ylab='lower delta bounds',main=paste("spg:  start: ",S[v],", for changing delta bounds"),cex=.5)
points(spgw[[v]][,2],spgw[[v]][,1],pch=2,col="red")
dev.off()

png(file=paste("spg countour start",S[v],".png"))
filled.contour(U, L, t(spg[[v]]), color = terrain.colors,
    plot.title = title(main = paste("spg:  start: ",S[v],", for changing delta bounds"),
    xlab = "upper delta bound", ylab = "lower delta bound"),
)
dev.off()
}

if (p<=4){k<-2}else{if (p<=9){k<-3}else{k<-4}}

png(file="bobyqa.png")
par(mfrow=c(k,k))
for (v in c(1:length(S)))
  {
plot(bobyqac[[v]][,2],bobyqac[[v]][,1],type="p",xlim=c(min(U),max(U)),ylim=c(min(L),max(L)),
xlab='',ylab='',main=paste("Start",round(S[v],3)),cex=.5)
points(bobyqaw[[v]][,2],bobyqaw[[v]][,1],pch=2,col="red")
}
mtext("bobyqa",at= 250, line = 22,cex=3)
dev.off()

png(file="bfgs.png")
par(mfrow=c(k,k))
for (v in c(1:length(S)))
  {
plot(bfgsc[[v]][,2],bfgsc[[v]][,1],type="p",xlim=c(min(U),max(U)),ylim=c(min(L),max(L)),
xlab='',ylab='',main=paste("Start",round(S[v],3)),cex=.5)
points(bfgsw[[v]][,2],bfgsw[[v]][,1],pch=2,col="red")
}
mtext("bfgs",at= 250, line = 22,cex=3)
dev.off()

png(file="nlminb.png")
par(mfrow=c(k,k))
for (v in c(1:length(S)))
  {
plot(nlminbc[[v]][,2],nlminbc[[v]][,1],type="p",xlim=c(min(U),max(U)),ylim=c(min(L),max(L)),
xlab='',ylab='',main=paste("Start",round(S[v],3)),cex=.5)
points(nlminbw[[v]][,2],nlminbw[[v]][,1],pch=2,col="red")
}
mtext("nlminb",at= 250, line = 22,cex=3)
dev.off()

png(file="rcgmin.png")
par(mfrow=c(k,k))
for (v in c(1:length(S)))
  {
plot(rcgminc[[v]][,2],rcgminc[[v]][,1],type="p",xlim=c(min(U),max(U)),ylim=c(min(L),max(L)),
xlab='',ylab='',main=paste("Start",round(S[v],3)),cex=.5)
points(rcgminw[[v]][,2],rcgminw[[v]][,1],pch=2,col="red")
}
mtext("rcgmin",at= 250, line = 22,cex=3)
dev.off()

png(file="rvmmin.png")
par(mfrow=c(k,k))
for (v in c(1:length(S)))
  {
plot(rvmminc[[v]][,2],rvmminc[[v]][,1],type="p",xlim=c(min(U),max(U)),ylim=c(min(L),max(L)),
xlab='',ylab='',main=paste("Start",round(S[v],3)),cex=.5)
points(rvmminw[[v]][,2],rvmminw[[v]][,1],pch=2,col="red")
}
mtext("rvmmin",at= 250, line = 22,cex=3)
dev.off()

png(file="spg.png")
par(mfrow=c(k,k))
for (v in c(1:length(S)))
  {
plot(spgc[[v]][,2],spgc[[v]][,1],type="p",xlim=c(min(U),max(U)),ylim=c(min(L),max(L)),
xlab='',ylab='',main=paste("Start",round(S[v],3)),cex=.5)
points(spgw[[v]][,2],spgw[[v]][,1],pch=2,col="red")
}
mtext("spg",at= 250, line = 22,cex=3)
dev.off()

source("http://wiki.cbr.washington.edu/qerm/sites/qerm/images/1/16/Filled.contour3.R")
source("http://wiki.cbr.washington.edu/qerm/sites/qerm/images/2/25/Filled.legend.R")

if (p<=4){
d<-vector("list", 5)
d[[1]]<-c(0.1,0.4,0.60,0.95)
d[[2]]<-c(0.5,0.8,0.60,0.95)
d[[3]]<-c(0.1,0.4,0.15,0.5)
d[[4]]<-c(0.5,0.8,0.15,0.5)
d[[5]]<-c(0.85,0.9,0.25,0.85)}else{if (p<=9){
d<-vector("list", 10)
d[[1]]<-c(0.03,0.27,0.69,0.94)
d[[2]]<-c(0.03,0.27,0.37,0.62)
d[[3]]<-c(0.03,0.27,0.06,0.31)
d[[4]]<-c(0.32,0.56,0.69,0.94)
d[[5]]<-c(0.32,0.56,0.37,0.62)
d[[6]]<-c(0.32,0.56,0.06,0.31)
d[[7]]<-c(0.60,0.84,0.69,0.94)
d[[8]]<-c(0.60,0.84,0.37,0.62)
d[[9]]<-c(0.60,0.84,0.06,0.31)
d[[10]]<-c(0.85,0.9,0.25,0.85)}else{
d<-vector("list", 17)
d[[1]]<-c(0.03,0.21,0.76,0.96)
d[[2]]<-c(0.03,0.21,0.52,0.72)
d[[3]]<-c(0.03,0.21,0.28,0.48)
d[[4]]<-c(0.03,0.21,0.04,0.24)
d[[5]]<-c(0.24,0.42,0.76,0.96)
d[[6]]<-c(0.24,0.42,0.52,0.72)
d[[7]]<-c(0.24,0.42,0.28,0.48)
d[[8]]<-c(0.24,0.42,0.04,0.24)
d[[9]]<-c(0.45,0.63,0.76,0.96)
d[[10]]<-c(0.45,0.63,0.52,0.72)
d[[11]]<-c(0.45,0.63,0.28,0.48)
d[[12]]<-c(0.45,0.63,0.04,0.24)
d[[13]]<-c(0.66,0.84,0.76,0.96)
d[[14]]<-c(0.66,0.84,0.52,0.72)
d[[15]]<-c(0.66,0.84,0.28,0.48)
d[[16]]<-c(0.66,0.84,0.04,0.24)
d[[17]]<-c(0.85,0.9,0.25,0.85)}}

library(RColorBrewer)
library(rgl)
nlevels=25
numcol<-50
levels=pretty(range(spg,na.rm=TRUE),nlevels,finite=TRUE)
k<-pretty(trunc(range(spg,na.rm=T)),nlevels)
wr.pal<-colorRampPalette(c("red","green"))
wr<-wr.pal(round(numcol*sum(k>=0)/length(k)))
bw.pal<-colorRampPalette(c("blue","light yellow"))
bw<-bw.pal(round(numcol*sum(k<=0)/length(k)))
cols<-c(bw,wr)
rcols<-cols[round((1:(length(levels)))*length(cols)/length(levels))]

png(file="bobyqa_contour.png")
plot.new()
for (v in c(1:length(S))){
par(new = "TRUE",plt = d[[v]],las = 1,cex.axis = 1)
filled.contour3(U, L, t(spg[[v]]),col=rcols,xlab = "",ylab = "",xlim = c(min(U),max(U)),ylim = c(min(L),max(L)),zlim = range(spg,na.rm=TRUE))
}
par(new = "TRUE",plt = d[[length(d)]],las = 1,cex.axis = 1)
filled.legend(U,L,spg[[1]],col = rcols,xlab = "",ylab = "",xlim = c(min(xintercepts),max(xintercepts)),ylim = c(min(slopes),max(slopes)),zlim = range(spg,na.rm=TRUE))
dev.off()

png(file="bfgs_contour.png")
plot.new()
for (v in c(1:length(S))){
par(new = "TRUE",plt = d[[v]],las = 1,cex.axis = 1)
filled.contour3(U, L, t(bfgs[[v]]),col=rcols,xlab = "",ylab = "",xlim = c(min(U),max(U)),ylim = c(min(L),max(L)),zlim = range(spg,na.rm=TRUE))
}
par(new = "TRUE",plt = d[[length(d)]],las = 1,cex.axis = 1)
filled.legend(U,L,bfgs[[1]],col = rcols,xlab = "",ylab = "",xlim = c(min(xintercepts),max(xintercepts)),ylim = c(min(slopes),max(slopes)),zlim = range(spg,na.rm=TRUE))
dev.off()

png(file="nlminb_contour.png")
plot.new()
for (v in c(1:length(S))){
par(new = "TRUE",plt = d[[v]],las = 1,cex.axis = 1)
filled.contour3(U, L, t(nlminb[[v]]),col=rcols,xlab = "",ylab = "",xlim = c(min(U),max(U)),ylim = c(min(L),max(L)),zlim = range(spg,na.rm=TRUE))
}
par(new = "TRUE",plt = d[[length(d)]],las = 1,cex.axis = 1)
filled.legend(U,L,nlminb[[1]],col = rcols,xlab = "",ylab = "",xlim = c(min(xintercepts),max(xintercepts)),ylim = c(min(slopes),max(slopes)),zlim = range(spg,na.rm=TRUE))
dev.off()

png(file="rcgmin_contour.png")
plot.new()
for (v in c(1:length(S))){
par(new = "TRUE",plt = d[[v]],las = 1,cex.axis = 1)
filled.contour3(U, L, t(rcgmin[[v]]),col=rcols,xlab = "",ylab = "",xlim = c(min(U),max(U)),ylim = c(min(L),max(L)),zlim = range(spg,na.rm=TRUE))
}
par(new = "TRUE",plt = d[[length(d)]],las = 1,cex.axis = 1)
filled.legend(U,L,rcgmin[[1]],col = rcols,xlab = "",ylab = "",xlim = c(min(xintercepts),max(xintercepts)),ylim = c(min(slopes),max(slopes)),zlim = range(spg,na.rm=TRUE))
dev.off()

png(file="rvmmin_contour.png")
plot.new()
for (v in c(1:length(S))){
par(new = "TRUE",plt = d[[v]],las = 1,cex.axis = 1)
filled.contour3(U, L, t(rvmmin[[v]]),col=rcols,xlab = "",ylab = "",xlim = c(min(U),max(U)),ylim = c(min(L),max(L)),zlim = range(spg,na.rm=TRUE))
}
par(new = "TRUE",plt = d[[length(d)]],las = 1,cex.axis = 1)
filled.legend(U,L,rvmmin[[1]],col = rcols,xlab = "",ylab = "",xlim = c(min(xintercepts),max(xintercepts)),ylim = c(min(slopes),max(slopes)),zlim = range(spg,na.rm=TRUE))
dev.off()

png(file="spg_contour.png")
plot.new()
for (v in c(1:length(S))){
par(new = "TRUE",plt = d[[v]],las = 1,cex.axis = 1)
filled.contour3(U, L, t(spg[[v]]),col=rcols,xlab = "",ylab = "",xlim = c(min(U),max(U)),ylim = c(min(L),max(L)),zlim = range(spg,na.rm=TRUE))
}
par(new = "TRUE",plt = d[[length(d)]],las = 1,cex.axis = 1)
filled.legend(U,L,spg[[1]],col = rcols,xlab = "",ylab = "",xlim = c(min(xintercepts),max(xintercepts)),ylim = c(min(slopes),max(slopes)),zlim = range(spg,na.rm=TRUE))
dev.off()
