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

bounds=list(delta=c(exp(-18.42068),exp(6.304337))) #Optimization bounds
meserr=0 #Measurement error
userstart=NULL #Starting position
badLnL=1000 #Value returned when the function to optimize is undefined
numcol<-50 #number of colors to use
precision=100 #Number of points per variable to compute
plot3d=FALSE #Surface plot. Note, doesn't work if part of the surface to be plotted is undefined

td<-treedata(geospiza.tree,geospiza.data[,1],sort=T)
ntax=length(td$phy$tip.label)

#Setting up the range upon which we want to plot the function
x1range<-seq(-5,10,length.out=precision) #beta
x2range<-seq(-2,8,length.out=precision) #delta

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

#--------------------------------
#---        FIT MODEL         ---
#--------------------------------

chdata<- geospiza.data[,1]# TIP data
tree<- td$phy# Tree
n<- length(chdata)

#----- MINIMIZE NEGATIVE LOG LIKELIHOOD

#How was this decided?
beta.start<-var(chdata)/max(branching.times(tree))

start=log(c(beta.start, 0.1))

if(!is.null(userstart)) {
	start<-log(userstart)
}

lower=log(bounds[1,c("beta","delta")])
upper=log(bounds[2,c("beta","delta")])

#--- WHAT do these functions DO?

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
        if (is.nan(determinantVCV)){determinantVCV=0}
	if (determinantVCV==0){ #old delta had bounds, so couldn't get very low values and so didn't get singular matrices. Now that can happen, so have to guard against it
		warning("Possibly singular variance-covariance matrix, so giving this particular parameter combination a very bad likelihood score (rather than crashing)")
		return(badLnL)
	}
	mu<-phylogMean(vv, chdata)
	mu<-rep(mu, n)
	-dmvnorm(chdata, mu, vv, log=T)
}

dplot<-FALSE

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

#Again help?

#Turns on tracking foo into foostore
dplot<-TRUE
foostore=list(x1=NULL,x2=NULL)
  print(start)

wellt <- c("spg", "Rcgmin", "Rvmmin", "bobyqa","L-BFGS-B",1,1,1,1,1,1,1)
well <- c("spg", "Rcgmin", "Rvmmin", "bobyqa","L-BFGS-B","nlminb","ucminf","Nelder-Mead","nlm","CG","BFGS","newuoa")
obj.list <- vector("list", length(well))
names(obj.list)<-well
for (i in c(1:length(obj.list))){
    if (well[i]==wellt[i]){
		obj.list[[i]] <- optimx(foo, p=start, lower=unlist(lower), upper=unlist(upper), 			            method=well[i])}else{
		obj.list[[i]] <- optimx(foo, p=start,method=well[i])}}

x<-foostore$x1
y<-foostore$x2
s <- seq(length(x)-1)# one shorter than data
#Finally I plot the figure
#image(x1range,x2range,res,col=topo.colors(250),useRaster=T)
#contour(x1range,x2range,res,add=TRUE,nlevels=25)
# arrows(x[s], y[s], x[s+1], y[s+1],length=.1,col=gray(seq(1,0,len=length(x))))
# points(o$par[1],o$par[2],col=2,pch=3)

#Contour Plot for 'L-BFGS-B'

nlevels=25
levels=pretty(range(res,na.rm=TRUE),nlevels,finite=TRUE)
k<-pretty(trunc(range(res,na.rm=T)),nlevels)
wr.pal<-colorRampPalette(c("light yellow","dark green"))
wr<-wr.pal(round(numcol*sum(k>=0)/length(k)))
bw.pal<-colorRampPalette(c("blue","light yellow"))
bw<-bw.pal(round(numcol*sum(k<=0)/length(k)))
cols<-c(bw,wr)
rcols<-cols[round((1:(length(levels)))*length(cols)/length(levels))]

png(file="contour.png")
filled.contour(x1range,x2range,res,levels=levels,col=rcols,  
		plot.axes={
			axis(1)
			axis(2)
      plot.title=title(main=paste("L-BFGS-B Delta bounds: ",round(bounds$delta[1],3),",",round(bounds$delta[2],3)),xlab="log beta",
      ylab="log delta")
			contour(x1range,x2range,res,add=TRUE,nlevels=length(levels))
			arrows(x[s], y[s], x[s+1], y[s+1],length=.1,col=gray(seq(0.5,0,len=length(x))),lwd=2)
			points(obj.list$'L-BFGS-B'$par[[1]][1],obj.list$'L-BFGS-B'$par[[1]][2],col=2,pch=3)
      points(x[1],y[1],col=2)
		}
)

f <- function(x){
xt <- table(x)
t <- as.numeric(names(xt[xt == max(xt)]))  
return(t)}

t1 <- rep(NA,length(obj.list)); t2 <- rep(NA,length(obj.list)); 
for (i in c(1:length(obj.list))){t1[i] <- round(obj.list[[i]]$par[[1]][1],5);t2[i] <- round(obj.list[[i]]$par[[1]][2],5)}
x1 <- f(t1)
y1 <- f(t2)

#Figure out zero point first, then it will work
points(x1-.75,y1, pch = 3, cex = 4, col = "yellow")
dev.off()
