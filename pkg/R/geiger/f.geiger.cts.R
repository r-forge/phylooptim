f <- function(x)
{ 
y <- x[,1]
cow <- as.numeric(x[,2])

obj.list <- vector("list", length(y))
names(obj.list)<-y

po <- seq(1,length(well),by=1)
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
          
		results<-list(lnl=-o$fvalues$fvalues, beta=exp(o$par$par[1]), delta=exp(o$par$par[2]),beta.bnd=c(lbb,ubb),delta.bnds=c(lbd,ubd),conv=o$conv$conv,time=time)

	results$aic<-2*k-2*results$lnl
	results$aicc<-2*k*(n-1)/(n-k-2)-2*results$lnl
	results$k<-k
	obj.list[[i]]<-results}
return(obj.list)
}
