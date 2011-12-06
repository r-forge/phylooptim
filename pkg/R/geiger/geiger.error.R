require(geiger)

`fitContinuous` <-
function(phy, data, data.names=NULL, model=c("BM", "OU", "lambda", "kappa", "delta", "EB", "white", "trend"), bounds=NULL,  meserr=NULL)
{
	
	# sort is T because sub-functions assume data are in
	# this particular order
	
	model<-match.arg(model)
	
	td<-treedata(phy, data, data.names, sort=T)

	ntax=length(td$phy$tip.label)

	if(is.null(meserr)) {
		me=td$data
		me[]=0
		meserr=me	
	} else if(length(meserr)==1) {
		me=td$data
		me[]=meserr
		meserr=me
	} else if(is.vector(meserr)) {
		if(!is.null(names(meserr))) {
			o<-match(rownames(td$data), names(meserr))
			if(length(o)!=ntax) stop("meserr is missing some taxa from the tree")
			meserr<-as.matrix(meserr[o,])
		} else {
			if(length(meserr)!=ntax) stop("No taxon names in meserr, and the number of taxa does not match the tree")
			me<-td$data
			me[]=meserr
			meserr=me
		}
	} else {
		if(!is.null(rownames(meserr))) {
			o<-match(rownames(td$data), rownames(meserr))
			meserr=meserr[o,]
		} else {
			if(sum(dim(meserr)!=dim(td$data))!=0)
				stop("No taxon names in meserr, and the number of taxa does not match the tree")
			print("No names in meserr; assuming that taxa are in the same order as tree")	
		}
	}

	#--------------------------------
    #---    PREPARE DATA LIST     ---
    #--------------------------------
	ds			<- list()
   		ds$tree 		<- td$phy          # TIP data 
    #--------------------------------
    #--- SET MODEL SPECIFICATIONS ---
    #--------------------------------
    cat("Fitting ", model, "model:\n")
    #-----------------------------
    #---  SET PARAMETER BOUNDS ---
    #-----------------------------
    #---- DEFAULT BOUNDS
    bounds.default			 <- matrix(c(0.00000001, 20, 0.0000001,1, 0.000001, 1, 0.00001, 2.999999, 0.0000001, 50, -3, 0, 0.0000000001, 100, -100, 100), nrow=8, ncol=2, byrow=TRUE)
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
   		} # END if list
   	}  # END user bound if loop
   	#--------------------------------
    #---   APPEND MODEL SETTINGS  ---
    #--------------------------------
  	ds$bounds <- data.frame(t(bounds))

  	ds$model  <- model
  	#--------------------------------
    #---        FIT MODEL         ---
    #--------------------------------
    result<-list()
    for(i in 1:ncol(td$data)) {
    	ds$data=td$data[,i]
    	ds$meserr=meserr[,i]
  		result[[i]]<-fitContinuousModel(ds, print=print)
  		if(!is.null(colnames(td$data))) names(result)[i]<-colnames(td$data)[i] else names(result)[i]<-paste("Trait", i, sep="")

  	}
  	result
}


`fitContinuousModel` <-
function(ds, print=TRUE)
{

wellt <- c("spg", "Rcgmin", "Rvmmin", "bobyqa","L-BFGS-B",1,1,1,1,1,1,1)
well <- c("spg", "Rcgmin", "Rvmmin", "bobyqa","L-BFGS-B","nlminb","ucminf","Nelder-Mead","nlm","CG","BFGS","newuoa")
obj.list <- vector("list", length(well))
names(obj.list)<-well

  for (i in c(1:length(obj.list))){

	bounds 	<- ds$bounds
	model 	<- ds$model
	n 		<- length(ds$data)

	#----- MINIMIZE NEGATIVE LOG LIKELIHOOD
	
	beta.start<-var(ds$data)/max(branching.times(ds$tree))


	out         <- NULL
	
	y			<- ds$data				# TIP data
	tree		<- ds$tree			# Tree
	meserr		<- ds$meserr
	n			<- length(y)
	

	#----------------------------------
	#-----       DEFAULT FIT      -----
	#----------------------------------
	if (model=="BM") {
		k<-2
	
		vcv<-vcv.phylo(tree)

		start=log(beta.start, 0.1)
                lower=log(bounds[1,c("beta","delta")])
                upper=log(bounds[2,c("beta","delta")])
		
		foo<-function(x) {
			vv<-exp(x)*vcv
			diag(vv)<-diag(vv)+meserr^2
			mu<-phylogMean(vv, y)
			mu<-rep(mu, n)
			-dmvnorm(y, mu, vv, log=T)
		}
		
#		o <- optimx(foo, p=start, lower=unlist(lower), upper=unlist(upper), 			            control=list(all.methods=TRUE),trace=0)

#        	o<-o[sort.list(unlist(o[,3])), ]

#		results<-list(lnl=-o[[2]][[i]][1], beta=exp(o[[1]][[i]][1]))
		o<-optim(foo, p=start, lower=lower, upper=upper, method="L")
		results<-list(lnl=-o$value, beta= exp(o$par))

	#----------------------------------
	#-----       LAMBDA ONLY      -----
	#----------------------------------
	} else if (model=="lambda"){
		k<-3
		
		start=log(c(beta.start, 0.5))
		lower=log(bounds[1,c("beta","lambda")])
		upper=log(bounds[2,c("beta","lambda")])
		
		
		foo<-function(x) {


			vcv<-vcv.phylo(tree)

			index			<-	matrix(TRUE, n,n)
			diag(index)		<- FALSE
			vcv[index] 	<- vcv[index]*exp(x[2])
			
			vv<-exp(x[1])*vcv

			
			diag(vv)<-diag(vv)+meserr^2
			mu<-phylogMean(vv, y)
			mu<-rep(mu, n)
			
			-dmvnorm(y, mu, vv, log=T)
		}

		o <- optimx(foo, p=start, lower=unlist(lower), upper=unlist(upper), 			            control=list(all.methods=TRUE),trace=0)

        	o<-o[sort.list(unlist(o[,3])), ]

		results<-list(lnl=-o[[2]][[i]], beta=exp(o[[1]][[i]][1]), lambda=exp(o[[1]][[i]][2]))

		#o<-optim(foo, p=start, lower=lower, upper=upper, method="L")	
		#results<-list(lnl=-o$value, beta= exp(o$par[1]), lambda=exp(o$par[2]))

	#----------------------------------
	#-----        KAPPA ONLY      -----
	#----------------------------------
	} else if (model=="kappa"){
		k<-3
		start=log(c(beta.start, 0.5))
		lower=log(bounds[1,c("beta","kappa")])
		upper=log(bounds[2,c("beta","kappa")])
				
		foo<-function(x) {
                  
		t<-kappaTree(tree, kappa=exp(x[2]))
		vcv<-vcv.phylo(t)

			
		vv<-exp(x[1])*vcv

			
		diag(vv)<-diag(vv)+meserr^2
		mu<-phylogMean(vv, y)
		mu<-rep(mu, n)
			
		-dmvnorm(y, mu, vv, log=T)
		}
                
		o <- optimx(foo, p=start, lower=unlist(lower), upper=unlist(upper), 			            control=list(all.methods=TRUE),trace=0)

        	o<-o[sort.list(unlist(o[,3])), ]

		results<-list(lnl=-o[[2]][[i]], beta=exp(o[[1]][[i]][1]), lambda=exp(o[[1]][[i]][2]))

		#o<-optim(foo, p=start, lower=lower, upper=upper, method="L")
		
		#results<-list(lnl=-o$value, beta= exp(o$par[1]), lambda=exp(o$par[2]))


	#----------------------------------
	#-----        DELTA ONLY      -----
	#----------------------------------	
	} else if (model=="delta"){
		

          k<-3
		start=log(c(beta.start, 0.5))
		lower=log(bounds[1,c("beta","delta")])
		upper=log(bounds[2,c("beta","delta")])

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
if (well[i]==wellt[i]){
		o <- optimx(foo, p=start, lower=unlist(lower), upper=unlist(upper), 			            method=well[i])}else{

		o <- optimx(foo, p=start,method=well[i])}
          
		results<-list(lnl=-o$fvalues$fvalues, beta=exp(o$par$par[1]), delta=exp(o$par$par[2]),beta.bnd=c(lower$beta,upper$beta),delta.bnds=c(lower$delta,upper$delta))

		#o<-optim(foo, p=start, lower=lower, upper=upper, method="L")
		
		#results<-list(lnl=-o$value, beta= exp(o$par[1]), delta=exp(o$par[2]))	
	#----------------------------------
	#-----        WHITE NOISE     -----
	#----------------------------------	
	} else if (model=="white"){
		
		k<-2
		start=c(mean(y), log(var(y)))
		lower=c(-Inf, log(bounds[1,"nv"]))
		upper=c(Inf, log(bounds[2, "nv"]))
		
		lnl.noise<- function (p, x, se)
		# p is the vector of parameters, tree is not needed
		# x and se are trait means and std errors
		{
  			## prep parameters
  			root<- p[1]	# trait value of root of tree (also optimum)
  			vs<- exp(p[2])		# white noise variance 
  			n<- length(x)
  			VV<- diag(vs, nrow=n)	
  			diag(VV)<- diag(VV) + se^2	# add sampling variance to diagonal
    
  			## logl
  			M<- rep(root,n)
  			-dmvnorm(x, M, VV, log=TRUE)
		}
		
		o <- optimx(lnl.noise, p=start, lower=unlist(lower), upper=unlist(upper), 			            control=list(all.methods=TRUE),trace=0)

        	o<-o[sort.list(unlist(o[,3])), ]

		results<-list(lnl=-o[[2]][[i]], mean=exp(o[[1]][[i]][1]), nv=exp(o[[1]][[i]][2]))

		#o<- optim(start, fn=lnl.noise, x=y, se=meserr, lower=lower, upper=upper, method="L")
		
		#results<-list(lnl=-o$value, mean= o$par[1], nv=exp(o$par[2]))	
	#----------------------------------
	#-----        TREND           -----
	#----------------------------------	
	} else if (model=="trend"){
		
		k<-3
		vcv<- vcv.phylo(tree)
  		ww<- lm(y ~ diag(vcv))
  		p0<- c(phylogMean(vcv, y), var(y)/max(branching.times(tree)), coef(ww)[2])
		if(is.na(p0[3])) {
			p0[3]<-0
			if(is.ultrametric(tree))
				cat("WARNING: Cannot estimate a trend with an ultrametric tree; lnl will be the same as the BM model")
		}
		lower=c(-Inf, log(bounds[1,"beta"]), bounds[1,"mu"])
		upper=c(Inf, log(bounds[2,"beta"]), bounds[2,"mu"])
		

		lnl.BMtrend<- function(p, vcv, x, se)
		# p is vector of parameters, tr is tree
		# x and se are vectors of trait means and standard errors
		{
  			## prep parameters
  			root<- p[1]	# trait value of root of tree
  			vs<- exp(p[2])		# BM variance 
  			ms<- p[3]		# BM trend
  			VV<- vs*vcv	
  			diag(VV)<- diag(VV) + se^2	# add sampling variance to diagonal
  			n<- length(x)
  
 			## logl
  			M<- root+ ms*diag(vcv)
  			- dmvnorm(x, M, VV, log=TRUE)
  		}

		o <- optimx(lnl.BMtrend, p=p0, lower=unlist(lower), upper=unlist(upper), control=list(all.methods=TRUE),trace=0)

        	o<-o[sort.list(unlist(o[,3])), ]

                names(o[[1]][[i]][3])<-NULL

		results<-list(lnl=-o[[2]][[i]], mean=exp(o[[1]][[i]][1]), beta=exp(o[[1]][[i]][2]), mu=o[[1]][[i]][3])

		#o<- optim(p0, fn=lnl.BMtrend, vcv=vcv, x=y, se=meserr, lower=lower, upper=upper, method="L")
		#names(o$par)<-NULL
		#results<-list(lnl=-o$value, mean= o$par[1], beta=exp(o$par[2]), mu=o$par[3])		
	#----------------------------------
	#-----        ALPHA ONLY      -----
	#----------------------------------			
	} else if (model=="OU"){
	## modified 12 dec 07 to call ouMatrix(x) instead of vcv.phylo(ouTree(x))

		k<-3
		
		start=log(c(beta.start, 0.5))
		lower=log(bounds[1,c("beta","alpha")])
		upper=log(bounds[2,c("beta","alpha")])
	
		vcvOrig<-vcv.phylo(tree)
		foo<-function(x) {
			vcv <- ouMatrix(vcvOrig, exp(x[2]))
			
			## t<-ouTree(tree, exp(x[2]))
			##vcv<-vcv.phylo(t)
			
			vv<-exp(x[1])*vcv
			diag(vv)<-diag(vv)+meserr^2
			
			mu<-phylogMean(vv, y)
			mu<-rep(mu, n)
			
			-dmvnorm(y, mu, vv, log=T)
		}
		
		outTries<-list()
		
		# First one: try near BM solution
		start=c(log(beta.start), -50)

		o <- optimx(foo, p=start, lower=unlist(lower), upper=unlist(upper), 			            control=list(all.methods=TRUE),trace=0)

        	o <- o[sort.list(unlist(o[,3])), ]

                outTries[[1]]<-o[i,]

		#outTries[[1]]<-optim(foo, p=start, lower=lower, upper=upper, method="L")
		
		# Second one: try one with very strong constraints
		tv<-var(y)
		start=log(c(tv*2000, 1000))

		o <- optimx(foo, p=start, lower=unlist(lower), upper=unlist(upper), 			            control=list(all.methods=TRUE),trace=0)

        	o<-o[sort.list(unlist(o[,3])), ]

		outTries[[2]]<-o[i,]

		#outTries[[2]]<-optim(foo, p=start, lower=lower, upper=upper, method="L")
	

	
		# Try ten random ones
		for(j in 1:10){
			while(1) {

				lower=c(runif(2, min=-20, max=-1))
				upper=lower+runif(2, min=0, max=10)
				start=c(runif(1, min=lower[1], max=upper[1]), runif(1, min=lower[2], max=upper[2]))

		o <- optimx(foo, p=start, lower=unlist(lower), upper=unlist(upper), 			            control=list(all.methods=TRUE),trace=0)

        	o<-o[sort.list(unlist(o[,3])), ]

				te<-try(outTries[[j+2]]<-o[i,])
				if(class(te)!="try-error") break
				}
				
		}
		
		# Try range of alphas
		atry<- -5:4
		stry<- log(tv*2*exp(atry))
		for(j in 1:10){
			while(1) {

				lower=c(-20, -20)
				upper=c(10, 10)
				start=c(stry[j], atry[j])

		o <- optimx(foo, p=start, lower=unlist(lower), upper=unlist(upper), 			            control=list(all.methods=TRUE),trace=0)

        	o<-o[sort.list(unlist(o[,3])), ]

				te<-try(outTries[[j+12]]<-o[i,])
				if(class(te)!="try-error") break
				}
				
		}
		
		
		
		ntries<-22
		ltry<-numeric(ntries)
		lsol<-matrix(nrow= ntries, ncol=2)
		for(j in 1:ntries) {
				ltry[j]<-outTries[[j]]$value
				lsol[j,]<-exp(outTries[[j]]$par)
			}

		ltd<-ltry-min(ltry)
		b<-min(which(ltry==min(ltry)))

		gc<-which(ltd<0.01)
		us<-lsol[gc,1]
		usc<-sum((us-min(us))>0.01)			
		out<-outTries[[b[1]]]	
		if(usc>1) {out$message="Warning: likelihood surface is flat."}
			
		if(out$convergence!=0) {out$message="Warning: may not have converged to a proper solution."}

		results<-list(lnl=-out$value, beta= exp(out$par[1]), alpha=exp(out$par[2]), convergence=out$convergence, message=out$message, k=k)


	#----------------------------------
	#-----        EB ONLY      -----
	#----------------------------------	
	} else if(model=="EB"){

		k<-3
		start=c(log(beta.start), 0.01)
		lower=c(log(bounds[1,"beta"]),bounds[1,"a"])
		upper=c(log(bounds[2,"beta"]),bounds[2,"a"])
		
		foo<-function(x) {
			t<-exponentialchangeTree(tree, a=x[2])

			vcv<-vcv.phylo(t)
			
			vv<-exp(x[1])*vcv
			diag(vv)<-diag(vv)+meserr^2
			
			mu<-phylogMean(vv, y)
			mu<-rep(mu, n)
			
			-dmvnorm(y, mu, vv, log=T)
		}

		o <- optimx(foo, p=start, lower=unlist(lower), upper=unlist(upper), 			            control=list(all.methods=TRUE),trace=0)

        	o<-o[sort.list(unlist(o[,3])), ]

		results<-list(lnl=-o[[2]][[i]], beta=exp(o[[1]][[i]][1]), a=o[[1]][[i]][2])

		#o<-optim(foo, p=start, lower=lower, upper=upper, method="L")
			
		#results<-list(lnl=-o$value, beta= exp(o$par[1]), a=o$par[2])	
	}
	
	results$aic<-2*k-2*results$lnl
	results$aicc<-2*k*(n-1)/(n-k-2)-2*results$lnl
	results$k<-k
	obj.list[[i]]<-results
}
return(obj.list)
}



phylogMean<-function(phyvcv, data) 
{
	o<-rep(1, length(data))
	ci<-solve(phyvcv)
	
	m1<-solve(t(o) %*% ci %*% o)
	m2<-t(o) %*% ci %*% data
	
	return(m1 %*% m2)
	
	}
	
ouMatrix <- function(vcvMatrix, alpha) 
{
## follows Hansen 1997; does not assume ultrametricity (AH 12 dec 07)
## vectorized by LJH
  vcvDiag<-diag(vcvMatrix)
  diagi<-matrix(vcvDiag, nrow=length(vcvDiag), ncol=length(vcvDiag))
  diagj<-matrix(vcvDiag, nrow=length(vcvDiag), ncol=length(vcvDiag), byrow=T)

  Tij = diagi + diagj - (2 * vcvMatrix)
    
  vcvRescaled = (1 / (2 * alpha)) * exp(-alpha * Tij) * (1 - exp(-2 * alpha * vcvMatrix))
  return(vcvRescaled) 
}
    	

require(optimx)
data(geospiza)
attach(geospiza)

td<-treedata(geospiza.tree,geospiza.data[,1],sort=T)
ntax=length(td$phy$tip.label)
chdata<- geospiza.data[,1]# TIP data
tree<- td$phy# Tree
n<- length(chdata)

#l<-fitContinuous(tree,chdata,model="delta",bounds=list(delta=c(.003,40)))


#Number of starting points
n<-50
#Max bound
M<-700
#Min bound
m<-400
#Start Value
s <- 1.000001e-08
#Lower Bound delta
#dLB <- 0.001
#Upper bound delta
#dUB <- 
#Lower bound beta
#bLB <- 1.000001e-08
#upper bound beta
#bUB <- 20
#jb <- seq(from=dLB, to=bUB, length.out=n)
j <- seq(from=m, to=M, length.out=n)
wellt <- c("spg", "Rcgmin", "Rvmmin", "bobyqa","L-BFGS-B",1,1,1,1,1,1,1)
well <- c("spg", "Rcgmin", "Rvmmin", "bobyqa","L-BFGS-B","nlminb","ucminf","Nelder-Mead","nlm","CG","BFGS","newuoa")
optimx <- vector("list", length(well))
names(optimx)<-well

for (i in c(1:length(optimx))){
      optimx[[i]]<-vector("list",m)
      names(optimx[[i]]) <- j
}

lt <- vector("list",length(optimx))
names(lt) <- well

fitContinuous.run<-function(){
        k<-matrix(j,1,n)
	k<-rbind(k,apply(k,2,function(x){
            fitContinuous(tree,
                          chdata,
                          model="delta",
                          bounds=list(delta=c(s,x),beta=c(s,x)))
          }))
        for (i in c(1:n)){
                optimx[[1]][[i]] <- k[2,][[i]]$Trait1$spg
                optimx[[2]][[i]] <- k[2,][[i]]$Trait1$Rcgmin
                optimx[[3]][[i]] <- k[2,][[i]]$Trait1$Rvmmin
                optimx[[4]][[i]] <- k[2,][[i]]$Trait1$bobyqa
                optimx[[5]][[i]] <- k[2,][[i]]$Trait1$'L-BFGS-B'
                optimx[[6]][[i]] <- k[2,][[i]]$Trait1$nlminb
                optimx[[7]][[i]] <- k[2,][[i]]$Trait1$ucminf
                optimx[[8]][[i]] <- k[2,][[i]]$Trait1$'Nelder-Mead'
                optimx[[9]][[i]] <- k[2,][[i]]$Trait1$nlm
                optimx[[10]][[i]] <- k[2,][[i]]$Trait1$CG
                optimx[[11]][[i]] <- k[2,][[i]]$Trait1$BFGS
                optimx[[12]][[i]] <- k[2,][[i]]$Trait1$newuoa}
            for(j in c(1:length(optimx))){
              if (well[j]==wellt[j]){
                lt[[j]]<-as.data.frame(t(rbind(unlist(k[1,]),as.data.frame(matrix(unlist(lapply(optimx[[j]],function(x) return(c(as.numeric(x$lnl),as.numeric(x$beta.bnd[1]),as.numeric(x$beta.bnd[2]),x$beta,as.numeric(x$delta.bnd[1]),as.numeric(x$delta.bnd[2]),x$delta)))),7,ncol(k))))))
	        colnames(lt[[j]])<-c("I","L","bLB","bUB","b","dLB","dUB","d")}else{
                lt[[j]]<-as.data.frame(t(rbind(unlist(k[1,]),as.data.frame(matrix(unlist(lapply(optimx[[j]],function(x) return(c(as.numeric(x$lnl),NA,NA,x$beta,NA,NA,x$delta)))),7,ncol(k))))))
	        colnames(lt[[j]])<-c("I","L","bLB","bUB","b","dLB","dUB","d")}}
	return(lt)
      }
#This fitContinuous.run does not include beta bounds.
fitContinuous.run<-function(){
        k<-matrix(j,1,n)
	k<-rbind(k,apply(k,2,function(x){
            fitContinuous(tree,
                          chdata,
                          model="delta",
                          bounds=list(delta=c(s,x)))
          }))
        for (i in c(1:n)){
                optimx[[1]][[i]] <- k[2,][[i]]$Trait1$spg
                optimx[[2]][[i]] <- k[2,][[i]]$Trait1$Rcgmin
                optimx[[3]][[i]] <- k[2,][[i]]$Trait1$Rvmmin
                optimx[[4]][[i]] <- k[2,][[i]]$Trait1$bobyqa
                optimx[[5]][[i]] <- k[2,][[i]]$Trait1$'L-BFGS-B'
                optimx[[6]][[i]] <- k[2,][[i]]$Trait1$nlminb
                optimx[[7]][[i]] <- k[2,][[i]]$Trait1$ucminf
                optimx[[8]][[i]] <- k[2,][[i]]$Trait1$'Nelder-Mead'
                optimx[[9]][[i]] <- k[2,][[i]]$Trait1$nlm
                optimx[[10]][[i]] <- k[2,][[i]]$Trait1$CG
                optimx[[11]][[i]] <- k[2,][[i]]$Trait1$BFGS
                optimx[[12]][[i]] <- k[2,][[i]]$Trait1$newuoa}
            for(j in c(1:length(optimx))){
              if (well[j]==wellt[j]){
                lt[[j]]<-as.data.frame(t(rbind(unlist(k[1,]),as.data.frame(matrix(unlist(lapply(optimx[[j]],function(x) return(c(as.numeric(x$lnl),as.numeric(x$beta.bnd[1]),as.numeric(x$beta.bnd[2]),x$beta,as.numeric(x$delta.bnd[1]),as.numeric(x$delta.bnd[2]),x$delta)))),7,ncol(k))))))}else{
                lt[[j]]<-as.data.frame(t(rbind(unlist(k[1,]),as.data.frame(matrix(unlist(lapply(optimx[[j]],function(x) return(c(as.numeric(x$lnl),NA,NA,x$beta,NA,NA,x$delta)))),7,ncol(k))))))
	        colnames(lt[[j]])<-c("I","L","bLB","bUB","b","dLB","dUB","d")}}
	return(lt)
      }

begin.time <-proc.time()
l2<-fitContinuous.run()

#Time in minutes
total.time <- as.numeric(proc.time()[3]-begin.time[3])/(60)

save.image("/home/michels/repository/phylooptim/pkg/R/geiger/geigererror.RData")
