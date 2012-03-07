getRateCats<-function(data, model)
{
	if(model=="ER") return(1)
	n<-nlevels(factor(data))
	if(model=="SYM") return(n*(n-1)/2)
	if(model=="ARD") return(n*(n-1))
	
}
	
if(treeTransform=="none") {
	f<-function(x) {
			likelihoodDiscrete(dv$phy, dv$data[,1], exp(x), model)
			}
	nep=0; pLow=-10; pHigh=log(1); pStart=NULL;		
}	

if(treeTransform=="lambda") {
	f<-function(x) {
			likelihoodDiscrete(dv$phy, dv$data[,1], exp(x[-1]), lambda=exp(x[1]), model=model)
			}	
	nep=1; pLow=-10; pHigh=log(1); pStart=0.1;
}

if(treeTransform=="delta") {
	f<-function(x) {
			likelihoodDiscrete(dv$phy, dv$data[,1], exp(x[-1]), delta=exp(x[1]), model=model)
			}
	nep=1; pLow=-10; pHigh=log(10);pStart=0.1;
}

if(treeTransform=="kappa") {
	f<-function(x) {
			likelihoodDiscrete(dv$phy, dv$data[,1], exp(x[-1]), kappa=exp(x[1]), model=model)
			}
	nep=1; pLow=-10; pHigh=log(1);pStart=0.1;
}

if(treeTransform=="linearChange") {
	f<-function(x) {
			likelihoodDiscrete(dv$phy, dv$data[,1], exp(x[-1]), endRate=exp(x[1]), linear=T, model=model)
			}
	nep=1; pLow=-10; pHigh=log(10);pStart=0.1;
}

if(treeTransform=="exponentialChange") {
	f<-function(x) {
			likelihoodDiscrete(dv$phy, dv$data[,1], exp(x[-1]), endRate=exp(x[1]), model=model)
			}
	nep=1; pLow=-10; pHigh=log(10);pStart=0.1;
}

if(treeTransform=="twoRate") {
	f<-function(x) {
			likelihoodDiscrete(dv$phy, dv$data[,1], exp(x[-(1:2)]), breakPoint=x[1], endRate=exp(x[2]), model=model)
			}
	mv<-max(branching.times(dv$phy))	
	nep=2; pLow=c(mv/1000, -10); pHigh=c(mv, 10);pStart=c(0.1, 0.1);
}
				
						
nRateCats<-getRateCats(dv$data[,1], model)

outTries <- vector("list", 20)
for (j in c(1:20)){outTries[[j]] <- list()}
totalbl<-sum(dv$phy$edge.length)
minQ=log(0.01/totalbl)
maxQ=log(1000/totalbl)
ntries<-20
ltry<-numeric(ntries)
ltry[]<-NA
lsol<-matrix(nrow= ntries, ncol=nRateCats+nep)
sp<-numeric(nRateCats)
qTries<-exp(-5:4)
	
if(nep==0) {
            lower=rep(minQ, nRateCats)
	    upper=rep(maxQ, nRateCats)
	   } else {
	    lower=c(pLow, rep(minQ, nRateCats))
	    upper=c(pHigh, rep(maxQ, nRateCats))
	   }
