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
conv <- vector("list",length(optimx))
for (i in c(1:length(optimx))){conv[[i]] <- rep(NA,m)}

k<-matrix(jj,1,m)
k<-rbind(k,apply(k,2,function(x){ouch.wrap.x(tree,trait,model=c("ou1"),x)}))

for (i in c(1:length(optimx))){for (j in c(1:m)){optimx[[i]][[j]] <- k[2,][[j]][[i]][[1]];time[[i]][j] <- k[2,][[j]][[i]][[2]];conv[[i]][j] <- k[2,][[j]][[i]][[3]]}}

for(j in c(1:length(optimx))){
  if (well[j]==wellt[j]){
    lt[[j]]<-as.data.frame(t(rbind(unlist(k[1,]),as.data.frame(matrix(unlist(lapply(optimx[[j]],function(x) return(c(0.001,1e50,as.numeric(x@theta),x@loglik,x@sigma,x@sqrt.alpha)))),6,ncol(k))),conv[[j]],time[[j]])))
    colnames(lt[[j]])<-c("I","lb","ub","T","L","S","A","conv","time")}else{
    lt[[j]]<-as.data.frame(t(rbind(unlist(k[1,]),as.data.frame(matrix(unlist(lapply(optimx[[j]],function(x) return(c(NA,NA,as.numeric(x@theta),x@loglik,x@sigma,x@sqrt.alpha)))),6,ncol(k))),conv[[j]],time[[j]])))
    colnames(lt[[j]])<-c("I","lb","ub","T","L","S","A","conv","time")}
}

return(lt)
}
