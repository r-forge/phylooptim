f <- function(x)
{ 

yy <- x

m <- length(yy[[1]])
jj <- yy

optimx <- vector("list", length(well))
names(optimx)<-well

for (i in c(1:length(optimx))){
     optimx[[i]]<-vector("list",m)
}

lt <- vector("list",length(optimx))
names(lt) <- well
for (i in c(1:length(lt))){lt[[i]] <- list()}
for (i in c(1:length(lt))){
     lt[[i]][[2]]<-vector("list",m)
}
for (i in c(1:length(lt))){
  for (j in c(1:m)){lt[[i]][[2]][[j]] <- vector("list",1)}
}

time <- vector("list",length(optimx))
for (i in c(1:length(optimx))){time[[i]] <- rep(NA,m)}
conv <- vector("list",length(optimx))
for (i in c(1:length(optimx))){conv[[i]] <- rep(NA,m)}
path <- vector("list",length(optimx))
for (i in c(1:length(optimx))){path[[i]] <- list()}

qqww <- vector("list",m)
for (i in c(1:length(qqww))){qqww[[i]] <- list()}
h <- 1
for (i in c(1:dim(jj[[1]])[1])){
  for (j in c(1:dim(jj[[1]])[2])){
    qqww[[h]][[1]] <- c(jj[[2]][i,j],jj[[1]][i,j])
    qqww[[h]][[2]] <- ouch.wrap2.x(tree,trait,model=c("ou1"),itn=iter,y=c(jj[[2]][i,j],jj[[1]][i,j]))
    h <- h+1
  }
}
sv <- NULL
for (i in c(1:m)){sv <- rbind(sv,qqww[[i]][[1]])}

for (i in c(1:length(optimx))){for (j in c(1:m)){optimx[[i]][[j]] <- qqww[[j]][[2]][[i]][[1]];time[[i]][j] <- qqww[[j]][[2]][[i]][[3]];conv[[i]][j] <- qqww[[j]][[2]][[i]][[4]];path[[i]][[j]] <- qqww[[j]][[2]][[i]][[2]]}}

for(j in c(1:length(optimx))){
  if (well[j]==wellt[j]){
    lt[[j]][[1]]<-as.data.frame(cbind(sv,t(as.data.frame(matrix(unlist(lapply(optimx[[j]],function(x) return(c(0.001,1e50,x@loglik,as.numeric(x@sqrt.alpha),as.numeric(x@sigma))))),5,m))),conv[[j]],time[[j]]))
    colnames(lt[[j]][[1]])<-c("sa.I","sig.I","lb","ub","L","SA","SIG","conv","time")
  for (i in c(1:m)){lt[[j]][[2]][[i]][[1]] <- path[[j]][[i]]}
  }else{
     lt[[j]][[1]]<-as.data.frame(cbind(sv,t(as.data.frame(matrix(unlist(lapply(optimx[[j]],function(x) return(c(NA,NA,x@loglik,as.numeric(x@sqrt.alpha),as.numeric(x@sigma))))),5,m))),conv[[j]],time[[j]]))
    colnames(lt[[j]][[1]])<-c("sa.I","sig.I","lb","ub","L","SA","SIG","conv","time")   
  for (i in c(1:m)){lt[[j]][[2]][[i]][[1]] <- path[[j]][[i]]}
   }
}

for (i in c(1:length(lt))){
     names(lt[[i]][[2]]) <- lt[[i]][[1]][,1]
}
for (i in c(1:length(lt))){
  for (j in c(1:m)){names(lt[[i]][[2]][[j]]) <- lt[[i]][[1]][,2][j]}
}

return(lt)
}
