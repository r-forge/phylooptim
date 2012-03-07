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

        k<-matrix(jj,1,m)
 k<-rbind(k,apply(k,2,function(x) {ace(as.vector(dv$data),dv$phy,type='continuous',ip=x,CI=FALSE)}))
            for (i in c(1:m)){
                optimx[[1]][[i]] <- k[2,][[i]]$spg
                optimx[[2]][[i]] <- k[2,][[i]]$Rcgmin
                optimx[[3]][[i]] <- k[2,][[i]]$Rvmmin
                optimx[[4]][[i]] <- k[2,][[i]]$bobyqa
                optimx[[5]][[i]] <- k[2,][[i]]$'L-BFGS-B'
                optimx[[6]][[i]] <- k[2,][[i]]$nlminb
                optimx[[7]][[i]] <- k[2,][[i]]$ucminf
                optimx[[8]][[i]] <- k[2,][[i]]$'Nelder-Mead'
                optimx[[9]][[i]] <- k[2,][[i]]$nlm
                optimx[[10]][[i]] <- k[2,][[i]]$CG
                optimx[[11]][[i]] <- k[2,][[i]]$BFGS
                optimx[[12]][[i]] <- k[2,][[i]]$newuoa}
            for(j in c(1:length(optimx))){
              if (well[j]==wellt[j]){
                lt[[j]]<-as.data.frame(t(rbind(unlist(k[1,]),as.data.frame(matrix(unlist(lapply(optimx[[j]],function(x) return(c(x$lb,x$ub,x$sigma2[1],x$loglik,x$time,x$conv)))),6,ncol(k))))))
               colnames(lt[[j]])<-c("I","lb","ub","S","L","time","conv")}else{
                lt[[j]]<-as.data.frame(t(rbind(unlist(k[1,]),as.data.frame(matrix(unlist(lapply(optimx[[j]],function(x) return(c(NA,NA,x$sigma2[1],x$loglik,x$time,x$conv)))),6,ncol(k))))))                  
         colnames(lt[[j]])<-c("I","lb","ub","S","L","time","conv")}}
 return(lt)
}
