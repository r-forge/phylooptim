f.geig.disc <- function(x)
{

yy <- x

m <- length(yy)
jj <- yy

optimx <- vector("list", length(well))
names(optimx)<-well

for (i in c(1:length(optimx))){
      optimx[[i]]<-vector("list",m)
}

lt <- vector("list",length(optimx))
names(lt) <- well

        k<-matrix(jj,1,m)
 k<-rbind(k,apply(k,2,function(x) {fitDiscrete(dv$phy,dv$data,model="ER",treeTransform="delta",s=x)}))
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
                lt[[j]]<-as.data.frame(t(rbind(as.data.frame(matrix(unlist(lapply(optimx[[j]],function(x) return(c(x$Trait1$lb[1],x$Trait1$ub[1],x$Trait1$start[1],x$Trait1$lb[2],x$Trait1$ub[2],x$Trait1$start[2],x$Trait1$lnl,x$Trait1$time,x$Trait1$conv)))),9,ncol(k))))))
               colnames(lt[[j]])<-c("lb1","ub1","S1","lb2","ub2","S2","L","time","conv")}
 return(lt)
}
