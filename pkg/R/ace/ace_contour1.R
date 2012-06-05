rm(list=ls())

name_all <- c("geo","aqui","mono","mam")
well <- c("spg", "Rcgmin", "Rvmmin", "bobyqa","L-BFGS-B","nlminb","ucminf","Nelder-Mead","nlm","CG","BFGS","newuoa")

#for (j in c(1:length(name_all))){
  j <- 1
  name <- name_all[[j]]
  load(paste("/home/michels/Hallowed/repository/phylooptim/pkg/R/ace/",name,"acects.RData",sep=""))
  xx <- list()
  yy <- list()
    for (i in c(1:length(well))){
      xx[[i]] <- l[[i]]$foostore$sigma2
      yy[[i]] <- l[[i]]$foostore$ace
    }
#     lx <- max(xx)-min(xx)
#     ly <- max(yy)-min(yy)
      x1range <- seq(0,1.2,length=100)
      x2range <- seq(3.7,4.7,length=100)
#     x1range <- seq(min(xx)+0.0000001,max(xx),length=100)
#     x2range <- seq(min(yy),max(yy),length=100)

      x <- as.vector(dv$data)
      phy <- dv$phy
      type <- 'continuous'
      numcol<-50 #number of colors to use
      if (!inherits(phy, "phylo"))
          stop('object "phy" is not of class "phylo".')
      if (is.null(phy$edge.length))
          stop("tree has no branch lengths")
      type <- match.arg(type, c("continuous", "discrete"))
      nb.tip <- length(phy$tip.label)
      nb.node <- phy$Nnode
      if (nb.node != nb.tip - 1)
          stop('"phy" is not rooted AND fully dichotomous.')
      if (length(x) != nb.tip)
          stop("length of phenotypic and of phylogenetic data do not match.")
      if (!is.null(names(x))) {
          if(all(names(x) %in% phy$tip.label))
            x <- x[phy$tip.label]
          else warning("the names of 'x' and the tip labels of the tree do not match: the former wfere ignored in the analysis.")
      }
      tip <- phy$edge[, 2] <= nb.tip
      dev.BM <- function(p) {
        if (p[1] < 0) return(1e100) # in case sigma² is negative
        x1 <- p[-1][phy$edge[, 1] - nb.tip]
        x2 <- numeric(length(x1))
        x2[tip] <- x[phy$edge[tip, 2]]
        x2[!tip] <- p[-1][phy$edge[!tip, 2] - nb.tip]
        -2 * (-sum((x1 - x2)^2/phy$edge.length)/(2*p[1]) -
        nb.node * log(p[1]))
      }

      res<-t(sapply(x1range,function(x,y=x2range) sapply(y,function(z) dev.BM(c(x,rep(z, length.out = nb.node))) )))

      if(any(res[!is.na(res)]<0)) {
      	res[!is.na(res) & res>=0 & res<1]<-0
      	res[!is.na(res) & res<0]<--1*log(-1*res[!is.na(res) & res<0])
      }
      res[!is.na(res) & res>0]<-log(res[!is.na(res) & res>0])

      s <- list()
      for (i in c(1:length(well))){s[[i]] <- seq(length(xx[[i]])-1)}

      nlevels=25
      levels=pretty(range(res,na.rm=TRUE),nlevels,finite=TRUE)
      k<-pretty(trunc(range(res,na.rm=T)),nlevels)
      wr.pal<-colorRampPalette(c("grey","black"))
      wr<-wr.pal(round(numcol*sum(k>=0)/length(k)))
      bw.pal<-colorRampPalette(c("white","grey"))
      bw<-bw.pal(round(numcol*sum(k<=0)/length(k)))
      cols<-c(bw,wr)
      rcols<-cols[round((1:(length(levels)))*length(cols)/length(levels))]
      
  xxsig <- list()
  yyace <- list()
    for (i in c(1:length(well))){
      kk <- floor(length(l[[i]]$foostore$loglik)/10)     
      newseq <- seq(0,kk*10,by=10)
      xxsig[[i]] <- rep(NA,length=kk)
      yyace[[i]] <- rep(NA,length=kk)
      for (k in (1:kk)){
          lnLTries <- na.omit(l[[i]]$foostore$loglik[(newseq[k]+1):newseq[k+1]])
          if (sum(attr(lnLTries,"na.action"))<sum(1:10)){
          sigma2x <- na.omit(l[[i]]$foostore$sigma2[(newseq[k]+1):newseq[k+1]])
          acey <- na.omit(l[[i]]$foostore$ace[(newseq[k]+1):newseq[k+1]])
          bestIndex<-which.min(lnLTries[1:20])
          xxsig[[i]][k] <- sigma2x[bestIndex]
          yyace[[i]][k] <- acey[bestIndex]}else{next}
      }
    }
      
      png(file="contourtotal.png")
      filled.contour(x1range,x2range,res,levels=levels,col=rcols,  
      		plot.axes={
      			axis(1)
      			axis(2)
            plot.title=title(main="All Optimizers",xlab="sigma2", ylab="ace")
      			contour(x1range,x2range,res,add=TRUE,nlevels=length(levels))
                        
     	                arrows(xxsig[[1]][s[[1]]], yyace[[1]][s[[1]]], xxsig[[1]][s[[1]]+1],
                              yyace[[1]][s[[1]]+1],length=.1,col=4,
                              type = "l",lwd=2)
                       points(l[[1]]$sigma2,l[[1]]$ace[1],col=2,pch=15,cex=3)
                       points(xx[[1]][1],yy[[1]][1],col="green",pch=19,cex=2)
                        
                        arrows(xxsig[[2]][s[[2]]], yyace[[2]][s[[2]]], xxsig[[2]][s[[2]]+1],
                               yyace[[2]][s[[2]]+1],length=.1,col=5,
                               type = "l",lwd=2)
                        points(l[[2]]$sigma2,l[[2]]$ace[1],col=2,pch=15,cex=3)
                        points(xx[[2]][1],yy[[2]][1],col="green",pch=19,cex=2)
                        
                        arrows(xxsig[[3]][s[[3]]], yyace[[3]][s[[3]]], xxsig[[3]][s[[3]]+1],
                               yyace[[3]][s[[3]]+1],length=.1,col=1,
                               type = "l",lwd=2)
                        points(l[[3]]$sigma2,l[[3]]$ace[1],col=2,pch=15,cex=3)
                        points(xx[[3]][1],yy[[3]][1],col="green",pch=19,cex=2)
                        
                        arrows(xxsig[[4]][s[[4]]], yyace[[4]][s[[4]]], xxsig[[4]][s[[4]]+1],
                               yyace[[4]][s[[4]]+1],length=.1,col=6,
                               type = "l",lwd=2)
                        points(l[[4]]$sigma2,l[[4]]$ace[1],col=2,pch=15,cex=3)
                        points(xx[[4]][1],yy[[4]][1],col="green",pch=19,cex=2)
                        
                        arrows(xxsig[[5]][s[[5]]], yyace[[5]][s[[5]]], xxsig[[5]][s[[5]]+1],
                               yyace[[5]][s[[5]]+1],length=.1,col=1,
                               type = "l",lwd=2)
                        points(l[[5]]$sigma2,l[[5]]$ace[1],col=2,pch=15,cex=3)
                        points(xx[[5]][1],yy[[5]][1],col="green",pch=19,cex=2)
                        
                        arrows(xxsig[[6]][s[[6]]], yyace[[6]][s[[6]]], xxsig[[6]][s[[6]]+1],
                               yyace[[6]][s[[6]]+1],length=.1,col=8,
                               type = "l",lwd=2)
                        points(l[[6]]$sigma2,l[[6]]$ace[1],col=2,pch=15,cex=3)
                        points(xx[[6]][1],yy[[6]][1],col="green",pch=19,cex=2)
                        
                        arrows(xxsig[[7]][s[[7]]], yyace[[7]][s[[7]]], xxsig[[7]][s[[7]]+1],
                               yyace[[7]][s[[7]]+1],length=.1,col=8,
                               type = "l",lwd=2)
     	                points(l[[7]]$sigma2,l[[7]]$ace[1],col=2,pch=15,cex=3)
                        points(xx[[7]][1],yy[[7]][1],col="green",pch=19,cex=2)
                        
                        arrows(xxsig[[8]][s[[8]]], yyace[[8]][s[[8]]], xxsig[[8]][s[[8]]+1],
                               yyace[[8]][s[[8]]+1],length=.1,col=7,
                               type = "l",lwd=2)
                        points(l[[8]]$sigma2,l[[8]]$ace[1],col=2,pch=15,cex=3)
                        points(xx[[8]][1],yy[[8]][1],col="green",pch=19,cex=2)
                        
                        arrows(xxsig[[9]][s[[9]]], yyace[[9]][s[[9]]], xxsig[[9]][s[[9]]+1],
                               yyace[[9]][s[[9]]+1],length=.1,col=6,
                               type = "l",lwd=2)
                        points(l[[9]]$sigma2,l[[9]]$ace[1],col=2,pch=15,cex=3)
                        points(xx[[9]][1],yy[[9]][1],col="green",pch=19,cex=2)
                        
                        arrows(xxsig[[10]][s[[10]]],yyace[[10]][s[[10]]],xxsig[[10]][s[[10]]+1],
                               yyace[[10]][s[[10]]+1],length=.1,col=5,
                               type = "l",lwd=2)
                        points(l[[10]]$sigma2,l[[10]]$ace[1],col=2,pch=15,cex=3)
                        points(xx[[10]][1],yy[[10]][1],col="green",pch=19,cex=2)
                        
                        arrows(xxsig[[11]][s[[11]]],yyace[[11]][s[[11]]],xxsig[[11]][s[[11]]+1],
                               yyace[[11]][s[[11]]+1],length=.1,col=8,
                               type = "l",lwd=2)
                        points(l[[11]]$sigma2,l[[11]]$ace[1],col=2,pch=15,cex=3)
                        points(xx[[11]][1],yy[[11]][1],col="green",pch=19,cex=2)                        
                        arrows(xxsig[[12]][s[[12]]],yyace[[12]][s[[12]]],xxsig[[12]][s[[12]]+1],
                               yyace[[12]][s[[12]]+1],length=.1,col=5,
                               type = "l",lwd=2)
                        points(l[[12]]$sigma2,l[[12]]$ace[1],col=2,pch=15,cex=3)
                        points(xx[[12]][1],yy[[12]][1],col="green",pch=19,cex=2)
      		}
      )
      dev.off()
    }
#rm(list=ls()[-1*c(which(ls()=="well"),which(ls()=="name_all"),which(ls()=="j"))])
#}
