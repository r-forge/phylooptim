# 
# 
# Author: Kurt Michels   kamichels@iplantcollaborative.org
#
###############################################################################

##DATA INPUT/GENERATION
#If one wants to load a dataset
#can be modified to load other files

## ace.R (2011-07-18)

##   Ancestral Character Estimation

## Copyright 2005-2011 Emmanuel Paradis and Ben Bolker

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

ace <- function(x, phy, type = "continuous", method = "ML", CI = TRUE,
                model = if (type == "continuous") "BM" else "ER",
                scaled = TRUE, kappa = 1, corStruct = NULL, ip = 0.1)
{
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
        else warning("the names of 'x' and the tip labels of the tree do not match: the former were ignored in the analysis.")
    }

wellt <- c("spg", "Rcgmin", "Rvmmin", "bobyqa","L-BFGS-B",1,1,1,1,1,1,1)
well <- c("spg", "Rcgmin", "Rvmmin", "bobyqa","L-BFGS-B","nlminb","ucminf","Nelder-Mead","nlm","CG","BFGS","newuoa")
obj.list <- vector("list", length(well))
names(obj.list)<-well

  for (i in c(1:length(well))){
    obj <- list()
    if (kappa != 1) phy$edge.length <- phy$edge.length^kappa
    if (type == "continuous") {
        switch(method, "REML" = {
            minusLogLik <- function(sig2) {
                if (sig2 < 0) return(1e100)
                V <- sig2 * vcv(phy)
                ## next three lines borrowed from dmvnorm() in 'mvtnorm'
                distval <- mahalanobis(x, center = mu, cov = V)
                logdet <- sum(log(eigen(V, symmetric = TRUE, only.values = TRUE)$values))
                (nb.tip * log(2 * pi) + logdet + distval)/2
            }
            mu <- rep(ace(x, phy, method="pic")$ace[1], nb.tip)
            out <- nlm(minusLogLik, 1, hessian = TRUE)

    out1 <- optimx(minusLogLik, p=1,
                  lower = 0.0001, upper = 1e50,
                  control=list(all.methods=TRUE, trace=0))

    out1<-out1[sort.list(unlist(out1[,3])), ]

            sigma2 <- out1[,1][[i]][-1]
            se_sgi2 <- sqrt(1/out$hessian)
            out <- nlm(function(p) minus.REML.BM(p),
                       p = rep(mu[1], nb.node), hessian = TRUE)

    out1 <- optimx(function(p) minus.REML.BM(p), p=rep(mu[1], nb.node),
                  lower = 0.0001, upper = 1e50,
                  control=list(all.methods=TRUE, trace=0))

    out1<-out1[sort.list(unlist(out1[,3])), ]

            obj$ace <- out1[,1][[i]][-1]
            obj$resloglik <- -out$minimum
            names(obj$ace) <- (nb.tip + 1):(nb.tip + nb.node)
            obj$sigma2 <- c(sigma2, se_sgi2)
            if (CI) {
                se <- sqrt(diag(solve(out$hessian)))
                tmp <- se * qt(0.025, nb.node)
                obj$CI95 <- cbind(obj$ace + tmp, obj$ace - tmp)
            }
        }, "pic" = {
            if (model != "BM")
                stop('the "pic" method can be used only with model = "BM".')
            ## See pic.R for some annotations.
            phy <- reorder(phy, "pruningwise")
            phenotype <- numeric(nb.tip + nb.node)
            phenotype[1:nb.tip] <- if (is.null(names(x))) x else x[phy$tip.label]
            contr <- var.con <- numeric(nb.node)
            ans <- .C("pic", as.integer(nb.tip), as.integer(nb.node),
                      as.integer(phy$edge[, 1]), as.integer(phy$edge[, 2]),
                      as.double(phy$edge.length), as.double(phenotype),
                      as.double(contr), as.double(var.con),
                      as.integer(CI), as.integer(scaled),
                      PACKAGE = "ape")
            obj$ace <- ans[[6]][-(1:nb.tip)]
            names(obj$ace) <- (nb.tip + 1):(nb.tip + nb.node)
            if (CI) {
                se <- sqrt(ans[[8]])
                tmp <- se * qnorm(0.025)
                obj$CI95 <- cbind(obj$ace + tmp, obj$ace - tmp)
            }
        }, "ML" = {
            if (model == "BM") {
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
                out <- nlm(function(p) dev.BM(p),
                           p = c(1, rep(mean(x), nb.node)), hessian = TRUE)

        out1 <- optimx(function(p) dev.BM(p), p=c(1, rep(mean(x), nb.node)),
                      lower = 0.0001, upper = 1e50,
                      control=list(all.methods=TRUE, trace=0))

        out1<-out1[sort.list(unlist(out1[,3])), ]

                obj$ace <- out1[,1][[i]][-1]
                obj$loglik <- -out$minimum / 2
                names(obj$ace) <- (nb.tip + 1):(nb.tip + nb.node)
                se <- sqrt(diag(solve(out$hessian)))
                obj$sigma2 <- c(out$estimate[1], se[1])
                if (CI) {
                    tmp <- se[-1] * qt(0.025, nb.node)
                    obj$CI95 <- cbind(obj$ace + tmp, obj$ace - tmp)
                }
            }
        }, "GLS" = {
            if (is.null(corStruct))
                stop('you must give a correlation structure if method = "GLS".')
            if (class(corStruct)[1] == "corMartins")
                M <- corStruct[1] * dist.nodes(phy)
            if (class(corStruct)[1] == "corGrafen")
                phy <- compute.brlen(attr(corStruct, "tree"),
                                     method = "Grafen",
                                     power = exp(corStruct[1]))
            if (class(corStruct)[1] %in% c("corBrownian", "corGrafen")) {
                dis <- dist.nodes(attr(corStruct, "tree"))
                MRCA <- mrca(attr(corStruct, "tree"), full = TRUE)
                M <- dis[as.character(nb.tip + 1), MRCA]
                dim(M) <- rep(sqrt(length(M)), 2)
            }
            varAY <- M[-(1:nb.tip), 1:nb.tip]
            varA <- M[-(1:nb.tip), -(1:nb.tip)]
            V <- corMatrix(Initialize(corStruct, data.frame(x)),
                           corr = FALSE)
            invV <- solve(V)
            o <- gls(x ~ 1, data.frame(x), correlation = corStruct)
            GM <- o$coefficients
            obj$ace <- drop(varAY %*% invV %*% (x - GM) + GM)
            names(obj$ace) <- (nb.tip + 1):(nb.tip + nb.node)
            if (CI) {
                se <- sqrt((varA - varAY %*% invV %*% t(varAY))[cbind(1:nb.node, 1:nb.node)])
                tmp <- se * qnorm(0.025)
                obj$CI95 <- cbind(obj$ace + tmp, obj$ace - tmp)
            }
        })
    } else { # type == "discrete"
        if (method != "ML")
            stop("only ML estimation is possible for discrete characters.")
        if (!is.factor(x)) x <- factor(x)
        nl <- nlevels(x)
        lvls <- levels(x)
        x <- as.integer(x)
        if (is.character(model)) {
            rate <- matrix(NA, nl, nl)
            if (model == "ER") np <- rate[] <- 1
            if (model == "ARD") {
                np <- nl*(nl - 1)
                rate[col(rate) != row(rate)] <- 1:np
            }
            if (model == "SYM") {
                np <- nl * (nl - 1)/2
                sel <- col(rate) < row(rate)
                rate[sel] <- 1:np
                rate <- t(rate)
                rate[sel] <- 1:np
            }
        } else {
            if (ncol(model) != nrow(model))
              stop("the matrix given as 'model' is not square")
            if (ncol(model) != nl)
              stop("the matrix 'model' must have as many rows as the number of categories in 'x'")
            rate <- model
            np <- max(rate)
        }
        index.matrix <- rate
        tmp <- cbind(1:nl, 1:nl)
        index.matrix[tmp] <- NA
        rate[tmp] <- 0
        rate[rate == 0] <- np + 1 # to avoid 0's since we will use this as numeric indexing

        liks <- matrix(0, nb.tip + nb.node, nl)
        TIPS <- 1:nb.tip
        liks[cbind(TIPS, x)] <- 1
        phy <- reorder(phy, "pruningwise")

        Q <- matrix(0, nl, nl)
        dev <- function(p, output.liks = FALSE) {
            if (any(is.nan(p)) || any(is.infinite(p))) return(1e50)
            ## from Rich FitzJohn:
            comp <- numeric(nb.tip + nb.node) # Storage...
            Q[] <- c(p, 0)[rate]
            diag(Q) <- -rowSums(Q)
            for (i  in seq(from = 1, by = 2, length.out = nb.node)) {
                j <- i + 1L
                anc <- phy$edge[i, 1]
                des1 <- phy$edge[i, 2]
                des2 <- phy$edge[j, 2]
                v.l <- matexpo(Q * phy$edge.length[i]) %*% liks[des1, ]
                v.r <- matexpo(Q * phy$edge.length[j]) %*% liks[des2, ]
                v <- v.l * v.r
                comp[anc] <- sum(v)
                liks[anc, ] <- v/comp[anc]
            }
            if (output.liks) return(liks[-TIPS, ])
            dev <- -2 * sum(log(comp[-TIPS]))
            if (is.na(dev)) Inf else dev
        }

if (well[i]==wellt[i]){

        out <- optimx(function(p) dev(p), p=rep(ip, length.out = np),
                      lower = rep(0, np), upper = rep(1e50, np),
                      method=well[i])
        lb <- 0;ub <- 1e50}else{
        
        out <- optimx(function(p) dev(p), p=rep(ip, length.out = np),
                      method=well[i])}

        obj$lb <- lb
        obj$ub <- ub
        obj$loglik <- -out$fvalues$fvalues/2
        obj$rates <- out$par$par
        oldwarn <- options("warn")
        options(warn = -1)

    #    h <- nlm(function(p) dev(p), p = obj$rates, iterlim = 1,
    #             stepmax = 0, hessian = TRUE)$hessian
        options(oldwarn)
    #    if (any(h == 0))
    #      warning("The likelihood gradient seems flat in at least one dimension (gradient null):\ncannot compute the standard-errors of the transition rates.\n")
    #    else obj$se <- sqrt(diag(solve(h)))
   #     obj$index.matrix <- index.matrix
        if (CI) {
            obj$lik.anc <- dev(obj$rates, TRUE)
            colnames(obj$lik.anc) <- lvls
        }
    }
    obj$call <- match.call()
    class(obj) <- "ace"
    obj.list[[i]]<-obj
}
obj.list
}

logLik.ace <- function(object, ...) object$loglik

deviance.ace <- function(object, ...) -2*object$loglik

AIC.ace <- function(object, ..., k = 2)
{
    if (is.null(object$loglik)) return(NULL)
    ## Trivial test of "type"; may need to be improved
    ## if other models are included in ace(type = "c")
    np <- if (!is.null(object$sigma2)) 1 else length(object$rates)
    -2*object$loglik + np*k
}

### by BB:
anova.ace <- function(object, ...)
{
    X <- c(list(object), list(...))
    df <- sapply(lapply(X, "[[", "rates"), length)
    ll <- sapply(X, "[[", "loglik")
    ## check if models are in correct order?
    dev <- c(NA, 2*diff(ll))
    ddf <- c(NA, diff(df))
    table <- data.frame(ll, df, ddf, dev,
                        pchisq(dev, ddf, lower.tail = FALSE))
    dimnames(table) <- list(1:length(X), c("Log lik.", "Df",
                                           "Df change", "Resid. Dev",
                                           "Pr(>|Chi|)"))
    structure(table, heading = "Likelihood Ratio Test Table",
              class = c("anova", "data.frame"))
}

print.ace <- function(x, digits = 4, ...)
{
    cat("\n    Ancestral Character Estimation\n\n")
    cat("Call: ")
    print(x$call)
    cat("\n")
    if (!is.null(x$loglik))
        cat("    Log-likelihood:", x$loglik, "\n\n")
    if (!is.null(x$resloglik))
        cat("    Residual log-likelihood:", x$resloglik, "\n\n")
    ratemat <- x$index.matrix
    if (is.null(ratemat)) { # to be improved
        class(x) <- NULL
        x$resloglik <- x$loglik <- x$call <- NULL
        print(x)
    } else {
        dimnames(ratemat)[1:2] <- dimnames(x$lik.anc)[2]
        cat("Rate index matrix:\n")
        print(ratemat, na.print = ".")
        cat("\n")
        npar <- length(x$rates)
        estim <- data.frame(1:npar, round(x$rates, digits), round(x$se, digits))
        cat("Parameter estimates:\n")
        names(estim) <- c("rate index", "estimate", "std-err")
        print(estim, row.names = FALSE)
        cat("\nScaled likelihoods at the root (type '...$lik.anc' to get them for all nodes):\n")
        print(x$lik.anc[1, ])
    }
}

library(optimx)

##For geospiza data
#library(geiger)
#data(geospiza)
#phy <- geospiza$geospiza.tree
#dv <- as.factor(geospiza$geospiza.data[,1]>4.2)
#names(dv) <- names(geospiza$geospiza.data)

#For aquilegia data
require(ape)
phy <- read.tree("Aquilegia.phy")
dv <- read.delim("Aquilegia-traits.txt", header=T)[,2]
names(dv) <- read.delim("Aquilegia-traits.txt", header=T)[,1]
l1 <- ace(dv,phy,type='discrete',ip=3)

#Number of starting points
m<- 5
j<-seq(from=1/100, to=10, length.out=m)
#j<-seq(from=0.72, to=4.31, length.out=m)
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

ace.run<-function(){
        k<-matrix(j,1,m)
	k<-rbind(k,apply(k,2,function(x) {ace(dv,phy,type='discrete',ip=x)}))
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
                lt[[j]]<-as.data.frame(t(rbind(unlist(k[1,]),as.data.frame(matrix(unlist(lapply(optimx[[j]],function(x) return(c(x$lb,x$ub,x$rates,x$loglik)))),4,ncol(k))))))
              	colnames(lt[[j]])<-c("I","lb","ub","P","L")}else{
                lt[[j]]<-as.data.frame(t(rbind(unlist(k[1,]),as.data.frame(matrix(unlist(lapply(optimx[[j]],function(x) return(c(NA,NA,x$rates,x$loglik,x$se)))),4,ncol(k))))))                  
	        colnames(lt[[j]])<-c("I","lb","ub","P","L")}}
	return(lt)
      }
        

begin.time <-proc.time()
l<-ace.run()

#Time in minutes
total.time <- as.numeric(proc.time()[3]-begin.time[3])/(60)

save.image("/home/michels/repository/phylooptim/pkg/R/ace/aceerror.RData")
