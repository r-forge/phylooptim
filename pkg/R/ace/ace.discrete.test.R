require(geiger)

data(geospiza)
name <- c("geo")
tree <- geospiza$geospiza.tree
a.trait <- data.frame(geospiza$geospiza.data,T=as.factor(geospiza$geospiza.data[,1]>4.2))

##Which column of data do you want?
kk <- 6

if (dim(a.trait)[2] == 1){a.trait <- data.frame(a.trait,well=rep(1,length(a.trait)))}
nc <- name.check(tree,a.trait)
if (nc[1]=="OK"){nc$Tree.not.data <- NULL}
tree <- drop.tip(tree,nc$Tree.not.data)
tree$edge.length[tree$edge.length<1e-5]=1e-5

a.trait <- a.trait[order(rownames(a.trait)),]
d <- data.frame(name=tree$tip.label,num=seq(1,length(tree$tip.label),by=1))
dd <- d[order(d[,1]) , ]
a.trait$new <- dd[,2]
a.trait <- a.trait[order(a.trait$new),]

dv <- treedata(tree,as.vector(a.trait[,kk]),sort=T)

x <- dv$data
phy <- dv$phy
type = "discrete"
method = "ML"
CI = FALSE
model = if (type == "continuous") "BM" else "ER"
scaled = TRUE
kappa = 1
corStruct = NULL
ip = 0.1 

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

  obj <- list()
    if (kappa != 1) phy$edge.length <- phy$edge.length^kappa

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
        out <- nlminb(rep(ip, length.out = np), function(p) dev(p),
                      lower = rep(0, np), upper = rep(1e50, np))
        obj$loglik <- -out$objective/2
        obj$rates <- out$par
        oldwarn <- options("warn")
        options(warn = -1)
        h <- nlm(function(p) dev(p), p = obj$rates, iterlim = 1,
                 stepmax = 0, hessian = TRUE)$hessian
        options(oldwarn)
        if (any(h == 0))
          warning("The likelihood gradient seems flat in at least one dimension (gradient null):\ncannot compute the standard-errors of the transition rates.\n")
        else obj$se <- sqrt(diag(solve(h)))
        obj$index.matrix <- index.matrix
        if (CI) {
            obj$lik.anc <- dev(obj$rates, TRUE)
            colnames(obj$lik.anc) <- lvls
        }
    }
ace(dv$data,dv$phy,type='discrete',ip=.1,CI=FALSE)

out$objective/-2

data(geospiza)
name <- c("geo")
tree <- geospiza$geospiza.tree
a.trait <- geospiza$geospiza.data

##Which column of data do you want?
kk <- 1

if (dim(a.trait)[2] == 1){a.trait <- data.frame(a.trait,well=rep(1,length(a.trait)))}
nc <- name.check(tree,a.trait)
if (nc[1]=="OK"){nc$Tree.not.data <- NULL}
tree <- drop.tip(tree,nc$Tree.not.data)
tree$edge.length[tree$edge.length<1e-5]=1e-5

a.trait <- a.trait[order(rownames(a.trait)),]
d <- data.frame(name=tree$tip.label,num=seq(1,length(tree$tip.label),by=1))
dd <- d[order(d[,1]) , ]
a.trait$new <- dd[,2]
a.trait <- a.trait[order(a.trait$new),]

dv <- treedata(tree,a.trait[,kk],sort=T)
tree <- dv$phy
trait <- data.frame(taxa=rownames(dv$data),as.vector(dv$data))

source("ouch.wrap.R")

ouch.wrap(tree,trait,model=c("ou1"))

rm(list=ls())

data(geospiza)
name <- c("geo")
tree <- geospiza$geospiza.tree
a.trait <- data.frame(geospiza$geospiza.data,T=as.factor(geospiza$geospiza.data[,1]>4.2))
treeTransform <- "delta"
model <- "ER"

##Which column of data do you want?
kk <- 6

if (dim(a.trait)[2] == 1){a.trait <- data.frame(a.trait,well=rep(1,length(a.trait)))}
nc <- name.check(tree,a.trait)
if (nc[1]=="OK"){nc$Tree.not.data <- NULL}
tree <- drop.tip(tree,nc$Tree.not.data)
tree$edge.length[tree$edge.length<1e-5]=1e-5

a.trait <- a.trait[order(rownames(a.trait)),]
d <- data.frame(name=tree$tip.label,num=seq(1,length(tree$tip.label),by=1))
dd <- d[order(d[,1]) , ]
a.trait$new <- dd[,2]
a.trait <- a.trait[order(a.trait$new),]

dv <- treedata(tree,as.vector(a.trait[,kk]),sort=T)

fitDiscrete(dv$phy,dv$data,model="ER",treeTransform="delta")

fitContinuous(dv$phy,dv$data,model="delta")

opt <- optim(
                   par=c(sqrt.alpha,sigma),
                   fn = function (par) {
                     ou.lik.fn(
                               tree=tree,
                               alpha=sym.par(par[seq(nalpha)]),
                               sigma=sym.par(par[nalpha+seq(nsigma)]),
                               beta=beta,
                               dat=dat
                               )$deviance
                   },
                   gr=NULL,
                   hessian=hessian,
                   method=method[4]
                   )
