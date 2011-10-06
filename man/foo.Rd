\name{foo}
\alias{foo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
foo(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x) 
{
    if (dplot) {
        foostore$x1 <<- c(foostore$x1, x[1])
        foostore$x2 <<- c(foostore$x2, x[2])
    }
    t <- deltaTree(tree, delta = exp(x[2]))
    vcv <- vcv.phylo(t)
    vv <- exp(x[1]) * vcv
    diag(vv) <- diag(vv) + meserr^2
    determinantVCV <- det(vv)
    if (determinantVCV == 0) {
        warning("Possibly singular variance-covariance matrix, so giving this particular parameter combination a very bad likelihood score (rather than crashing)")
        return(badLnL)
    }
    mu <- phylogMean(vv, chdata)
    mu <- rep(mu, n)
    -dmvnorm(chdata, mu, vv, log = T)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
