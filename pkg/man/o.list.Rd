\name{o.list}
\alias{o.list}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(o.list)}
\format{
  The format is:
List of 1
 $ :List of 1
  ..$ :'data.frame':	6 obs. of  10 variables:
  .. ..$ par    :List of 6
  .. .. ..$ : num [1:2] -3.88 3.05
  .. .. ..$ : num [1:2] -3.88 3.05
  .. .. ..$ : num [1:2] -3.88 3.05
  .. .. ..$ : num [1:2] -3.88 3.05
  .. .. ..$ : num [1:2] -3.88 3.05
  .. .. ..$ : num [1:2] -3.88 3.05
  .. ..$ fvalues:List of 6
  .. .. ..$ : num -10.4
  .. .. ..$ : num -10.4
  .. .. ..$ : num -10.4
  .. .. ..$ : num -10.4
  .. .. ..$ : num -10.4
  .. .. ..$ : num -10.4
  .. ..$ method :List of 6
  .. .. ..$ : chr "spg"
  .. .. ..$ : chr "L-BFGS-B"
  .. .. ..$ : chr "Rvmmin"
  .. .. ..$ : chr "Rcgmin"
  .. .. ..$ : chr "bobyqa"
  .. .. ..$ : chr "nlminb"
  .. ..$ fns    :List of 6
  .. .. ..$ : num 18
  .. .. ..$ : Named int 13
  .. .. .. ..- attr(*, "names")= chr "function"
  .. .. ..$ : num 30
  .. .. ..$ : num 43
  .. .. ..$ : int 62
  .. .. ..$ : Named int 13
  .. .. .. ..- attr(*, "names")= chr "function"
  .. ..$ grs    :List of 6
  .. .. ..$ : logi NA
  .. .. ..$ : Named int 13
  .. .. .. ..- attr(*, "names")= chr "gradient"
  .. .. ..$ : num 15
  .. .. ..$ : num 24
  .. .. ..$ : logi NA
  .. .. ..$ : Named int 31
  .. .. .. ..- attr(*, "names")= chr "gradient"
  .. ..$ itns   :List of 6
  .. .. ..$ : num 15
  .. .. ..$ : NULL
  .. .. ..$ : NULL
  .. .. ..$ : NULL
  .. .. ..$ : NULL
  .. .. ..$ : int 11
  .. ..$ conv   :List of 6
  .. .. ..$ : num 0
  .. .. ..$ : int 0
  .. .. ..$ : num 0
  .. .. ..$ : num 0
  .. .. ..$ : num 0
  .. .. ..$ : int 0
  .. ..$ KKT1   :List of 6
  .. .. ..$ : logi TRUE
  .. .. ..$ : logi TRUE
  .. .. ..$ : logi TRUE
  .. .. ..$ : logi TRUE
  .. .. ..$ : logi TRUE
  .. .. ..$ : logi TRUE
  .. ..$ KKT2   :List of 6
  .. .. ..$ : logi TRUE
  .. .. ..$ : logi TRUE
  .. .. ..$ : logi TRUE
  .. .. ..$ : logi TRUE
  .. .. ..$ : logi TRUE
  .. .. ..$ : logi TRUE
  .. ..$ xtimes :List of 6
  .. .. ..$ : Named num 0.371
  .. .. .. ..- attr(*, "names")= chr "user.self"
  .. .. ..$ : Named num 0.486
  .. .. .. ..- attr(*, "names")= chr "user.self"
  .. .. ..$ : Named num 0.582
  .. .. .. ..- attr(*, "names")= chr "user.self"
  .. .. ..$ : Named num 0.905
  .. .. .. ..- attr(*, "names")= chr "user.self"
  .. .. ..$ : Named num 0.467
  .. .. .. ..- attr(*, "names")= chr "user.self"
  .. .. ..$ : Named num 0.357
  .. .. .. ..- attr(*, "names")= chr "user.self"
  .. ..- attr(*, "details")=List of 6
  .. .. ..$ :List of 14
  .. .. .. ..$ par        : num [1:2] -3.88 3.05
  .. .. .. ..$ value      : num -10.4
  .. .. .. ..$ convergence: int 0
  .. .. .. ..$ message    : chr "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"
  .. .. .. ..$ conv       : int 0
  .. .. .. ..$ fevals     : Named int 13
  .. .. .. .. ..- attr(*, "names")= chr "function"
  .. .. .. ..$ gevals     : Named int 13
  .. .. .. .. ..- attr(*, "names")= chr "gradient"
  .. .. .. ..$ kkt1       : logi TRUE
  .. .. .. ..$ kkt2       : logi TRUE
  .. .. .. ..$ ngatend    : num [1:2] 6.15e-06 -1.05e-06
  .. .. .. ..$ nhatend    : num [1:2, 1:2] 6.5 0.39 0.39 1.37
  .. .. .. ..$ evnhatend  : num [1:2] 6.53 1.34
  .. .. .. ..$ systime    : Named num 0.486
  .. .. .. .. ..- attr(*, "names")= chr "user.self"
  .. .. .. ..$ method     : chr "L-BFGS-B"
  .. .. ..$ :List of 15
  .. .. .. ..$ par        : num [1:2] -3.88 3.05
  .. .. .. ..$ convergence: int 0
  .. .. .. ..$ message    : chr "relative convergence (4)"
  .. .. .. ..$ conv       : int 0
  .. .. .. ..$ value      : num -10.4
  .. .. .. ..$ fevals     : Named int 13
  .. .. .. .. ..- attr(*, "names")= chr "function"
  .. .. .. ..$ gevals     : Named int 31
  .. .. .. .. ..- attr(*, "names")= chr "gradient"
  .. .. .. ..$ niter      : int 11
  .. .. .. ..$ kkt1       : logi TRUE
  .. .. .. ..$ kkt2       : logi TRUE
  .. .. .. ..$ ngatend    : num [1:2] -9.53e-09 1.50e-08
  .. .. .. ..$ nhatend    : num [1:2, 1:2] 6.5 0.39 0.39 1.37
  .. .. .. ..$ evnhatend  : num [1:2] 6.53 1.34
  .. .. .. ..$ systime    : Named num 0.357
  .. .. .. .. ..- attr(*, "names")= chr "user.self"
  .. .. .. ..$ method     : chr "nlminb"
  .. .. ..$ :List of 17
  .. .. .. ..$ par         : num [1:2] -3.88 3.05
  .. .. .. ..$ value       : num -10.4
  .. .. .. ..$ gradient    : num 1.12e-05
  .. .. .. ..$ fn.reduction: num 12.7
  .. .. .. ..$ convergence : num 0
  .. .. .. ..$ message     : chr "Successful convergence"
  .. .. .. ..$ conv        : num 0
  .. .. .. ..$ fevals      : num 18
  .. .. .. ..$ gevals      : logi NA
  .. .. .. ..$ niter       : num 15
  .. .. .. ..$ kkt1        : logi TRUE
  .. .. .. ..$ kkt2        : logi TRUE
  .. .. .. ..$ ngatend     : num [1:2] 5.16e-07 -1.13e-05
  .. .. .. ..$ nhatend     : num [1:2, 1:2] 6.5 0.39 0.39 1.37
  .. .. .. ..$ evnhatend   : num [1:2] 6.53 1.34
  .. .. .. ..$ systime     : Named num 0.371
  .. .. .. .. ..- attr(*, "names")= chr "user.self"
  .. .. .. ..$ method      : chr "spg"
  .. .. ..$ :List of 15
  .. .. .. ..$ par        : num [1:2] -3.88 3.05
  .. .. .. ..$ value      : num -10.4
  .. .. .. ..$ counts     : num [1:2] 43 24
  .. .. .. ..$ convergence: num 0
  .. .. .. ..$ message    : chr "Rcgmin seems to have converged"
  .. .. .. ..$ conv       : num 0
  .. .. .. ..$ fevals     : num 43
  .. .. .. ..$ gevals     : num 24
  .. .. .. ..$ kkt1       : logi TRUE
  .. .. .. ..$ kkt2       : logi TRUE
  .. .. .. ..$ ngatend    : num [1:2] -1.78e-07 1.32e-07
  .. .. .. ..$ nhatend    : num [1:2, 1:2] 6.5 0.39 0.39 1.37
  .. .. .. ..$ evnhatend  : num [1:2] 6.53 1.34
  .. .. .. ..$ systime    : Named num 0.905
  .. .. .. .. ..- attr(*, "names")= chr "user.self"
  .. .. .. ..$ method     : chr "Rcgmin"
  .. .. ..$ :List of 16
  .. .. .. ..$ par        : num [1:2] -3.88 3.05
  .. .. .. ..$ value      : num -10.4
  .. .. .. ..$ counts     : num [1:2] 30 15
  .. .. .. ..$ convergence: num 0
  .. .. .. ..$ message    : chr "Converged"
  .. .. .. ..$ bdmsk      : num [1:2] 1 1
  .. .. .. ..$ conv       : num 0
  .. .. .. ..$ fevals     : num 30
  .. .. .. ..$ gevals     : num 15
  .. .. .. ..$ kkt1       : logi TRUE
  .. .. .. ..$ kkt2       : logi TRUE
  .. .. .. ..$ ngatend    : num [1:2] -1.28e-06 -2.08e-07
  .. .. .. ..$ nhatend    : num [1:2, 1:2] 6.5 0.39 0.39 1.37
  .. .. .. ..$ evnhatend  : num [1:2] 6.53 1.34
  .. .. .. ..$ systime    : Named num 0.582
  .. .. .. .. ..- attr(*, "names")= chr "user.self"
  .. .. .. ..$ method     : chr "Rvmmin"
  .. .. ..$ :List of 16
  .. .. .. ..$ par      : num [1:2] -3.88 3.05
  .. .. .. ..$ fval     : num -10.4
  .. .. .. ..$ feval    : int 62
  .. .. .. ..$ ierr     : int 0
  .. .. .. ..$ msg      : chr "Normal exit from bobyqa"
  .. .. .. ..$ conv     : num 0
  .. .. .. ..$ fevals   : int 62
  .. .. .. ..$ gevals   : logi NA
  .. .. .. ..$ value    : num -10.4
  .. .. .. ..$ kkt1     : logi TRUE
  .. .. .. ..$ kkt2     : logi TRUE
  .. .. .. ..$ ngatend  : num [1:2] -2.12e-07 -5.22e-08
  .. .. .. ..$ nhatend  : num [1:2, 1:2] 6.5 0.39 0.39 1.37
  .. .. .. ..$ evnhatend: num [1:2] 6.53 1.34
  .. .. .. ..$ systime  : Named num 0.467
  .. .. .. .. ..- attr(*, "names")= chr "user.self"
  .. .. .. ..$ method   : chr "bobyqa"
  .. .. .. ..- attr(*, "class")= chr [1:2] "bobyqa" "minqa"
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(o.list)
## maybe str(o.list) ; plot(o.list) ...
}
\keyword{datasets}
