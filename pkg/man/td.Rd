\name{td}
\alias{td}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(td)}
\format{
  The format is:
List of 2
 $ phy :List of 4
  ..$ edge       : int [1:24, 1:2] 14 15 16 17 18 19 20 21 22 22 ...
  ..$ edge.length: num [1:24] 0.0492 0.0686 0.134 0.1035 0.0355 ...
  ..$ Nnode      : int 12
  ..$ tip.label  : chr [1:13] "fuliginosa" "fortis" "magnirostris" "conirostris" ...
  ..- attr(*, "class")= chr "phylo"
  ..- attr(*, "order")= chr "cladewise"
 $ data: num [1:13, 1] 4.4 4.35 4.22 4.26 4.24 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:13] "fuliginosa" "fortis" "magnirostris" "conirostris" ...
  .. ..$ : NULL
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
data(td)
## maybe str(td) ; plot(td) ...
}
\keyword{datasets}
