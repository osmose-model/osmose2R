\name{osmose2R}
\alias{osmose2R}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Get OSMOSE outputs into R}
\description{
This function parse all output directory structure (as in OSMOSE 3.0) and read all the data into R.}
\usage{
osmose2R(path = NULL)
}
\arguments{
  \item{path}{The path to the OSMOSE's output directory.}
}
\details{
In interactive mode, without arguments, the function will prompt a window to choose the output directory.
}
\value{
An object of class \emph{osmose} is created, containing all the outputs and model information. Information can be extracted with the proper methods for the class osmose. 
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Ricardo Oliveros-Ramos
}
\note{
This function is in develop stage.
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
}
\keyword{ osmose }
