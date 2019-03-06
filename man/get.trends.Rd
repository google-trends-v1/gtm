\name{get.trends}
\alias{get.trends}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Getting Google Trends
%%  ~~function to do ... ~~
}
\description{The function downloads time series from Google Trends based on a list of queries provided by the user, and stores them as .csv files either under the current working directory or elsewhere as speficied otherwise.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
get.trends(queries=NULL, geo="US", time="all", path=getwd())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{quries}{a vector of quries to get Google Trends for.
%%     ~~Describe \code{x} here~~
}
  \item{geo}{character, indicating the geographic location for the queries.
}
  \item{time}{character, indicating the time span for the queries.
}
  \item{path}{string, denoting the directory under which the csv files should be stored.}
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
Zehua Wu
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
list<-c("unemployment", "us unemployment")
dir="~/Documents"
get.trends(queries=list, geo="US", time="all", path=dir)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
