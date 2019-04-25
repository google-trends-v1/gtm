\name{clust.factor}
\alias{clust.factor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Dynamic Factors Extraction based on Time Series Clustering
%%  ~~function to do ... ~~
}
\description{This function extracts dynamic factors based on time series clustering. Time series are first clustered based on partitional, hierarchical, tadpole or fuzzy methods, and one dynamic factor is then extracted within each cluster, based on either a two-step method as proposed by Doz, Gianone & Reichlin (2011) or a simple aggregation method.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
clust.factor(data=NULL, fac.num=3, method="two-step", clustor.type="partitional", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{a multivariate time series or a numeric matrix to which the clustering and extraction are applied.}
  \item{fac.num}{number of desired clusters, which is also the number of desired factors.}
  \item{method}{method to be applied for extraction of dynamic factors. Options include \code{"two-step"} and \code{"aggregation"}. Defualt is \code{"two-step"}.}
  \item{clustyor.type}{the type of clustering method to be applied. Options include \code{"partitional"}, \code{"hierarchical"}, \code{"tadpole"}, and \code{"fuzzy"}. Default is \code{"partitional"}. }
  \item{\dots}{extra arguments which are passed to the \code{\link{tsclust}} function.}
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
\author{Zehua Wu
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
