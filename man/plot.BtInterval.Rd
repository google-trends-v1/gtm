\name{plot.BtInterval}
\alias{plot.BtInterval}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Visualization of Bootstrapped Predictive Intervals
%%  ~~function to do ... ~~
}
\description{This function plots the bootstrapped predictive intervals constructed by the \code{\link{bt.interval}} function.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
## S3 method for class 'BtInterval'
plot(.data, out.true=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{.data}{an object of class "BtInterval" which should be returned by the \code{\link{bt.interval}} function.
%%     ~~Describe \code{.data} here~~
}
  \item{out.true}{logical, whether the realized values should be plotted. If \code{TRUE}, the realized values will be plotted as an additional blue dashed line.
%%     ~~Describe \code{out.true} here~~
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

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
