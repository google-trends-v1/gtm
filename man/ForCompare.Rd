\name{ForCompare}
\alias{ForCompare}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generate Model Comparison Table
%%  ~~function to do ... ~~
}
\description{This model generates a comparison table that reports the MSE, Success Ratio and MSE Ratio for a group of models being compared.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ForCompare(..., benchmark=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{...}{one or more output from functions \code{\link{maeforecast}}, \code{\link{Metrics}}, and \code{\link{Bagging}} to generate the comparison table for. If omitted, the \code{MSERatio} will not be computed.}
  \item{benchmark}{either an output from the above listed functions, or a numeric indicating the MSE for the benchmark model, or an integer indicating which model listed in \code{\dots} should be treated as the benchmark.}
%%     ~~Describe \code{x} here~~
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{

This function returns a data table containing three columns.

\item{MSE}{the mean squared error of point forecasts for each model being compared.}
\item{SRatio}{the success ratio of the directional forecasts for each model being compared.}
\item{MSERatio}{the ratio of each model's MSE against that of a benchmark.}
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
AR.For<-maeforecast(mydata, w_size=72, window="recursive",
        model="ar")
Lasso.For<-maeforecast(mydata, w_size=72, window="recursive",
        model="lasso")
ForCompare(AR.For, Lasso.For)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
