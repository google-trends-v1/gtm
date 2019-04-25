\name{ForCompare}
\alias{ForCompare}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generate Model Comparison Table
%%  ~~function to do ... ~~
}
\description{This model generates a comparison table that reports the MSE, Success Ratio, MSE Ratio, and some test statistics for a group of models being compared.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ForCompare(..., benchmark.index=NULL, test=c("unweighted", "binary"),
          h=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{...}{one or more output from functions \code{\link{maeforecast}}, \code{\link{Metrics}}, and \code{\link{Bagging}} to generate the comparison table for.}
  \item{benchmark.index}{an integer indicating which model listed in \code{\dots} should be treated as the benchmark. If omitted, the \code{MSERatio} and \code{DMW} will not be computed.}
  \item{test}{statistical test p values to be reported. Options include \code{"weighted"} (weighted directional forecast test) and \code{"binary"} (unweighted directional forecast test). See \code{\link{Directional_NW}} for details.}
  \item{h}{a numeric indicating the forecast horizon used in the models.}
%%     ~~Describe \code{x} here~~
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{

This function returns a data frame with the following potential columns

\item{MSE}{the mean squared error of point forecasts for each model being compared.}
\item{SRatio}{the success ratio of the directional forecasts for each model being compared.}
\item{MSERatio}{the ratio of each model's MSE against that of a benchmark.}
\item{DMW}{the p values returned from DMW tests against a benchmark indicated by \code{benchmark.index}. The null hypothesis is that the model being compared has the same forecast accuracy as the benchmark; the alternative hupothesis is that the model being compared is better than the benchmark.}
\item{Weighted}{p value from weighted directional forecast test.}
\item{Unweighted}{p value from unweighted directional forecast test.}
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
