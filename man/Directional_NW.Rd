\name{Directional_NW}
\alias{Directional_NW}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Weighted and Unweighted Directional Forecast Test with Newey-West HAC Estimator
%%  ~~function to do ... ~~
}
\description{This function implements formal statistical tests to check whether the directional forecasts are able to improve upon driftless random walk forecasts. Two test options include weighted and unweighted directional forecasts, and the test statistic is calculated with the Newey-West HAC estimator to control for serial correlation.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Directional_DW(forecasts, p=1, weighted=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{forecasts}{an object of class "Maeforecast" or "MaeBagging". Can be returned by functions \code{\link{maeforecast}}, \code{\link{Bagging}}, and \code{\link{Metrics}}.}
  \item{p}{the largest lag used in the forecasting model. Default is \code{1}.}
  \item{weighted}{logical. If \code{TRUE}, weighted directional forecast test will be applied; otherwise the test will be applied upon unweighted directional forecasts. Default is \code{TRUE}.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{test.statistic}{the test statistic.}
\item{pvalue}{the p value.}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Kim, Young Ju, Zhipeng Liao, and Aaron Tornell. 2014. “Speculators Positions and Exchange Rate
Forecasts: Beating Random Walk Models.” Unpublished Manuscript UCLA.
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
\examples{
AR.For<-maeforecast(mydata, w_size=72, window="recursive",
        model="ar")
Directional_NW(AR.FOR, p=1, weighted=T)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
