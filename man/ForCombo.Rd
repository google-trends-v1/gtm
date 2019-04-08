\name{ForCombo}
\alias{ForCombo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Forecast Combination based on Mean or Median
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ForCombo(..., type="mean")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{...}{one or more output from functions \code{\link{maeforecast}}, \code{\link{Metrics}}, and \code{\link{Bagging}} for which the combined forecasts are calculated.
%%     ~~Describe \code{x} here~~
}
  \item{type}{a character indicating whether the mean or median should be used for the combination. Default is \code{"mean"}.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{Forecasts }{data matrix, containing the point forecasts, realized values, forecast errors, signs of the forecasts and realized values, and success in predicting the signs.}
  \item{MSE }{numeric, mean squred error of the point forecasts.}
  \item{SRatio}{numeric, success ratio of the point forecasts. Success is claimed when the point forecasts and realized values have the same sign.}
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
ForCombo(AR.For, Lasso.For, type="mean")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
