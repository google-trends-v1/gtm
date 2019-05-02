\name{Metrics}
\alias{Metrics}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{MSE and Success Ratio of Forecasts
%%  ~~function to do ... ~~
}
\description{
This function calculates mean squared error (MSE) and the success ratio of a given set of forecasts.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Metrics(pred=NULL, true=NULL, h=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pred}{vector of forecasts.
%%     ~~Describe \code{x} here~~
}
  \item{true}{vector of true values.}
  \item{h}{forecasting horizon used in making the forecasts. Default is \code{0}. }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
This function returns the same outputs as does the \code{\link{maeforecast}} function. More specifically, it returns a list containing the following objects:
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

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
