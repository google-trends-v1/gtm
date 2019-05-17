\name{RecursiveMSE}
\alias{RecursiveMSE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Recursive MSE
%%  ~~function to do ... ~~
}
\description{This function calculates the recursive MSEs of the forecasts made by two models, one of which is the benchmark. Recursive MSE is calculated as the mean squared errors of the point forecasts by adding one point forecast at a time, starting with the first out-of-sample forecast. One can use this to see how the MSEs of the two models evolve relative to each other as more out-of-sample forecasts are made.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
RecursiveMSE(forecast.main, forecast.benchmark)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{forecast.main}{an onject of class "Maeforecast". The forecasts that are to be compared against the benchmark.}
  \item{forecast.benchmark}{an object of class "Maeforecast". The forecasts made by the benchmark model. Note that the two sets of forecasts should be made with the same window, window size and forecasting horizon.
%%     ~~Describe \code{forecast.benchmark} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{This function returns an object of class "RecMSE" that contains the recursive MSEs of the two sets of forecasts, as well as their differences and ratios.
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
\code{\link{plot.RecMSE}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
