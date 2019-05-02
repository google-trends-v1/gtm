\name{plot.Maeforecast}
\alias{plot.Maeforecast}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot Results from Maeforecast Functions
%%  ~~function to do ... ~~
}
\description{This function plots the point forecasts along side the realized values using \code{ggplot}.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
## S3 method for class 'Maeforecast'
plot(forecasts, start=NULL, frequency='month',
    forecast.lab="Forecasts", true.lab="Realized",
    x.lab="Time", y.lab="Value")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{forecasts}{an object of class "Maeforecast". Can be returned by functions \code{\link{maeforecast}} and \code{\link{Bagging}}.}
  \item{start}{the starting time of the forecasts. Should have the format \code{"\%Y-\%m-\%d"}. If omitted, time will be shown as integers starting from 1.}
  \item{frequency}{if \code{start} is not omitted, \code{frequency} is used to contruct the time sequence, which indicates the frequency of the time series. Options include \code{"year"}, \code{"month"}, \code{"week"}, and \code{"day"}.}
  \item{forecast.lab}{character, label of the forecasts. Default is \code{"Forecasts"}.}
  \item{true.lab}{character, label of the realized values. Default is \code{"Realized"}.}
  \item{x.lab}{character, label of the x-axis. Default is \code{"Time"}.}
  \item{y.labe}{character, label of the y-axis. Default is \code{"Value"}.}
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
AR.For<-maeforecast(mydata, w_size=72, window="recursive",
        model="ar")
plot(AR.For, start="2010-02-01", frequency="month")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
