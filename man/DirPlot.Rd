\name{DirPlot}
\alias{DirPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Visualization of the Directional Forecasts and the Original Time Series
%%  ~~function to do ... ~~
}
\description{This function plots the original time series, and draws dots that represent the directional forecasts. The color and shape of the dots represent direction up or direction down, and the size of the dots represents whether or not the directional forecast is correct.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
DirPlot(forecasts, original, size.range, title,
        legend.position="bottom")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{forecasts}{an object of class "Maeforecast". Can be returned by functions \code{\link{maeforecast}} and \code{\link{Bagging}}.}
  \item{original}{a numeric vector or a \code{ts} or \code{zoo} object that contains the original time series of the dependent variable. Note that this should not be the log-differenced series used in making the forecasts.}
  \item{size.range}{a numeric vector of length 2 that contains desirable sizes for the incorrect and then the correct directional forecast dots. Note that the order matters here.}
  \item{title}{the title of the plot, which can be omitted.}
  \item{legend.position}{the position of legends (\code{"none"}, \code{"left"}, \code{"right"}, \code{"bottom"}, \code{"top"}, or two-element numeric vector)}
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
