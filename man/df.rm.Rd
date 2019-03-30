\name{df.rm}
\alias{df.rm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Remove Time Series Components for Variables in a Data Frame
%%  ~~function to do ... ~~
}
\description{This function removes given time series component(s) from selected columns in a data frame or data matrix. Each column in the data frame will be automatically converted to a time series object, and a decomposition process will then be applied to it. Users can choose to remove the seasonal component, the trend, or both components from the time series. A new data frame will be returned.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
df.rm(df, index, rm="both", frequency=12)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{a data frame or matrix upon which the removal process is applied.}
  \item{index}{integers connected by a colon indicating which columns in the data frame should be considered. If ommitted, all columns are considered by default.}
  \item{rm}{a character indicating which time series component(s) should be removed. Options include \code{"seasonal"}, \code{"trend"}, and \code{"both"}. Defualt is \code{"both"}.}
  \item{frequency}{an integer indicating the frequency of the time series. Default is \code{12} for montly observations.}
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
\code{\link{df.clean}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data("US_Unemployment")
mydata<-df.rm(df=US_Unemployment, index=2:ncol(US_Unemployment))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
