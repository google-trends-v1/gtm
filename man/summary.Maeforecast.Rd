\name{summary.Maeforecast}
\alias{summary.Maeforecast}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summarize Results from Maeforecast Functions
%%  ~~function to do ... ~~
}
\description{This function summarizes the results from Maeforecast Functions. It prints out the MSE and the success ratio, as well as some information of the model used in making the forecasts.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
## S3 method for class 'Maeforecast'
summary(x, digits=7)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object of class "Maeforecast". Should be returned by functions \code{\link{maeforecast}} and \code{\link{Bagging}}.}
  \item{digits}{integer indicating the number of decimal places. Defualt is \code{7}.}
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
summary(AR.For)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
