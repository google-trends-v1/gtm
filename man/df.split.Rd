\name{df.split}
\alias{df.split}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Splitting a Multivariate Time Series Data Frame into Individual Time Series
%%  ~~function to do ... ~~
}
\description{This function splits a multivariate time series data frame into individual data frames containing each time series and a common variable shared by them. The individual data frames are stored and returned in a list. A common example of the "common variable" is Date.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
df.split(dataframe=NULL, by=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataframe}{the multivariate time series dataframe to be splitted.
%%     ~~Describe \code{x} here~~
}
  \item{by}{character; the common variable to be shared by the individual data frames.}
}
\details{Sometimes, the merge function in R returns zero observations when dealing with data frames containing large amounts of columns. A good way to go about this is to split those large data frames into individual data frames and merge them with the reduce function.
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
Zehua Wu
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{df.merge}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data("US_Unemployment")
data.list<-df.split(dataframe=US_Unemployment, by="Date")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
