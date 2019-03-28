\name{df.clean}
\alias{df.clean}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Remove Columns in a Data Frame Based on a Threshold
%%  ~~function to do ... ~~
}
\description{This function removes columns in a data frame if the values they contain are on average smaller/greater than a threshold specified by the user.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
df.clean(df=NULL, threshold=NULL, smaller=TRUE, index=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{a data frame or matrix upon which the removal process is applied.
%%     ~~Describe \code{x} here~~
}
  \item{threshold}{numeric; for each column in the data frame, the mean is compared against the threshold.}
  \item{smaller}{logical; default is TRUE. If set FALSE, columns with a mean greater than the threshold will be removed.}
  \item{index}{integers connected by a colon indicating which columns in the data frame should be considered. If ommitted, all columns are considered by default}
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
x<-rnorm(n=10, mean=0, sd=1)
y<-rnorm(n=10, mean=5, sd=1)
z<-rnorm(n=10, mean=10, sd=1)
data<-cbind(x, y, z)
df.clean(data, threshold=5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
