\name{df.clean}
\alias{df.clean}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Remove Columns in a Data Frame Based on a Threshold
%%  ~~function to do ... ~~
}
\description{This function removes variables in a data frame if the values they contain have averages that are smaller/greater than a threshold, or if they are highly correlated with other variables.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
df.clean(df=NULL, type="correlation", threshold=NULL,
        smaller=TRUE, index=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{a data frame or matrix upon which the removal process is applied.}
  \item{type}{a character indicating what type of removal criterion should be applied. Options include \code{"correlation"} and \code{"mean"}.}
  \item{threshold}{numeric; if \code{type="mean"}, for each column in the data frame, the mean is compared against the threshold; if \code{type="correlation"}, the pair-wise correlation coefficients among the variables are correlated.}
  \item{smaller}{logical; for the use of \code{type="mean"} only. Default is \code{TRUE}. If set \code{FALSE}, columns with a mean greater than the threshold will be removed.}
  \item{index}{integers connected by a colon indicating which columns in the data frame should be considered. If ommitted, all columns are considered by default.}
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
\code{\link{df.rm}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(US_Unemployment)
mydata<-df.clean(US_Unemployment, threshold=0.95,
        index=3:ncol(US_Unemployment))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
