\name{df.merge}
\alias{df.merge}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Recursively Merging Data from a List
%%  ~~function to do ... ~~
}
\description{This function recursively merges data frames in a list into a multivariate dataset. The function automatically detects error in the merging process and reports where the error occurs.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
df.merge(df.list=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df.list}{a list that contains data frames to be merged. Those data frames should share one and only one identical column.
%%     ~~Describe \code{x} here~~
}
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
dir="~/Documents/Mydata"
dfs<-load.data(path=dir, pattern=NULL, merge=FALSE)
data.merged<-df.merge(dfs)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
