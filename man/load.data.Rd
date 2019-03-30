\name{load.data}
\alias{load.data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Loading Data from a Directory
%%  ~~function to do ... ~~
}
\description{The function automatically loads .csv files, either under the current working directory or the directory specified by the user, as data frames. If a pattern is specified, only files with names matching the pattern will be loaded. The argument merge is used to indicate whether the data frames should be stored in a single list or merged together in one data frame.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
load.data(pattern=NULL, path=getwd(), merge=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{path}{string, indicating the directory under which the files to be loaded are stored. If omitted, files will be loaded from the current working directory.}
    \item{pattern}{character, indicating the pattern to be mathced in the file names. If specified, only files with names matching the pattern will be loaded.
%%     ~~Describe \code{x} here~~
}
  \item{merge}{logical, \code{FALSE} by default. If \code{TRUE}, data frames will be merged based on a common variable name.
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
\author{
Zehua Wu
%%  ~~who you are~~
}
\note{
R has trouble processing variable names that contain certain special characters, and this might interfere with the merging function. If this happens, the function will stop running and reporting where the error takes place.
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{df.merge}}, \code{\link{df.split}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
dir="~/Documents/Mydata"
data<-load.data(path=dir, pattern=NULL, merge=TRUE)
write.csv(data, "Merged.csv")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
