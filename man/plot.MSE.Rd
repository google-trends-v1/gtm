\name{plot.RecMSE}
\alias{plot.RecMSE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot the Ratio or Difference between Recursive MSEs of two Models
%%  ~~function to do ... ~~
}
\description{This function plots the ratio of the difference between the recursive MSEs of two models.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
## S3 method for class 'RecMSE'
plot(.data, type="ratio", title=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{.data}{an object of class "RecMSE" returned by the function \code{\link{RecursiveMSE}}.}
  \item{type}{a character indicating whether the MSE ratio or the MSE difference should be plotted.}
  \item{title}{the title of the plot. If omitted, the title will be "Recursive MSE"}
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
