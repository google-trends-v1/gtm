\name{maeforecast.lasso}
\alias{maeforecast.lasso}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Economic Forecasting with High-Dimensional Data (Lasso)
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
maeforecast.lasso(data, w_size, window="recursive", y.index=1,
            standardize=TRUE, lambda)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{a data frame or a matrix; the first column should contain the time series variable for which the forecasts are to be made. Other columns should contain the covariates.}
  \item{w_size}{numeric, indicating the index where the forecasting should begin. If the first point forecast should be made at the 73th observation, for example, \code{w_size} should be set to be 72.}
  \item{window}{character, indicating the forecasting scheme to be applied. Options include \code{"recursive"}, \code{"rolling"}, and \code{"fixed"}.}
  \item{y.index}{numeric, indicating the column position of the time series for which the forecasts are made (Y). Defualt is \code{1}.}
  \item{standardize}{logical, indicating whether the data matrix should be scaled before the model is fitting, for the use of variable selection/shrinkage models. Default is \code{TRUE}.}
  \item{lambda}{optional user-supplied lambda sequence; default is NULL, and the function \code{glmnet} chooses its own sequence.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{Forecasts }{data matrix, containing the point forecasts, realized values, forecast errors, signs of the forecasts and realized values, and success in predicting the signs.}
  \item{MSE }{numeric, mean squred error of the point forecasts.}
  \item{SRatio}{numeric, success ratio of the point forecasts. Success is claimed when the point forecasts and realized values have the same sign.}
  \item{Variables}{list, containing the predictors selected by the shrinkage model at every iteration.}
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
\code{\link{maeforecast}}, \code{\link{maeforecast.alasso}}, \code{\link{maeforecast.postlasso}}, \code{\link{maeforecast.postalasso}}, \code{\link{maeforecast.ridge}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
