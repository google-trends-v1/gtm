\name{maeforecast.postlasso}
\alias{maeforecast.postlasso}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Economic Forecasting with High-Dimensional Data (Post-Lasso)
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
maeforecast.postlasso(data, w_size, window="recursive", y.index=1,
            t.select, h=0, standardize=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{a data frame or a matrix; the first column should contain the time series variable for which the forecasts are to be made. Other columns should contain the covariates.}
  \item{w_size}{numeric, indicating the index where the forecasting should begin. If the first point forecast should be made at the 73th observation, for example, \code{w_size} should be set to be 72.}
  \item{window}{character, indicating the forecasting scheme to be applied. Options include \code{"recursive"}, \code{"rolling"}, and \code{"fixed"}.}
  \item{y.index}{numeric, indicating the column position of the time series for which the forecasts are made (Y). Defualt is \code{1}.}
  \item{t.select}{number of covariates to be included. If omitted, every covariate will be included. Otherwise, a regression between the dependant variable, its lag and each covariate will be run and a statistical test will be applied for the significance of the covariate's coefficient. The covariates will then be ranked based on their test statistics, and \code{t.select} number of them with the highest test statistics will be included in the machine learning algorithms. Note that the forecasting horizon \code{h} is considered. This selection process is only applied at the beggining of the algorithms; that is, the list of covariates will not be updated with the rolling or recursive windows.}
  \item{h}{forecasting horizon. Default is \code{0}. }
  \item{standardize}{logical, indicating whether the data matrix should be scaled before the model is fitting, for the use of variable selection/shrinkage models. Default is \code{TRUE}.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{Forecasts }{data matrix, containing the point forecasts, realized values, forecast errors, signs of the forecasts and realized values, and success in predicting the signs.}
  \item{MSE }{numeric, mean squred error of the point forecasts.}
  \item{SRatio}{numeric, success ratio of the point forecasts. Success is claimed when the point forecasts and realized values have the same sign.}
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
\code{\link{maeforecast}}, \code{\link{maeforecast.lasso}}, \code{\link{maeforecast.alasso}}, \code{\link{maeforecast.postalasso}}, \code{\link{maeforecast.ridge}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
