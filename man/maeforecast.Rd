\name{maeforecast}
\alias{maeforecast}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Economic Forecasting with High-Dimensional Data
%%  ~~function to do ... ~~
}
\description{This function calculates one-step-ahead point forecasts with recursive, rolling or fixed windows in a data-rich environment. Supported machine learning algorithms include Lasso, Ridge, Adaptive Lasso, Adaptive Elastic Net, as well as their interactions with simple AR(1) models. Dynamic factor models, Regression Trees and Random Forest algorithms are also supported. Out-of-sample forecasts are returned alongside mean squred errors and success ratios.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
maeforecast(data, model="ar", w_size, window="recursive", ...)

maeforecast.ar(data, w_size, window="recursive")
maeforecast.lasso(data, w_size, window="recursive")
maeforecast.arlasso(data, w_size, window="recursive")
maeforecast.ridge(data, w_size, window="recursive")
maeforecast.alasso(data, w_size, window="recursive")
maeforecast.aralasso(data, w_size, window="recursive")
maeforecast.arnet(data, w_size, window="recursive", pred)
maeforecast.dfm(data, w_size, window="recursive", factor.num)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{a data frame or a matrix; the first column should contain the time series variable for which the forecasts are to be made. Other columns should contain the covariates.}
  \item{model}{character, indicating which time series model or machine learning algorithm should be implemented. Options include \code{"ar"}, \code{"lasso"}, \code{"arlasso"} (AR+Lasso), \code{"ridge"}, \code{"alasso"} (Adaptive Lasso), \code{"aralasso"} (AR+Ada-Lasso), \code{"arnet"} (AR+Ada-ElasticNet), \code{"dfm"}(Dynamic Factor Model). Default is \code{"ar"}. Note that \code{maeforecast.ar()} is equivalent to \code{maeforecast(model="ar")}, and so on.}
  \item{w_size}{numeric, indicating the index where the forecasting should begin. If the first point forecast should be made at the 73th observation, for example, \code{w_size} should be set to be 72.}
  \item{window}{character, indicating the forecasting scheme to be applied. Options include \code{"recursive"}, \code{"rolling"}, and \code{"fixed"}.}
  \item{pred}{numeric, indicating the number of predicators being considered in the Adaptive Elastic Net model. Default is set to be equal to the number of observations.}
  \item{factor.num}{numeric, indicating the number of dynamic factors to be extracted from the covariates in the Dynamic Factor Model. Default is \code{3}.}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{Forecasts }{data matrix, containing the point forecasts, realized values, forecast errors, signs of the forecasts and realized values, and success in predicting the signs.}
  \item{MSE }{numeric, mean squred error of the point forecasts.}
  \item{SRatio}{numeric, success ratio of the point forecasts. Success is claimed when the point forecasts and realized values have the same sign.}
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
