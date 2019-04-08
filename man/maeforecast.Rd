\name{maeforecast}
\alias{maeforecast}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Economic Forecasting with High-Dimensional Data
%%  ~~function to do ... ~~
}
\description{This function calculates one-step-ahead point forecasts with recursive, rolling or fixed windows in a data-rich environment. Supported machine learning algorithms include Lasso, Ridge, Adaptive Lasso, Adaptive Elastic Net, as well as AR(1) model as a potential benchmark. Dynamic factor models, Regression Trees and Random Forest algorithms are also supported. Out-of-sample forecasts are returned alongside mean squred errors and success ratios.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
maeforecast(data, model="ar", w_size, window="recursive", ...)

maeforecast.ar(data, w_size, window="recursive")
maeforecast.lasso(data, w_size, window="recursive", standardize=TRUE)
maeforecast.postlasso(data, w_size, window="recursive", standardize=TRUE)
maeforecast.ridge(data, w_size, window="recursive", standardize=TRUE)
maeforecast.alasso(data, w_size, window="recursive", standardize=TRUE)
maeforecast.postalasso(data, w_size, window="recursive", standardize=TRUE)
maeforecast.postnet(data, w_size, window="recursive", pred, standardize=TRUE)
maeforecast.rf(data, w_size, window="recursive", ntree=500, replace=TRUE)
maeforecast.dfm(data, w_size, window="recursive", factor.num)
maeforecast.dfm(data, w_size, window="recursive", factor.num, method="two-step")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{a data frame or a matrix; the first column should contain the time series variable for which the forecasts are to be made. Other columns should contain the covariates.}
  \item{model}{character, indicating which model should be used to make the forecasts. Default is an AR(1) model. Note that \code{maeforecast(model="ar")} is equivalent to \code{maeforecast.ar()} and so on. See Details for the full list of supported models.}
  \item{w_size}{numeric, indicating the index where the forecasting should begin. If the first point forecast should be made at the 73th observation, for example, \code{w_size} should be set to be 72.}
  \item{window}{character, indicating the forecasting scheme to be applied. Options include \code{"recursive"}, \code{"rolling"}, and \code{"fixed"}.}
  \item{pred}{numeric, indicating the number of predicators being considered in the Adaptive Elastic Net model. Default is set to be equal to the number of observations.}
  \item{factor.num}{numeric, indicating the number of dynamic factors to be extracted from the covariates in the Dynamic Factor Model. Default is \code{3}.}
  \item{standardize}{logical, indicating whether the data matrix should be scaled before the model is fitting, for the use of variable selection/shrinkage models. Default is \code{TRUE}.}
  \item{ntree}{numeric, number of trees to grow. Default is \code{500}.}
  \item{replace}{logical, indicating whether sampling of cases should be done with replacement. Default is \code{TRUE}.}
  \item{method}{character, indicating which method should be used to extract dynamic factors. If \code{"aggregation"}, the covariate time series are first clustered based on partitional method, and a simple aggreagtion is applied to each cluster. If \code{"two-step"}, one factor is extracted within each cluster based on Doz, Giannone & Reichlin (2011).}

}
\details{
Supported models include \code{"ar"}, \code{"lasso"}, \code{"postlasso"} (Post-Lasso), \code{"ridge"}, \code{"alasso"} (Adaptive Lasso), \code{"postalasso"} (Post-AdaptiveLasso), \code{"postnet"} (Post-Adaptive ElasticNet), \code{"rf"} (Random Forests), \code{"dfm"} & \code{"dfm2"} (Dynamic Factor Models).

For shrinkage models, the penalty factor lambda is selected automatically by 10-fold cross-validation.
%%  ~~ If necessary, more details than the description above ~~
}

\note{
%%  ~~further notes~~
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


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{maeforecast.ar}}, \code{\link{maeforecast.lasso}}, \code{\link{maeforecast.postlasso}}, \code{\link{maeforecast.ridge}}, \code{\link{maeforecast.alasso}}, \code{\link{maeforecast.postalasso}}, \code{\link{maeforecast.rf}}, \code{\link{maeforecast.postnet}}, \code{\link{maeforecast.dfm}}, \code{\link{maeforecast.dfm2}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
