\name{maeforecast.dfm2}
\alias{maeforecast.dfm2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Economic Forecasting with High-Dimensional Data (Dynamic Factor Model with Clustering)
%%  ~~function to do ... ~~
}
\description{Like its counterpart \code{\link{maeforecast.dfm}}, this function makes out-of-sample forecasts based on a dynamic factor model. The difference is in the way dynamic factors are extracted. The \code{maeforecast.dfm2} function first implements a clustering process to the covariate time series based on the partitional method, and one dynamic factor is then extracted within each cluster either based on the two-step method proposed by Doz, Gianone & Reichlin (2011) or by aggregation.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
maeforecast.dfm2(data, w_size, window = "recursive",
                  factor.num = 3, method="two-step")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{a data frame or a matrix; the first column should contain the time series variable for which the forecasts are to be made. Other columns should contain the covariates.
%%     ~~Describe \code{data} here~~
}
  \item{w_size}{numeric, indicating the index where the forecasting should begin. If the first point forecast should be made at the 73th observation, for example, \code{w_size} should be set to be 72.
%%     ~~Describe \code{w_size} here~~
}
  \item{window}{character, indicating the forecasting scheme to be applied. Options include \code{"recursive"}, \code{"rolling"}, and \code{"fixed"}.
%%     ~~Describe \code{window} here~~
}
  \item{factor.num}{numeric, indicating the number of dynamic factors to be extracted from the covariates in the Dynamic Factor Model. Default is \code{3}.
%%     ~~Describe \code{factor.num} here~~
}
  \item{method}{character, indicating which method should be used to extract dynamic factors. If \code{"aggregation"}, the covariate time series are first clustered based on partitional method, and a simple aggreagtion is applied to each cluster. If \code{"two-step"}, one factor is extracted within each cluster based on Doz, Giannone & Reichlin (2011).}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{Forecasts }{data matrix, containing the point forecasts, realized values, forecast errors, signs of the forecasts and realized values, and success in predicting the signs.}
  \item{MSE }{numeric, mean squred error of the point forecasts.}
  \item{SRatio}{numeric, success ratio of the point forecasts. Success is claimed when the point forecasts and realized values have the same sign.}
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

\seealso{\code{\link{maeforecast.dfm}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line