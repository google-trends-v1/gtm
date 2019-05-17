\name{Bagging}
\alias{Bagging}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Bagging Forecasts for Multivariate Time Series
%%  ~~function to do ... ~~
}
\description{This function generates forecasts for multivariate time series based on bagging methods. The bootstrapping process is done using block bootstrapping, and one-step-ahead forecasts are made for each bootstrapping with the \code{maeforecast.simplified} function. Forecasts of the bootstrapped time series then are aggregated to produce a single set of forecasts.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Bagging(data, boot, sim="fixed", l=12L, norm=TRUE, n.sim=NROW(data),
        model="ar", w_size, window="recursive", y.index=1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{a data frame or a matrix; the first column should contain the time series variable for which the forecasts are to be made. Other columns should contain the covariates.
%%     ~~Describe \code{data} here~~
}
  \item{boot}{number of bootstrapped versions to generate.
%%     ~~Describe \code{boot} here~~
}

  \item{sim}{The type of simulation required to generate the replicate time series. Defualt is \code{"fixed"}. Refer to \code{\link{tsboot}} for a full list of possible inputs.
  }
  \item{l}{If sim is "fixed" then l is the fixed block length used in generating the replicate time series. If sim is "geom" then l is the mean of the geometric distribution used to generate the block lengths. l should be a positive integer less than the length of tseries. This argument is not required when sim is "model" but it is required for all other simulation types. Defualt is \code{12L}.}
  \item{endcorr}{A logical variable indicating whether end corrections are to be applied when sim is "fixed". When sim is "geom", endcorr is automatically set to TRUE; endcorr is not used when sim is "model" or "scramble".}
  \item{norm}{A logical argument indicating whether normal margins should be used for phase scrambling. If norm is FALSE then margins corresponding to the exact empirical margins are used.}
  \item{n.sim}{The length of the simulated time series. Defualt is the length of the original time series.}
  \item{model}{character, indicating which model should be used to make the forecasts. Default is an AR(1) model. Supported models include \code{"ar"}, \code{"lasso"}, \code{"postlasso"} (Post-Lasso), \code{"ridge"}, \code{"alasso"} (Adaptive Lasso), \code{"postalasso"} (Post-AdaptiveLasso), \code{"postnet"} (Post-Adaptive ElasticNet), \code{"rf"} (Random Forests), \code{"dfm"}, \code{"dfm2"} (Dynamic Factor Models). Note that a random walk model cannot be used in this function.
%%     ~~Describe \code{FUN} here~~
}

  \item{w_size}{numeric, indicating the index where the forecasting should begin. If the first point forecast should be made at the 73th observation, for example, \code{w_size} should be set to be 72.}
  \item{window}{character, indicating the forecasting scheme to be applied. Options include \code{"recursive"}, \code{"rolling"}, and \code{"fixed"}.}
  \item{y.index}{numeric, indicating the column position of the time series for which the forecasts are made (Y). Defualt is \code{1}.}

  \item{\dots}{extra arguments supported by the \code{\link{maeforecast}} function.
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
This function returns an object of class "MaeBagging" that contains the following components:
  \item{Forecasts }{data matrix, containing the point forecasts, realized values, forecast errors, signs of the forecasts and realized values, and success in predicting the signs.}
  \item{MSE }{numeric, mean squred error of the point forecasts.}
  \item{SRatio}{numeric, success ratio of the point forecasts. Success is claimed when the point forecasts and realized values have the same sign.}
  \item{Data}{the data as used in the model.}
  \item{Model}{some specifics about the model used. }
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
For more detailed description of the arguments used in the bootstrapping process, refer to \code{\link{tsboot}}.
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{maeforecast}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
