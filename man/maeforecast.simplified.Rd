\name{maeforecast.simplified}
\alias{maeforecast.simplified}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Economic Forecasting with High-Dimensional Data
%%  ~~function to do ... ~~
}
\description{This function is almost identical to the \code{\link{maeforecast}} function excpet that it only returns the out-of-sample point forecasts.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
maeforecast.simplified(data, model="ar", w_size, window="recursive",
            y.index=1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{a data frame or a matrix; the first column should contain the time series variable for which the forecasts are to be made. Other columns should contain the covariates.}
  \item{model}{character, indicating which model should be used to make the forecasts. Default is an AR(1) model. Note that \code{maeforecast(model="ar")} is equivalent to \code{maeforecast.ar()} and so on. See Details for the full list of supported models.}
  \item{w_size}{numeric, indicating the index where the forecasting should begin. If the first point forecast should be made at the 73th observation, for example, \code{w_size} should be set to be 72.}
  \item{window}{character, indicating the forecasting scheme to be applied. Options include \code{"recursive"}, \code{"rolling"}, and \code{"fixed"}.}
  \item{y.index}{numeric, indicating the column position of the time series for which the forecasts are made (Y). Defualt is \code{1}.}
  \item{\dots}{other arguments that may be used. Refer to \code{\link{maeforecast}} for a full list.}

}
\details{
Supported models include \code{"ar"}, \code{"lasso"}, \code{"postlasso"} (Post-Lasso), \code{"ridge"}, \code{"alasso"} (Adaptive Lasso), \code{"postalasso"} (Post-AdaptiveLasso), \code{"postnet"} (Post-Adaptive ElasticNet), \code{"rf"} (Random Forests), \code{"dfm"} & \code{"dfm2"} (Dynamic Factor Models).

%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  A vector of out-of-sample point forecasts.
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
\code{\link{maeforecast}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
