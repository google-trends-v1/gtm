maeforecast.simplified<-function(data=NULL, model="ar", w_size=NULL, window="recursive", y.index=1, h=0, ...){
  if(model %in% c("ar", "lasso", "postlasso", "ridge",
                  "alasso", "postalasso", "postnet",
                  "dfm", "dfm2", "rf", "rw")==FALSE){
    stop("Unsupported model type. Refer to help(maeforecast) for a list of supported models.")
  }
  FUN<-paste("maeforecast.", model, sep="")
  FUN<-match.fun(FUN)
  results<-FUN(data=data, w_size=w_size, window=window, y.index=y.index, h=h, ...)$Forecasts$Forecasts
  return(results)
}
