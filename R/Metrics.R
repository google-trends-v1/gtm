Metrics<-function(pred=NULL, true=NULL, h=0){
  if(is.null(pred)|is.null(true)){
    stop("Arguments 'pred' and 'true' cannot be ommitted.")
  }
  forecasts<-data.frame(Forecasts=pred, Realized=true)
  forecasts$Errors<-forecasts$Realized-forecasts$Forecasts
  mse<-mean(na.omit(forecasts$Errors)^2)
  if(h==0|h==1){
    forecasts$Success <- ifelse(sign(forecasts$Forecasts) == sign(forecasts$Realized), 1, 0)
    success_ratio = sum(forecasts$Success)/nrow(forecasts)
    forecasts$Forecasted_Direction <- sign(forecasts$Forecasts)
    forecasts$Forecasted_Direction <- ifelse(sign(forecasts$Forecasted_Direction) == 1, 1, 0)
    forecasts$True_Direction <- sign(forecasts$Realized)
    forecasts$True_Direction <- ifelse(sign(forecasts$True_Direction) == 1, 1, 0)
  }else{
    foredir<-vector()
    truedir<-vector()
    foredir[1:(h-1)]<-NA
    truedir[1:(h-1)]<-NA
    for(i in h:length(forecasts$Forecasts)){
      foredir[i]<-cumsum(forecasts$Forecasts[(i+1-h):i])[h]
      truedir[i]<-cumsum(forecasts$Realized[(i+1-h):i])[h]
    }
    forecasts$Success<-ifelse(sign(foredir)==sign(truedir), 1, 0)
    success_ratio = sum(na.omit(forecasts$Success))/nrow(na.omit(forecasts))
    forecasts$Forecasted_Direction <- ifelse(sign(foredir) == 1, 1, 0)
    forecasts$True_Direction <- ifelse(sign(truedir) == 1, 1, 0)
  }
  results<-list(Forecasts=forecasts, MSE=mse, SRatio=success_ratio)
  return(results)
}
