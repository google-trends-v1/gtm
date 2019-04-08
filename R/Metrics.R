Metrics<-function(pred=NULL, true=NULL){
  if(is.null(pred)|is.null(true)){
    stop("Arguments 'pred' and 'true' cannot be ommitted.")
  }
  forecasts<-data.frame(Forecasts=pred, Realized=true)
  forecasts$Errors<-forecasts$Forecasts-forecasts$Realized
  mse<-mean(na.omit(forecasts$Errors)^2)
  forecasts$Success <- ifelse(sign(forecasts$Forecasts) == sign(forecasts$Realized), 1, 0)
  success_ratio = sum(forecasts$Success)/nrow(forecasts)
  forecasts$Forecasted_Direction <- sign(forecasts$Forecasts)
  forecasts$Forecasted_Direction <- ifelse(sign(forecasts$Forecasted_Direction) == 1, 1, 0)
  forecasts$True_Direction <- sign(forecasts$Realized)
  forecasts$True_Direction <- ifelse(sign(forecasts$True_Direction) == 1, 1, 0)
  results<-list(Forecasts=forecasts, MSE=mse, SRatio=success_ratio)
  class(results)<-"Maeforecast"
  return(results)
}
