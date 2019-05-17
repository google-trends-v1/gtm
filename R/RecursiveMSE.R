RecursiveMSE<-function(forecast.main=NULL, forecast.benchmark=NULL){
  if(class(forecast.main)!="Maeforecast"|class(forecast.benchmark)!="Maeforecast"){
    stop("Inputs must be of class 'Maeforecast'.")
  }
  if(length(forecast.main$Forecasts$Forecasts)!=length(forecast.benchmark$Forecasts$Forecasts)){
    stop("Lengths of forecasts differ.")
  }
  MSE.main<-vector()
  MSE.benchmark<-vector()
  for(i in 1:length(forecast.main$Forecasts$Forecasts)){
    MSE.main[i]<-mean(na.omit(forecast.main$Forecasts$Errors[1:i])^2)
    MSE.benchmark[i]<-mean(na.omit(forecast.benchmark$Forecasts$Errors[1:i])^2)
  }
  results<-data.frame(MSE.main=MSE.main, MSE.benchmark=MSE.benchmark)
  results$MSE.diff<-MSE.main-MSE.benchmark
  results$MSE.ratio<-MSE.main/MSE.benchmark
  class(results)<-"RecMSE"
  return(results)
}
