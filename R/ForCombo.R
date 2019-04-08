ForCombo<-function(..., type="mean"){
  suppressMessages(require(matrixStats))
  model.list<-list(...)
  forecasts<-matrix(ncol=length(model.list), nrow=nrow(model.list[[1]]$Forecasts))
  for(i in 1:length(model.list)){
    forecasts[,i]<-model.list[[i]]$Forecasts$Forecasts
  }
  if(type=="mean"){
    pred<-rowMeans(forecasts)
  }else if(type=="median"){
    pred<-matrixStats::rowMedians(forecasts)
  }else{
    stop("Type not supported.")
  }
  true<-model.list[[1]]$Forecasts$Realized
  results<-Metrics(pred=pred, true=true)
  return(results)
}
