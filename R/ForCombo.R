ForCombo<-function(..., type="mean"){
  suppressMessages(require(matrixStats))
  model.list<-list(...)
  for(i in 1:length(model.list)){
    if(!class(model.list[[i]])%in%c("Maeforecast", "MaeBagging")){
      stop(paste("Object number ", i, " is not of class 'Maeforecast' or 'MaeBagging'."))
    }
  }
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
  results<-Metrics(pred=pred, true=true, h=model.list[[1]]$Model$Horizon)
  if(type=="mean"){
    Model="Mean"
  }else if(type=="median"){
    Model="Median"
  }
  results$Model<-list(Model=Model, Window=model.list[[1]]$Model$Window,
                      Size=model.list[[1]]$Model$Size, Horizon=model.list[[1]]$Model$Horizon,
                      Preselection="Check individual models.")
  class(results)<-"Maeforecast"
  return(results)
}
