Bagging<-function(data=NULL, boot=NULL, model="ar", w_size=NULL, sim="fixed", l=12L, endcorr=TRUE, norm=TRUE, n.sim=NROW(data), window="recursive", y.index=1, h=0, ...){
  if(is.null(boot)){
    stop("Please indicate the number of bootstrapped versions to generate.")
  }
  suppressMessages(require(boot))
  results<-boot::tsboot(tseries=data, statistic=maeforecast.simplified, R=boot,
                  sim=sim, l=l, n.sim=n.sim, norm=norm, model=model, w_size=w_size, window=window,
                  endcorr=endcorr, y.index=y.index, h=h, ...)$t
  results<-t(results)
  forecasts<-rowMeans(results)
  trues<-data[(w_size+1+h):nrow(data),1]
  results<-Metrics(pred=forecasts, true=trues, h=h)
  model<-paste("Bagging.", model, sep="")
  results$Model<-list(Model=model, Window=window, Size=w_size, Horizon=h, Bootstrap=boot, Block=l, Simulation=sim, Length=n.sim)
  class(results)<-"Maeforecast"
  return(results)
}
