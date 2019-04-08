Bagging<-function(data=NULL, boot=NULL, model="ar", w_size=NULL, sim="fixed", l=12L, endcorr=TRUE, norm=TRUE, n.sim=NROW(data), window="recursive", ...){
  if(is.null(boot)){
    stop("Please indicate the number of bootstrapped versions to generate.")
  }
  suppressMessages(require(boot))
  results<-tsboot(tseries=data, statistic=maeforecast.simplified, R=boot,
                  sim=sim, l=l, n.sim=n.sim, norm=norm, model=model, w_size=w_size, window=window,
                  endcorr=endcorr, ...)$t
  results<-t(results)
  forecasts<-rowMeans(results)
  trues<-data[(w_size+1):nrow(data),1]
  return(Metrics(pred=forecasts, true=trues))
}
