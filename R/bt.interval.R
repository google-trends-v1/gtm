bt.interval<-function(.data=NULL, boot=100, forecast="model"){
  if(class(.data)!="Maeforecast"){
    stop("The object passed to the argument .data has to be of class 'Maeforecast'.")
  }
  sim="fixed"
  l=12L
  endcorr=TRUE
  normal=TRUE
  Data<-.data$Data
  Model=.data$Model$Model
  window=.data$Model$Window
  w_size=.data$Model$Size
  h=.data$Model$Horizon
  y.index=.data$Model$Index
  Model.mega<-.data$Model

  suppressMessages(require(boot))
  suppressMessages(require(ggfortify))
  suppressMessages(require(forecast))

  if(Model=="Lasso"){
    results<-boot::tsboot(tseries=Data, statistic=maeforecast.simplified, R=boot,
                          sim=sim, l=l, n.sim=nrow(Data), norm=norm, model="lasso",
                          w_size=w_size, window=window, endcorr=endcorr, y.index=y.index,
                          h=h, standardize=.data$Model$Standardize, t.select=.data$Model$Preselection,
                          Lambda=.data$Model$Lambda, t.update=.data$Model$Update)$t
  }else if(Model=="Post Lasso"){
    results<-boot::tsboot(tseries=Data, statistic=maeforecast.simplified, R=boot,
                          sim=sim, l=l, n.sim=nrow(Data), norm=norm, model="postlasso",
                          w_size=w_size, window=window, endcorr=endcorr, y.index=y.index,
                          h=h, standardize=.data$Model$Standardize, t.select=.data$Model$Preselection,
                          Lambda=.data$Model$Lambda, t.update=.data$Model$Update)$t
  }else if(Model=="Ridge"){
    results<-boot::tsboot(tseries=Data, statistic=maeforecast.simplified, R=boot,
                          sim=sim, l=l, n.sim=nrow(Data), norm=norm, model="ridge",
                          w_size=w_size, window=window, endcorr=endcorr, y.index=y.index,
                          h=h, standardize=.data$Model$Standardize, t.select=.data$Model$Preselection,
                          Lambda=.data$Model$Lambda, t.update=.data$Model$Update)$t
  }else if(Model=="Adaptive Lasso"){
    results<-boot::tsboot(tseries=Data, statistic=maeforecast.simplified, R=boot,
                          sim=sim, l=l, n.sim=nrow(Data), norm=norm, model="alasso",
                          w_size=w_size, window=window, endcorr=endcorr, y.index=y.index,
                          h=h, standardize=.data$Model$Standardize, t.select=.data$Model$Preselection,
                          lambda.ridge=.data$Model$RidgeLambda, lambda.lasso=.data$Model$LassoLambda,
                          t.update=.data$Model$Update)$t
  }else if(Model=="Post Adaptive Lasso"){
    results<-boot::tsboot(tseries=Data, statistic=maeforecast.simplified, R=boot,
                          sim=sim, l=l, n.sim=nrow(Data), norm=norm, model="postalasso",
                          w_size=w_size, window=window, endcorr=endcorr, y.index=y.index,
                          h=h, standardize=.data$Model$Standardize, t.select=.data$Model$Preselection,
                          lambda.ridge=.data$Model$RidgeLambda, lambda.lasso=.data$Model$LassoLambda,
                          t.update=.data$Model$Update)$t
  }else if(Model=="Post Adaptive ElasticNet"){
    results<-boot::tsboot(tseries=Data, statistic=maeforecast.simplified, R=boot,
                          sim=sim, l=l, n.sim=nrow(Data), norm=norm, model="postnet",
                          w_size=w_size, window=window, endcorr=endcorr, y.index=y.index,
                          h=h, standardize=.data$Model$Standardize, t.select=.data$Model$Preselection,
                          pred=.data$Model$Predictors, alphas=.data$Model$Alphas,
                          t.update=.data$Model$Update)$t
  }else if(Model=="Dynamic Factor Model 2"){
    results<-boot::tsboot(tseries=Data, statistic=maeforecast.simplified, R=boot,
                          sim=sim, l=l, n.sim=nrow(Data), norm=norm, model="dfm2",
                          w_size=w_size, window=window, endcorr=endcorr, y.index=y.index,
                          h=h,t.select=.data$Model$Preselection, factor.num=.data$Model$Factors,
                          method=.data$Model$Method, clustor.type=.data$Model$Clustor,
                          t.update=.data$Model$Update)$t
  }else if(Model=="Dynamic Factor Model 1"){
    results<-boot::tsboot(tseries=Data, statistic=maeforecast.simplified, R=boot,
                          sim=sim, l=l, n.sim=nrow(Data), norm=norm, model="dfm",
                          w_size=w_size, window=window, endcorr=endcorr, y.index=y.index,
                          h=h,t.select=.data$Model$Preselection, factor.num=.data$Model$Factors,
                          t.update=.data$Model$Update)$t
  }else if(Model=="Random Forest"){
    results<-boot::tsboot(tseries=Data, statistic=maeforecast.simplified, R=boot,
                          sim=sim, l=l, n.sim=nrow(Data), norm=norm, model="rf",
                          w_size=w_size, window=window, endcorr=endcorr, y.index=y.index,
                          h=h,t.select=.data$Model$Preselection, ntree=.data$Model$Trees,
                          replace=.data$Model$Replace, t.update=.data$Model$Update)$t
  }else if(Model=="AR"){
    results<-boot::tsboot(tseries=Data, statistic=maeforecast.simplified, R=boot,
                          sim=sim, l=l, n.sim=nrow(Data), norm=norm, model="ar",
                          w_size=w_size, window=window, endcorr=endcorr, y.index=y.index,
                          h=h)$t
  }else if(Model=="Random Walk"){
    results<-boot::tsboot(tseries=Data, statistic=maeforecast.simplified, R=boot,
                          sim=sim, l=l, n.sim=nrow(Data), norm=norm, model="rw",
                          w_size=w_size, window=window, endcorr=endcorr, y.index=y.index,
                          h=h)$t
  }else if(Model=="ARIMAX"){
    results<-boot::tsboot(tseries=Data, statistic=maeforecast.simplified, R=boot,
                          sim=sim, l=l, n.sim=nrow(Data), norm=norm, model="arimax",
                          w_size=w_size, window=window, endcorr=endcorr, y.index=y.index,
                          h=h)$t
  }
  results<-t(results)
  forecasts<-.data$Forecasts$Forecasts
  lower<-apply(results, 1, quantile, prob=0.025)
  upper<-apply(results, 1, quantile, prob=0.975)

  if(forecast=="model"){
    forecasts=forecasts
  }else if(forecast=="mean"){
    forecasts=rowMeans(results)
  }

  Intervals<-data.frame(Lower=lower, Forecast=forecasts, Upper=upper)
  colnames(results)<-paste("Bootstrap.", seq(from=1, to=boot, by=1), sep="")
  Results<-list(Interval=Intervals, Model=Model.mega, Data=Data, Bootstrapped=results)
  class(Results)<-"BtInterval"
  return(Results)
}
