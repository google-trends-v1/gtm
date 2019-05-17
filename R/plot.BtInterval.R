plot.BtInterval<-function(.data=NULL, out.true=F){
  Results<-.data$Interval
  w_size=.data$Model$Size
  h=.data$Model$Horizon
  y.index=.data$Model$Index
  Model=.data$Model$Model
  Data<-.data$Data
  start=w_size+h+1
  interval <- structure(list(
    mean = ts(Results$Forecast, start=start),
    lower = ts(Results$Lower, start=start),
    upper = ts(Results$Upper, start=start),
    level=95), class="forecast")
  plot<-ggplot2::autoplot(ts(Data[,y.index][1:(w_size+h)]))+
    ggtitle(paste(Model, " 95% Out-of-Sample Predictive Interval, h=", h,  sep="")) +
    xlab("Time") + ylab("Value") +
    autolayer(interval, series="95% Interval")+
    theme(legend.position="none")
  if(out.true){
    true<-data.frame(Index=start:length(Data[,y.index]),
                     Value=Data[,y.index][start:length(Data[,y.index])])
    plot<-plot+
      geom_line(data=true, mapping=aes(x=Index, y=Value),
                color="blue", alpha=0.6, linetype="longdash")
  }
  return(plot)
}
