OriginalPlot<-function(.data=NULL, original=NULL){
  if(!class(.data)%in%c("Maeforecast", "BtInterval")){
    stop("The object passed to the argument forecasts has to be of class 'Maeforecast'.")
  }
  if(class(original)=="numeric"){
    original<-ts(original)
  }
  if(!class(original) %in% c("ts", "zoo")){
    stop("The object passed to the argument original has to be of class 'ts' or 'zoo'.")
  }

  if(class(.data)=="Maeforecast"){
    w_size=.data$Model$Size
    h=.data$Model$Horizon
    model=.data$Model$Model
    original=zoo(Panel$VALUE, order.by=Panel$Date)
    Date<-index(original)
    Value<-as.numeric(original)
    ForeValue<-exp(cumsum(.data$Forecasts$Forecasts))*Value[(w_size+1+h)]
    ForeValue<-zoo(ForeValue, order.by=Date[(w_size+2+h):length(Date)])
    Values<-merge(original, ForeValue)
    names(Values)<-c("Realized", "Forecasts")

    title<-paste(model, " Out-of-Sample Forecasts vs. Original Time Series, h=", h, sep="")

    plot<-ggplot2::ggplot(Values)+
      geom_line(aes(x=Date, y=as.numeric(Realized)))+
      geom_line(aes(x=Date, y=as.numeric(Forecasts)), color="red", linetype="longdash")+
      ylab("Values")+
      ggtitle(title)
  }else{
    w_size=.data$Model$Size
    h=.data$Model$Horizon
    model=.data$Model$Model
    boots=.data$Bootstrapped
    start=w_size+2+h
    original=zoo(Panel$VALUE, order.by=Panel$Date)
    Date<-index(original)
    Value<-as.numeric(original)
    ForeValue<-exp(cumsum(.data$Interval[,2]))*Value[(w_size+1+h)]
    Boot<-matrix(ncol=ncol(boots), nrow=nrow(boots))
    for(i in 1:ncol(boots)){
      Boot[,i]<-exp(cumsum(boots[,i]))*Value[(w_size+1+h)]
    }
    ForeLower<-apply(Boot, 1, quantile, prob=0.025)
    ForeUpper<-apply(Boot, 1, quantile, prob=0.975)
    interval <- structure(list(
      mean = ts(ForeValue, start=start),
      lower = ts(ForeLower, start=start),
      upper = ts(ForeUpper, start=start),
      level=95), class="forecast")

    plot<-ggplot2::autoplot(ts(as.numeric(original)))+
      ggtitle(paste(model, " Out-of-Sample Forecasts vs. Original Time Series, h=", h,  sep="")) +
      xlab("Time") + ylab("Value") +
      autolayer(interval, series="95% Interval")+
      theme(legend.position="none")
  }
  return(plot)
}
