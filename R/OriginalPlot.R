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

  suppressMessages(require(zoo))
  suppressMessages(require(ggplot2))

  if(class(.data)=="Maeforecast"){
    w_size=.data$Model$Size
    h=.data$Model$Horizon
    model=.data$Model$Model
    Date<-index(original)
    Value<-as.numeric(original)
    ForeValue<-vector()

    for(i in 1:(length(Date)-w_size)){
      ForeValue[i]<-exp(.data$Forecasts$Forecasts[i])*Value[(w_size+i+h)]
    }

    ForeValue<-zoo(ForeValue, order.by=Date[(w_size+2+h):length(Date)])
    Values<-merge(zoo(original, Date), ForeValue)
    names(Values)<-c("Realized", "Forecasts")

    title<-paste(model, " Out-of-Sample Forecasts vs. Original Time Series, h=", h, sep="")

    plot<-ggplot2::ggplot()+
      geom_line(aes(x=Date, y=as.numeric(Values$Realized)))+
      geom_line(aes(x=Date, y=as.numeric(Values$Forecasts)), color="red", linetype="longdash", na.rm=T)+
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
    ForeValue<-vector()
    ForeLower<-vector()
    ForeUpper<-vector()

    for(i in 1:(length(Date)-w_size-1)){
      ForeValue[i]<-exp(.data$Interval$Forecast[i])*Value[(w_size+i+h+1)]
      ForeLower[i]<-exp(.data$Interval$Lower[i])*Value[(w_size+i+h+1)]
      ForeUpper[i]<-exp(.data$Interval$Upper[i])*Value[(w_size+i+h+1)]
    }


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
