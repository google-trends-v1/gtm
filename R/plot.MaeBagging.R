plot.MaeBagging<-function(forecasts=NULL, start=NULL, frequency='month', forecast.lab="Forecasts", true.lab="Realized", x.lab="Time", y.lab="Value", title=NULL){
  if(class(forecasts)!="MaeBagging"){
    stop("Argument forecasts has to be of class 'Maeforecast'.")
  }
  if(is.null(start)){
    dates=seq(from=1, by=1, length.out=length(forecasts$Forecasts$Forecasts))
  }else{
    dates=seq(from = as.Date(start), by=frequency, length.out=length(forecasts$Forecasts$Forecasts))
  }
  suppressMessages(require(reshape2))
  suppressMessages(require(ggplot2))

  df <- data.frame(date=dates, predicts=as.numeric(forecasts$Forecasts$Forecasts),
                   true=as.numeric(forecasts$Forecasts$Realized))
  colnames(df) <-c("Time", forecast.lab, true.lab)
  df <- reshape2::melt(df, id.vars = "Time", value.name="Value", variable.name="Variable")
  plot<-ggplot2::ggplot(data = df) +
    geom_line(aes(x = Time, y = Value, linetype = Variable, color = Variable))+
    xlab(x.lab)+ylab(y.lab)
  if(!is.null(title)){
    plot<-plot+ggtitle(title)
  }
  return(plot)
}
