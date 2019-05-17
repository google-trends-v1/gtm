DirPlot<-function(forecasts=NULL, original=NULL, size.range=NULL, title=NULL, legend.position="bottom"){
  if(class(forecasts)!="Maeforecast"){
    stop("The object passed to the argument forecasts has to be of class 'Maeforecast'.")
  }
  if(class(original)=="numeric"){
    original<-ts(original)
  }
  if(!class(original) %in% c("ts", "zoo")){
    stop("The object passed to the argument original has to be of class 'ts' or 'zoo'.")
  }
  if(is.null(size.range)){
    size.range<-c(3, 7)
  }

  suppressMessages(require(ggplot2))
  options(warn=-1)

  w_size=forecasts$Model$Size
  h=forecasts$Model$Horizon
  model=forecasts$Model$Model
  Date<-index(original)
  Value<-as.numeric(original)
  ForeDir<-as.factor(as.character(forecasts$Forecasts$Forecasted_Direction))
  levels(ForeDir)<-c("Down", "Up")
  Success<-as.factor(as.character(forecasts$Forecasts$Success))
  levels(Success)<-c("No", "Yes")
  Points<-data.frame(ForecastDirection=ForeDir, Success=Success,
                     Date=Date[(w_size+2+h):length(Date)],
                     Value=Value[(w_size+2+h):length(Value)])
  Points<-na.omit(Points)
  plot<-ggplot2::ggplot(data=Points)+
    geom_line(mapping=aes(x=Date, y=Value))+
    geom_point(data=Points, mapping=aes(color=ForecastDirection, size=Success,
                                        x=Date, y=Value, shape=ForecastDirection))+
    scale_color_discrete(name="Forecasted Direction")+
    scale_shape_discrete(name="Forecasted Direction")+
    scale_size_discrete(name="Success", range=c(3, 7))+
    theme(legend.position=legend.position)

  if(is.null(title)){
    title<-paste(model, " Model Out-of-Sample Directional Forecasts, h=", h, sep="")
  }else{
    title=title
  }

  plot<-plot+ggtitle(title)
  return(plot)
  options(warn=0)
}
