plot.RecMSE<-function(.data=NULL, type="ratio", title=NULL){
  if(type=="ratio"){
    x<-.data$MSE.ratio
    lab<-"MSE Ratio (Main/Benchmark)"
    intercept=1
  }else if(type=="diff"){
    x<-.data$MSE.diff
    lab<-"MSE Difference (Main-Benchmark)"
    intercept=0
  }else{
    stop("Type not supported. It should either be 'ratio' or 'diff'.")
  }
  if(is.null(title)){
    title<-"Recursive MSE"
  }

  suppressMessages(require(ggplot2))

  x<-data.frame(Index=1:length(x), Value=x)
  plot<-ggplot2::ggplot(x)+
    geom_line(aes(x=Index, y=Value))+
    xlab("Forecast Stopping Point")+
    ylab(lab)+
    ggtitle(title)+
    geom_hline(yintercept=intercept, color="red")
  return(plot)
}
