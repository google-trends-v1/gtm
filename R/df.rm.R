df.rm<-function(df=NULL, index=NULL, rm="both", frequency=12){
  suppressMessages(require(stats))
  if(is.null(df)){
    return("Dataframe cannot be omitted.")
  }
  if(rm %in% c("seasonal", "trend", "both")==F){
    return("Component not removable.")
  }
  if(is.null(index)){
    index<-1:dim(df)[2]
  }
  df.new<-df
  for(i in index){
    components<-decompose(ts(df[,i], frequency=frequency))
    seasonal<-as.numeric(components$seasonal)
    trend<-as.numeric(components$trend)
    trend[is.na(trend)]<-0
    if(rm=="seasonal"){
      df.new[,i]<-df[,i]-seasonal
    }else if(rm=="trend"){
      df.new[,i]<-df[,i]-trend
    }else{
      df.new[,i]<-df[,i]-trend-seasonal
    }
  }
  return(df.new)
}
