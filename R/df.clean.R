df.clean<-function(df=NULL, threshold=NULL, smaller=TRUE, index=NULL){
  rm<-c()
  if(is.null(index)){
    index=1:dim(df)[2]
  }else{
    index=index
  }
  if(isTRUE(smaller)){
    for(i in index){
      if(mean(as.numeric(df[,i]))<=threshold){
        rm<-c(rm, i)
      }else{
        next
      }
    }
  }else{
    for(i in index){
      if(mean(as.numeric(df[,i]))>=threshold){
        rm<-c(rm, i)
      }else{
        next
      }
    }
  }
  if(length(rm)!=0){
    df.new=subset(df, select=-rm)
  }else{
    df.new=df
  }
  return(df.new)
}
