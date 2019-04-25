df.clean<-function(df=NULL, type="correlation", threshold=NULL, smaller=TRUE, index=NULL){
  rm<-c()
  if(is.null(index)){
    index=1:dim(df)[2]
  }else{
    index=index
  }
  if(type=="mean"){
    if(smaller){
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
  }else if(type=="correlation"){
    df.new<-df[, index]
    tmp <- cor(df.new)
    tmp[upper.tri(tmp)] <- 0
    diag(tmp) <- 0
    df.new <- df.new[,!apply(tmp,2,function(x) any(x>=threshold))]
    df.new <-cbind(df[,-index], df.new)
    return(df.new)
  }else{
    stop("Type not supported.")
  }
}
