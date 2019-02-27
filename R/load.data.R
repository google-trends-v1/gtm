load.data<-function(pattern=NULL, path=getwd(), merge=FALSE){
  setwd(path)
  list<-list.files(path)
  if(is.null(pattern)){
    list<-paste(path, "/", list, sep="")
  }else{
    list<-list[grep(pattern=pattern, x=list)]
    list<-paste(path, "/", list, sep="")
  }
  data.list<-list()
  df.names<-paste("df", 1:length(list), sep="")
  for(i in 1:length(list)){
    data<-read.csv(list[i])
    data.list<-append(data.list, list(assign(df.names[i], data)))
  }
  if(merge==FALSE){
    names(data.list)<-df.names
    return(data.list)
  }else{
    data.merged<-data.list[[1]]
    critical.len<-length(data.list[[1]][,1])
    for(i in 2:length(data.list)){
      data.merged<-merge(data.merged, data.list[[i]])
      if(length(data.merged[,1])<critical.len){
        return(cat("Error happens at ", list[i], "."))
      }
    }
    return(data.merged)
  }
}
