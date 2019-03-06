df.split<-function(dataframe=NULL, by=NULL){
  common<-dataframe[match(by, colnames(dataframe))]
  unique<-dataframe[which(colnames(dataframe)!=by)]
  df.names<-colnames(unique)
  data.list<-list()
  for(i in 1:length(colnames(unique))){
    data<-as.data.frame(cbind(common, unique[which(colnames(unique)==colnames(unique)[i])]))
    row.names(data)<-NULL
    data.list<-append(data.list, list(data))
  }
  names(data.list)<-df.names
  return(data.list)
}
