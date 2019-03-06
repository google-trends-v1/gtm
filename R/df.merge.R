df.merge<-function(df.list=NULL){
  data.merged<-df.list[[1]]
  critical.len<-length(row.names(df.list[[1]]))
  for(i in 2:length(df.list)){
    data.merged<-merge(data.merged, df.list[[i]])
    if(length(row.names(data.merged))<critical.len){
      return(cat("Error happens at ", names(df.list)[i], "."))
    }
  }
  return(data.merged)
}
