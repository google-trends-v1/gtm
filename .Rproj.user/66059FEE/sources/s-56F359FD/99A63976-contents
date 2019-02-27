library(gtrendsR)

get.trends<-function(queries=NULL, geo="US", time="all", path=getwd()){
  for(i in 1:length(queries)){
    keyword=queries[i]
    data=gtrends(keyword=keyword, geo=geo, time=time)[[1]][,1:4]
    if(is.null(data)){
      next
    }
    colnames(data)[2]<-data[2,3]
    data<-data[,1:2]
    file.name=paste(path, "/", keyword, ".csv", sep="")
    write.csv(x=data, file=file.name, row.names=F)
  }
}


