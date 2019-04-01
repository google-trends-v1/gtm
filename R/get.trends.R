get.trends<-function(queries=NULL, geo="US", time="all", path=getwd()){
  suppressMessages(require(gtrendsR))
  if(is.null(queries)|class(queries)!="character"){
    stop("Queries have to be provided as a vector of characters.")
  }else{
    num.queries=as.numeric(length(queries))
    num.files=0
    failures=vector()
    for(i in 1:length(queries)){
    keyword=queries[i]
    if(class(try(df<-gtrends(keyword=keyword, geo=geo, time=time)[[1]][,1:4], T))[1]!="try-error"){
      data<-df
    }else{
      keyword<-paste("'", keyword, "'", " ", sep="")
      failures<-append(failures, keyword)
      next
    }
    if(is.null(data)){
      keyword<-paste("'", keyword, "'", " ", sep="")
      failures<-append(failures, keyword)
      next
    }
    colnames(data)[2]<-attr(factor(data$keyword), "level")
    data<-data[,1:2]
    file.name=paste(path, "/", keyword, ".csv", sep="")
    write.csv(x=data, file=file.name, row.names=F)
    num.files=num.files+1
  }
  num.omit=num.queries-num.files
  cat(num.queries, " queries submitted in total.", "\n",num.files, " queries processed and saved in ", "'", path, "'", ".", "\n", num.omit, " queries omitted. ", "\n", sep="")
  if(num.omit!=0){
    cat("These are: ")
    cat(failures)
  }
  }
}



