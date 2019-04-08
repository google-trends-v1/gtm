ForCompare<-function(..., benchmark=NULL){
  model.list<-list(...)
  MSE<-vector()
  SRatio<-vector()
  for(i in 1:length(model.list)){
    MSE[i]<-model.list[[i]]$MSE
    SRatio[i]<-model.list[[i]]$SRatio
  }
  table<-data.frame(MSE=MSE, SRatio=SRatio)
  if(is.null(benchmark)){
    return(table)
  }else{
    if(class(benchmark)=="Maeforecast"){
      table$MSERatio<-table$MSE/benchmark$MSE
    }else if(class(benchmark)=="numeric"){
      table$MSERatio<-table$MSE/benchmark
    }else if(class(benchmark)=="integer"){
      table$MSERatio<-table$MSE/model.list[[benchmark]]$MSE
    }
    return(table)
  }
}
