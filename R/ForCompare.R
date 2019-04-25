ForCompare<-function(..., benchmark.index=NULL, test=c("weighted", "binary"), h=1){
  model.list<-list(...)
  MSE<-vector()
  SRatio<-vector()
  if(class(benchmark.index)=="integer"){
    DMW<-vector()
    MSERatio<-vector()
    e1=as.numeric(model.list[[benchmark.index]]$Forecasts$Errors)
    MSEbench<-model.list[[benchmark.index]]$MSE
    for(i in 1:length(model.list)){
      MSE[i]<-model.list[[i]]$MSE
      SRatio[i]<-model.list[[i]]$SRatio
      MSERatio[i]<-MSE[i]/MSEbench
      if(i==as.numeric(benchmark.index)){
        DMW[i]<-NA
      }else{
        DMW[i]<-forecast::dm.test(e1=e1, e2=model.list[[i]]$Forecasts$Errors, "greater", h=h)$p.value
      }
    }
    table<-data.frame(MSE=MSE, SRatio=SRatio, MSERatio=MSERatio, DMW=DMW)
  }else if(is.null(benchmark.index)){
    for(i in 1:length(model.list)){
      MSE[i]<-model.list[[i]]$MSE
      SRatio[i]<-model.list[[i]]$SRatio
    }
    table<-data.frame(MSE=MSE, SRatio=SRatio)
  }else{
    stop("The argument benchmark.index should either be omitted or of the classe 'integer'. ")
  }
  if("weighted"%in%test){
    weighted.p<-vector()
    for(i in 1:length(model.list)){
      weighted.p[i]<-Directional_NW(forecasts=model.list[[i]], p=1, weighted=T)$pvalue
    }
    table$Weighted<-weighted.p
  }
  if("binary"%in%test){
    unweighted.p<-vector()
    for(i in 1:length(model.list)){
      unweighted.p[i]<-Directional_NW(forecasts=model.list[[i]], p=1, weighted=F)$pvalue
    }
    table$Unweighted<-unweighted.p
  }
  return(table)
}
