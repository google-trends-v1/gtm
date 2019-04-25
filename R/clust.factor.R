clust.factor<-function(data=NULL, fac.num=NULL, method="two-step", clustor.type="partitional"){
  if(is.null(data)){
    stop("Empty data.")
  }
  if(is.null(fac.num)){
    stop("Desired number of factors not specified.")
  }
  suppressMessages(require(dtwclust))
  if(method=="aggregation"){
    Partition<-dtwclust::tsclust(series=t(data), k=fac.num, type=clustor.type)
    clusters<-sort(unique(Partition@cluster))
    factors<-matrix(nrow=dim(data)[1], ncol=fac.num)
    for(i in clusters){
      members<-data.frame(data)[,which(Partition@cluster==i)]
      factors[,i]<-base::rowMeans(as.matrix(members, nrow=dim(data)[1]), na.rm=T)
    }
  }else if(method=="two-step"){
    if("dynfactoR" %in% rownames(installed.packages())==F){
      stop("Package dynfactorR is not installed. To install, call library(devtools) and then call install_github('rbagd/dynfactoR').")
    }
    suppressMessages(require(dynfactoR))
    Partition<-dtwclust::tsclust(series=t(data), k=fac.num, type=clustor.type)
    clusters<-sort(unique(Partition@cluster))
    factors<-matrix(nrow=dim(data)[1], ncol=fac.num)
    for(i in clusters){
      members<-data.frame(data)[,which(Partition@cluster==i)]
      factors[,i]<-dynfactoR::dfm(as.matrix(members, nrow=dim(data)[1]), r=1, q=1, p=1)$twostep
    }
  }else{
    stop("Mthod type not supported.")
  }
  return(factors)
}
