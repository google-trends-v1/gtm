t.rank<-function(data=NULL, y.index=1L, h=0){
  if(is.null(data)){
    stop("Have to provide values for data and w_size.")
  }

  suppressMessages(require(forecast))
  suppressMessages(require(lmtest))

  Data = as.matrix(data)
  X = matrix(Data[1:(dim(Data)[1]-h),-y.index],nrow = (dim(Data)[1]-h))
  Y = matrix(Data[(1+h):dim(Data)[1],y.index],nrow = (dim(Data)[1]-h))
  index<-1:dim(X)[2]
  z<-vector()
  for(i in index){
    test.mod<-forecast::Arima(ts(Y), xreg=X[,i], order=c(1,0,0), include.mean=F)
    z[i]<-abs(as.numeric(lmtest::coeftest(test.mod)[2,3]))
  }
  rank<-as.data.frame(cbind(index, z))
  colnames(rank)<-c("index", "z")
  rank<-rank[order(rank$z, decreasing=T),]
  row.names(rank)<-NULL
  colnames(rank)<-c("X Index", "Z Score")
  return(rank)
}
