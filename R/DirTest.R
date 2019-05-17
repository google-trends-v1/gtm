Directional_NW<-function(forecasts=NULL, h=0, p=1, weighted=TRUE){
  if(!class(forecasts)%in%c("Maeforecast", "MaeBagging")){
    stop("Argument forecasts has to be of class 'Maeforecast' or 'MaeBagging'.")
  }
  nw <- function(y,qn){
    T <- length(y)
    ybar <- rep(1,T) * ((sum(y))/T)
    dy <- y-ybar
    G0 <- t(dy) %*% dy/T
    for (j in 1:qn){
      gamma <- t(dy[(j+1):T]) %*% dy[1:(T-j)]/T
      G0 <- G0+(gamma+t(gamma))*(1-abs(j/(qn+1)))
    }
    return(as.numeric(G0))
  }

  if(weighted){
    if(h==0 | h==1){
      For<-na.omit(forecasts$Forecasts$Forecasted_Direction)
      True<-na.omit(forecasts$Forecasts$Realized)
    }else{
      For<-na.omit(forecasts$Forecasts$Forecasted_Direction)
      True<-vector()
      True[1:(h-1)]<-NA
      for(i in h:length(forecasts$Forecasts$Realized)){
        True[i]<-cumsum(forecasts$Forecasts$Realized[(i+1-h):i])[h]
      }
      True<-na.omit(True)
    }
    weighted.dir<-as.numeric(For)*as.numeric(True)
    P <- length(weighted.dir)
    varfroll.adj <- nw(y=weighted.dir, qn=p)
    CW <- sqrt(P)*(mean(weighted.dir))/sqrt(varfroll.adj)
    pv <- 1-pnorm(CW,0,1)
    results=list(test.statistic=CW, pvalue=pv)
  }else{
    For<-na.omit(forecasts$Forecasts$Forecasted_Direction)-mean(na.omit(forecasts$Forecasts$Forecasted_Direction))
    True<-na.omit(forecasts$Forecasts$True_Direction)-mean(na.omit(forecasts$Forecasts$True_Direction))
    unweighted.dir<-as.numeric(For)*as.numeric(True)
    P<-length(unweighted.dir)
    varfroll.adj<-nw(y=unweighted.dir, qn=p)
    CW<-sqrt(P)*(mean(unweighted.dir))/sqrt(varfroll.adj)
    pv<-1-pnorm(CW,0,1)
    results=list(test.statistic=CW, pvalue=pv)
  }
  return(results)
}
