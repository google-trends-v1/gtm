Directional_NW<-function(forecasts=NULL, p=1, weighted=TRUE){
  if(class(forecasts)!="Maeforecast"){
    stop("Argument forecasts has to be of class 'Maeforecast'.")
  }
  nw <- function(y,qn){
    T <- length(y)
    ybar <- rep(1,T) * ((sum(y))/T)
    dy <- y-ybar
    G0 <- t(dy) %*% dy/T
    for (j in 1:qn){
      gamma <- t(dy[(j+1):T]) %*% dy[1:T-j]/T
      G0 <- G0+(gamma+t(gamma))*(1-abs(j/(qn+1)))
    }
    return(as.numeric(G0))
  }

  if(weighted){
    weighted.dir<-forecasts$Forecasts$Forecasted_Direction*forecasts$Forecasts$Realized
    P <- length(weighted.dir)
    varfroll.adj <- nw(y=weighted.dir, qn=p)
    CW <- sqrt(P)*(mean(weighted.dir))/sqrt(varfroll.adj)
    pv <- 1-pnorm(CW,0,1)
    results=list(test.statistic=CW, pvalue=pv)
  }else{
    For<-forecasts$Forecasts$Forecasted_Direction-mean(forecasts$Forecasts$Forecasted_Direction)
    True<-forecasts$Forecasts$True_Direction-mean(forecasts$Forecasts$True_Direction)
    unweighted.dir<-For*True
    P<-length(unweighted.dir)
    varfroll.adj<-nw(y=unweighted.dir, qn=p)
    CW<-sqrt(P)*(mean(unweighted.dir))/sqrt(varfroll.adj)
    pv<-1-pnorm(CW,0,1)
    results=list(test.statistic=CW, pvalue=pv)
  }
  return(results)
}
