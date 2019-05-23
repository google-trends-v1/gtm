ForCompare<-function(..., benchmark.index=NULL, test=c("weighted", "binary"), test.statistic=T, h=0){
  cw <- function(e.m1,e.m2,yf.m1,yf.m2){
    nw <- function(y,qn){
      T <- length(y)
      ybar <- rep(1,T) * ((sum(y))/T)
      dy <- y-ybar
      G0 <- t(dy) %*% dy/T
      for (j in 1:qn){
        gamma <- t(dy[(j+1):T]) %*% dy[1:T-j]/(T-1)
        G0 <- G0+(gamma+t(gamma))*(1-abs(j/qn))
      }
      return(as.numeric(G0))
    }
    P <- length(e.m1)
    froll.adj <- e.m1^2-(e.m2^2-(yf.m1-yf.m2)^2)
    varfroll.adj <- nw(froll.adj,1)
    CW <- sqrt(P)*(mean(froll.adj))/sqrt(varfroll.adj)
    pv <- 1-pnorm(CW,0,1)
    results=list(test=CW,pvalue=pv)
    return(results)
  }

  if(class(list(...)[[1]])=="list"){
    model.list<-list(...)[[1]]
  }else{
    model.list<-list(...)
  }

  for(i in 1:length(model.list)){
    if(!class(model.list[[i]])%in%c("Maeforecast", "MaeBagging")){
      stop(paste("Object number ", i, " is not of class 'Maeforecast' or 'MaeBagging'."))
    }
  }
  MSE<-vector()
  SRatio<-vector()
  names<-vector()
  for(i in 1:length(model.list)){
    names[i]<-model.list[[i]]$Model$Model
  }
  if(class(benchmark.index)=="integer"){
    DMWP<-vector()
    DMWT<-vector()
    CWP<-vector()
    CWT<-vector()
    MSERatio<-vector()
    e1=as.numeric(model.list[[benchmark.index]]$Forecasts$Errors)
    yf1=as.numeric(model.list[[benchmark.index]]$Forecasts$Forecasts)
    MSEbench<-model.list[[benchmark.index]]$MSE
    for(i in 1:length(model.list)){
      MSE[i]<-model.list[[i]]$MSE
      SRatio[i]<-model.list[[i]]$SRatio
      MSERatio[i]<-MSE[i]/MSEbench
      if(i==as.numeric(benchmark.index)){
        DMWP[i]<-NA
        DMWT[i]<-NA
        CWP[i]<-NA
        CWT[i]<-NA
      }else{

        if(class(dmw<-try(forecast::dm.test(e1=e1, e2=model.list[[i]]$Forecasts$Errors, "greater", h=1)$p.value, silent=T))=="try-error"){
          DMWP[i]<-NaN
          DMWT[i]<-NaN
        }else{
          DMWP[i]<-dmw
          DMWT[i]<-forecast::dm.test(e1=e1, e2=model.list[[i]]$Forecasts$Errors, "greater", h=1)$statistic

        }

        CWP[i]<-cw(e.m1=e1,
                  e.m2=model.list[[i]]$Forecasts$Errors,
                  yf.m1=yf1,
                  yf.m2=model.list[[i]]$Forecasts$Forecasts)$pvalue

        CWT[i]<-cw(e.m1=e1,
                   e.m2=model.list[[i]]$Forecasts$Errors,
                   yf.m1=yf1,
                   yf.m2=model.list[[i]]$Forecasts$Forecasts)$test
      }
    }
    if(test.statistic){
      table<-data.frame(Model=names, MSE=MSE, SRatio=SRatio, MSERatio=MSERatio,
                        DMWP=DMWP, DMWT=DMWT,
                        CWP=CWP, CWT=CWT)
    }else{
      table<-data.frame(Model=names, MSE=MSE, SRatio=SRatio, MSERatio=MSERatio,
                        DMWP=DMWP, CWP=CWP)
    }

  }else if(is.null(benchmark.index)){
    for(i in 1:length(model.list)){
      MSE[i]<-model.list[[i]]$MSE
      SRatio[i]<-model.list[[i]]$SRatio
    }
    table<-data.frame(Model=names, MSE=MSE, SRatio=SRatio)
  }else{
    stop("The argument benchmark.index should either be omitted or of the classe 'integer'. ")
  }
  if("weighted"%in%test){
    weighted.p<-vector()
    for(i in 1:length(model.list)){
      weighted.p[i]<-Directional_NW(forecasts=model.list[[i]], p=1, weighted=T, h=h)$pvalue
    }
    table$Weighted<-weighted.p
  }
  if("binary"%in%test){
    unweighted.p<-vector()
    for(i in 1:length(model.list)){
      unweighted.p[i]<-Directional_NW(forecasts=model.list[[i]], p=1, weighted=F, h=h)$pvalue
    }
    table$Unweighted<-unweighted.p
  }
  return(table)
}
