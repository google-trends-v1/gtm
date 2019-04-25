maeforecast.lasso<-function(data=NULL, w_size=NULL, window="recursive", h=0, standardize=TRUE, lambda=NULL, y.index=1){
  if(is.null(data)|is.null(w_size)){
    stop("Have to provide values for data and w_size.")
  }
  if(window %in% c("recursive", "rolling", "fixed")==F){
    stop("Unsupported forecasting sheme. Has to be either 'recursive', 'rolling' or 'fixed'.")
  }

  Data = data
  X = matrix(Data[1:(dim(Data)[1]-h),-y.index],nrow = (dim(Data)[1]-h))
  Y = matrix(Data[(1+h):dim(Data)[1],y.index],nrow = (dim(Data)[1]-h))
  w_size = w_size
  n_windows = dim(Y)[1] - w_size
  predicts<-c()
  trues<-c()
  errors<-c()
  covariates<-list()

  suppressMessages(require(glmnet))

  if(window=="recursive"){
    for(i in 1:n_windows){
      lasso.mod<-glmnet(x = X[1:(w_size + i - 1),], y = Y[1:(w_size + i - 1),],
                        alpha=1,
                        standardize=standardize)

      lasso_cv<-cv.glmnet(x = X[1:(w_size + i - 1),], y = Y[1:(w_size + i - 1),],
                          type.measure="mse",
                          nfold=10,
                          alpha=1,
                          standardize=standardize,
                          lambda=lambda)


      best_lasso_coef <- coef(lasso_cv, s=lasso_cv$lambda.1se)
      best_lasso_coef = as.numeric(best_lasso_coef) [-1]
      nonzero_index<-which(best_lasso_coef!=0)
      covariates[i]<-list(nonzero_index)

      lasso_predict<-predict(lasso.mod, s=lasso_cv$lambda.1se, newx=matrix(X[(w_size+i),], ncol=ncol(X)))
      y_real<-Y[w_size+i,]
      e<-y_real-lasso_predict

      predicts[i]<-lasso_predict
      trues[i]<-y_real
      errors[i]<-e
    }
  }else if(window=="rolling"){
    for(i in 1:n_windows){
      lasso.mod<-glmnet(x = X[i:(w_size + i - 1),], y = Y[i:(w_size + i - 1),],
                        alpha=1, standardize=standardize)

      lasso_cv<-cv.glmnet(x = X[i:(w_size + i - 1),], y = Y[i:(w_size + i - 1),],
                          type.measure="mse",
                          nfold=10,
                          alpha=1,
                          standardize=standardize,
                          lambda=lambda)

      best_lasso_coef <- coef(lasso_cv, s = lasso_cv$lambda.1se)
      best_lasso_coef = as.numeric(best_lasso_coef) [-1]
      nonzero_index<-which(best_lasso_coef!=0)
      covariates[i]<-list(nonzero_index)

      lasso_predict<-predict(lasso.mod, s=lasso_cv$lambda.1se, newx=matrix(X[(w_size+i),], ncol=ncol(X)))
      y_real<-Y[w_size+i,]
      e<-y_real-lasso_predict

      predicts[i]<-lasso_predict
      trues[i]<-y_real
      errors[i]<-e
    }
  }else{
    lasso.mod<-glmnet(x = X[1:w_size,], y = Y[1:w_size,],
                      alpha=1, standardize=standardize)

    lasso_cv<-cv.glmnet(x = X[1:w_size,], y = Y[1:w_size,],
                        type.measure="mse",
                        nfold=10,
                        alpha=1,
                        standardize=standardize,
                        lambda=lambda)

    best_lasso_coef <- coef(lasso_cv, s = lasso_cv$lambda.1se)
    best_lasso_coef = as.numeric(best_lasso_coef) [-1]
    nonzero_index<-which(best_lasso_coef!=0)
    covariates<-list(nonzero_index)

    for(i in 1:n_windows){
      lasso_predict<-predict(lasso.mod, s=lasso_cv$lambda.1se, newx=matrix(X[(w_size+i),], ncol=ncol(X)))
      y_real<-Y[w_size+i,]
      e<-y_real-lasso_predict

      predicts[i]<-lasso_predict
      trues[i]<-y_real
      errors[i]<-e
    }
  }
  results<-Metrics(pred=predicts, true=trues)
  results$Variables<-covariates
  return(results)
}












maeforecast.postlasso<-function(data=NULL, w_size=NULL, window="recursive", h=0, standardize=TRUE, lambda=NULL, y.index=1){
  if(is.null(data)|is.null(w_size)){
    stop("Have to provide values for data and w_size.")
  }
  if(window %in% c("recursive", "rolling", "fixed")==F){
    stop("Unsupported forecasting sheme. Has to be either 'recursive', 'rolling' or 'fixed'.")
  }

  Data = data
  X = matrix(Data[1:(dim(Data)[1]-h),-y.index],nrow = (dim(Data)[1]-h))
  Y = matrix(Data[(1+h):dim(Data)[1],y.index],nrow = (dim(Data)[1]-h))
  w_size = w_size
  n_windows = dim(Y)[1] - w_size
  predicts<-c()
  trues<-c()
  errors<-c()
  covariates<-list()

  suppressMessages(require(glmnet))
  suppressMessages(require(forecast))

  if(window=="recursive"){
    for(i in 1:n_windows){
      lasso.mod<-glmnet(x = X[1:(w_size + i - 1),], y = Y[1:(w_size + i - 1),],
                        alpha=1, standardize=standardize)

      lasso_cv<-cv.glmnet(x = X[1:(w_size + i - 1),], y = Y[1:(w_size + i - 1),],
                          type.measure="mse",
                          nfold=10,
                          alpha=1,
                          standardize=standardize,
                          lambda=lambda)


      best_lasso_coef <- coef(lasso_cv, s = lasso_cv$lambda.1se)
      best_lasso_coef = as.numeric(best_lasso_coef) [-1]
      nonzero_index<-which(best_lasso_coef!=0)
      covariates[i]<-list(nonzero_index)

      if (length(nonzero_index) !=0){
        nonzero_index <-as.vector(nonzero_index)
        xregs <- X[1:(w_size + i - 1), nonzero_index]
        newxregs <-t(matrix(X[w_size + i, nonzero_index]))

        Data_uni = ts(Y, frequency=12)
        AR_lasso <-Arima(Data_uni[1:(w_size + i - 1)], order=c(1,0,0), xreg =xregs)
        AR_lasso_predict <- predict(AR_lasso, newxreg = newxregs, n.ahead =1)$pred
        y_real =Data_uni[w_size + i]

        predicts[i] <- AR_lasso_predict
        trues[i] <- y_real
        e <- y_real - AR_lasso_predict
        errors[i] <- e
      } else {
        Data_uni = ts(Y, frequency=12)
        AR_lasso <-Arima(Data_uni[1:(w_size + i - 1)], order=c(1,0,0))
        AR_lasso_predict <- predict(AR_lasso, n.ahead =1)$pred
        y_real =Data_uni[w_size + i]

        predicts[i] <- AR_lasso_predict
        trues[i] <- y_real
        e <- y_real - AR_lasso_predict
        errors[i] <- e
      }
    }
  }else if(window=="rolling"){
    for(i in 1:n_windows){
      lasso.mod<-glmnet(x = X[i:(w_size + i - 1),], y = Y[i:(w_size + i - 1),],
                        alpha=1, standardize=standardize)

      lasso_cv<-cv.glmnet(x = X[i:(w_size + i - 1),], y = Y[i:(w_size + i - 1),],
                          type.measure="mse",
                          nfold=10,
                          alpha=1,
                          standardize=standardize,
                          lambda=lambda)

      best_lasso_coef <- coef(lasso_cv, s = lasso_cv$lambda.1se)
      best_lasso_coef = as.numeric(best_lasso_coef) [-1]
      nonzero_index<-which(best_lasso_coef!=0)
      covariates[i]<-list(nonzero_index)

      if (length(nonzero_index) !=0){
        nonzero_index <-as.vector(nonzero_index)
        xregs <- X[i:(w_size + i - 1), nonzero_index]
        newxregs <-t(matrix(X[w_size + i, nonzero_index]))

        Data_uni = ts(Y, frequency=12)
        AR_lasso <-Arima(Data_uni[i:(w_size + i - 1)], order=c(1,0,0), xreg =xregs)
        AR_lasso_predict <- predict(AR_lasso, newxreg = newxregs, n.ahead =1)$pred
        y_real =Data_uni[w_size + i]

        predicts[i] <- AR_lasso_predict
        trues[i] <- y_real
        e <- y_real - AR_lasso_predict
        errors[i] <- e
      } else {
        Data_uni = ts(Y, frequency=12)
        AR_lasso <-Arima(Data_uni[i:(w_size + i - 1)], order=c(1,0,0))
        AR_lasso_predict <- predict(AR_lasso, n.ahead =1)$pred
        y_real =Data_uni[w_size + i]

        predicts[i] <- AR_lasso_predict
        trues[i] <- y_real
        e <- y_real - AR_lasso_predict
        errors[i] <- e
      }
    }
  }else{
    lasso.mod<-glmnet(x = X[1:w_size,], y = Y[1:w_size,],
                      alpha=1, standardize=standardize)

    lasso_cv<-cv.glmnet(x = X[1:w_size,], y = Y[1:w_size,],
                        type.measure="mse",
                        nfold=10,
                        alpha=1,
                        standardize=standardize,
                        lambda=lambda)

    best_lasso_coef <- coef(lasso_cv, s = lasso_cv$lambda.1se)
    best_lasso_coef = as.numeric(best_lasso_coef) [-1]
    nonzero_index<-which(best_lasso_coef!=0)
    covariates<-list(nonzero_index)

    if (length(nonzero_index) !=0){
      nonzero_index <-as.vector(nonzero_index)
      Data_uni = ts(Y, frequency=12)
      xregs <- X[1:w_size, nonzero_index]
      AR_lasso <-Arima(Data_uni[1:w_size], order=c(1,0,0), xreg =xregs)
      for(i in 1:n_windows){
        newxregs <-t(matrix(X[w_size + i, nonzero_index]))

        AR_lasso_predict <- predict(AR_lasso, newxreg = newxregs, n.ahead =1)$pred
        y_real =Data_uni[w_size + i]

        predicts[i] <- AR_lasso_predict
        trues[i] <- y_real
        e <- y_real - AR_lasso_predict
        errors[i] <- e
      }
    }else {
      Data_uni = ts(Y, frequency=12)
      AR_model <- Arima(Data_uni[1:w_size], order=c(1,0,0))
      predicts<-as.numeric(farima(model=AR_model, test=Data_uni[(w_size+1):length(Data)])$forecasts)
      trues<-Data_uni[(w_size+1):length(Data)]
      errors<-trues-predicts
    }
  }
  results<-Metrics(pred=predicts, true=trues)
  results$Variables<-covariates
  return(results)
}









maeforecast.ridge<-function(data=NULL, w_size=NULL, window="recursive", h=0, standardize=TRUE, lambda=NULL, y.index=1){
  if(is.null(data)|is.null(w_size)){
    stop("Have to provide values for data and w_size.")
  }
  if(window %in% c("recursive", "rolling", "fixed")==F){
    stop("Unsupported forecasting sheme. Has to be either 'recursive', 'rolling' or 'fixed'.")
  }

  Data = data
  X = matrix(Data[1:(dim(Data)[1]-h),-y.index],nrow = (dim(Data)[1]-h))
  Y = matrix(Data[(1+h):dim(Data)[1],y.index],nrow = (dim(Data)[1]-h))
  w_size = w_size
  n_windows = dim(Y)[1] - w_size
  predicts<-c()
  trues<-c()
  errors<-c()

  suppressMessages(require(glmnet))

  if(window=="recursive"){
    for(i in 1:n_windows){
      ridge.mod<-glmnet(x = X[1:(w_size + i - 1),], y = Y[1:(w_size + i - 1),],
                        alpha=0, standardize=standardize)

      ridge_cv<-cv.glmnet(x = X[1:(w_size + i - 1),], y = Y[1:(w_size + i - 1),],
                          type.measure="mse",
                          nfold=10,
                          alpha=0,
                          standardize=standardize,
                          lambda=lambda)

      ridge_predict<-predict(ridge.mod, s=ridge_cv$lambda.1se, newx=matrix(X[(w_size+i),], ncol=ncol(X)))
      y_real<-Y[w_size+i,]
      e<-y_real-ridge_predict

      predicts[i] <- ridge_predict
      trues[i] <- y_real
      errors[i] <- e
    }
  }else if (window=="rolling"){
    for(i in 1:n_windows){
      ridge.mod<-glmnet(x = X[i:(w_size + i - 1),], y = Y[i:(w_size + i - 1),],
                        alpha=0, standardize=standardize)

      ridge_cv<-cv.glmnet(x = X[i:(w_size + i - 1),], y = Y[i:(w_size + i - 1),],
                          type.measure="mse",
                          nfold=10,
                          alpha=0,
                          standardize=standardize,
                          lambda=lambda)

      ridge_predict<-predict(ridge.mod, s=ridge_cv$lambda.1se, newx=matrix(X[(w_size+i),], ncol=ncol(X)))
      y_real<-Y[w_size+i,]
      e<-y_real-ridge_predict

      predicts[i] <- ridge_predict
      trues[i] <- y_real
      errors[i] <- e
    }
  }else{
    ridge.mod<-glmnet(x = X[1:w_size,], y = Y[1:w_size,],
                      alpha=0, standardize=standardize)

    ridge_cv<-cv.glmnet(x = X[1:w_size,], y = Y[1:w_size,],
                        type.measure="mse",
                        nfold=10,
                        alpha=0,
                        standardize=standardize,
                        lambda=lambda)

    for(i in 1:n_windows){
      ridge_predict<-predict(ridge.mod, s=ridge_cv$lambda.1se, newx=matrix(X[(w_size+i),], ncol=ncol(X)))
      y_real<-Y[w_size+i,]
      e<-y_real-ridge_predict

      predicts[i] <- ridge_predict
      trues[i] <- y_real
      errors[i] <- e
    }
  }
  results<-Metrics(pred=predicts, true=trues)
  return(results)
}









maeforecast.alasso<-function(data=NULL, w_size=NULL, window="recursive", h=0, standardize=TRUE, lambda.ridge=NULL, lambda.lasso=NULL, y.index=1){
  if(is.null(data)|is.null(w_size)){
    stop("Have to provide values for data and w_size.")
  }
  if(window %in% c("recursive", "rolling", "fixed")==F){
    stop("Unsupported forecasting sheme. Has to be either 'recursive', 'rolling' or 'fixed'.")
  }

  Data = data
  X = matrix(Data[1:(dim(Data)[1]-h),-y.index],nrow = (dim(Data)[1]-h))
  Y = matrix(Data[(1+h):dim(Data)[1],y.index],nrow = (dim(Data)[1]-h))
  w_size = w_size
  n_windows = dim(Y)[1] - w_size
  predicts<-c()
  trues<-c()
  errors<-c()
  covariates<-list()

  suppressMessages(require(glmnet))

  if(window=="recursive"){
    for (i in 1:n_windows) {
      ridge_cv <- cv.glmnet(x = X[1:(w_size + i - 1),], y = Y[1:(w_size + i - 1),],
                             type.measure = "mse",
                             nfold = 10,
                             alpha = 0,
                            standardize=standardize,
                            lambda=lambda.ridge)
      best_ridge_coef <- as.numeric(coef(ridge_cv, s = ridge_cv$lambda.1se))[-1]

      alasso.mod <- glmnet(x = X[1:(w_size + i - 1),], y = Y[1:(w_size + i - 1),],
                            alpha = 1,
                            penalty.factor = 1/abs(best_ridge_coef),
                           standardize=standardize)

      alasso_cv <- cv.glmnet(x = X[1:(w_size + i - 1),], y = Y[1:(w_size + i - 1),],
                              type.measure = "mse",
                              nfold = 10,
                              alpha = 1,
                              penalty.factor = 1/abs(best_ridge_coef),
                              keep = TRUE,
                             standardize=standardize,
                             lambda=lambda.lasso)

      best_alasso_coef <- coef(alasso_cv, s = alasso_cv$lambda.1se)
      best_alasso_coef = as.numeric(best_alasso_coef) [-1]
      nonzero_index<-which(best_alasso_coef!=0)
      covariates[i]<-list(nonzero_index)

      alasso_predict<-predict(alasso.mod, s=alasso_cv$lambda.1se, newx=matrix(X[(w_size+i),], ncol=ncol(X)))
      y_real<-Y[w_size+i,]
      e<-y_real-alasso_predict

      predicts[i]<- alasso_predict
      trues[i]<- y_real
      errors[i]<- e
    }
  }else if(window=="rolling"){
    for (i in 1:n_windows) {
      ridge_cv <- cv.glmnet(x = X[i:(w_size + i - 1),], y = Y[i:(w_size + i - 1),],
                            type.measure = "mse",
                            nfold = 10,
                            alpha = 0,
                            standardize=standardize,
                            lambda=lambda.ridge)
      best_ridge_coef <- as.numeric(coef(ridge_cv, s = ridge_cv$lambda.1se))[-1]

      alasso.mod <- glmnet(x = X[i:(w_size + i - 1),], y = Y[i:(w_size + i - 1),],
                           alpha = 1,
                           penalty.factor = 1/abs(best_ridge_coef),
                           standardize=standardize)

      alasso_cv <- cv.glmnet(x = X[i:(w_size + i - 1),], y = Y[i:(w_size + i - 1),],
                             type.measure = "mse",
                             nfold = 10,
                             alpha = 1,
                             penalty.factor = 1/abs(best_ridge_coef),
                             keep = TRUE,
                             standardize=standardize,
                             lambda=lambda.lasso)

      best_alasso_coef <- coef(alasso_cv, s = alasso_cv$lambda.1se)
      best_alasso_coef = as.numeric(best_alasso_coef) [-1]
      nonzero_index<-which(best_alasso_coef!=0)
      covariates[i]<-list(nonzero_index)

      alasso_predict<-predict(alasso.mod, s=alasso_cv$lambda.1se, newx=matrix(X[(w_size+i),], ncol=ncol(X)))
      y_real<-Y[w_size+i,]
      e<-y_real-alasso_predict

      predicts[i]<- alasso_predict
      trues[i]<- y_real
      errors[i]<- e
    }
  }else{
      ridge_cv <- cv.glmnet(x = X[1:w_size,], y = Y[1:w_size,],
                            type.measure = "mse",
                            nfold = 10,
                            alpha = 0,
                            standardize=standardize,
                            lambda=lambda.ridge)
      best_ridge_coef <- as.numeric(coef(ridge_cv, s = ridge_cv$lambda.1se))[-1]

      alasso.mod <- glmnet(x = X[1:w_size,], y = Y[1:w_size,],
                           alpha = 1,
                           penalty.factor = 1/abs(best_ridge_coef),
                           standardize=standardize)

      alasso_cv <- cv.glmnet(x = X[1:w_size,], y = Y[1:w_size,],
                             type.measure = "mse",
                             nfold = 10,
                             alpha = 1,
                             penalty.factor = 1/abs(best_ridge_coef),
                             keep = TRUE,
                             standardize=standardize,
                             lambda=lambda.lasso)

      best_alasso_coef <- coef(alasso_cv, s = alasso_cv$lambda.1se)
      best_alasso_coef = as.numeric(best_alasso_coef) [-1]
      nonzero_index<-which(best_alasso_coef!=0)
      covariates<-list(nonzero_index)

      for(i in 1:n_windows){
        alasso_predict<-predict(alasso.mod, s=alasso_cv$lambda.1se, newx=matrix(X[(w_size+i),], ncol=ncol(X)))
        y_real<-Y[w_size+i,]
        e<-y_real-alasso_predict

        predicts[i]<- alasso_predict
        trues[i]<- y_real
        errors[i]<- e
      }
    }
  results<-Metrics(pred=predicts, true=trues)
  results$Variables<-covariates
  return(results)

}






















maeforecast.postalasso<-function(data=NULL, w_size=NULL, window="recursive", h=0, standardize=TRUE, lambda.ridge=NULL, lambda.lasso=NULL, y.index=1){
  if(is.null(data)|is.null(w_size)){
    stop("Have to provide values for data and w_size.")
  }
  if(window %in% c("recursive", "rolling", "fixed")==F){
    stop("Unsupported forecasting sheme. Has to be either 'recursive', 'rolling' or 'fixed'.")
  }

  Data = data
  X = matrix(Data[1:(dim(Data)[1]-h),-y.index],nrow = (dim(Data)[1]-h))
  Y = matrix(Data[(1+h):dim(Data)[1],y.index],nrow = (dim(Data)[1]-h))
  w_size = w_size
  n_windows = dim(Y)[1] - w_size
  predicts<-c()
  trues<-c()
  errors<-c()
  covariates<-list()

  suppressMessages(require(forecast))
  suppressMessages(require(glmnet))

  farima<-function(model=NULL, test=NULL){
    suppressMessages(require(stats))
    suppressMessages(require(forecast))
    suppressMessages(require(zoo))
    if(class(model)[1]!="ARIMA"){
      stop("The function only works with model fitted by Arima.")
    }else{
      forecasts<-vector()
      trainData<-as.numeric(model[["x"]])
      test.start<-length(trainData)+1
      testData<-as.numeric(test)
      fullData<-append(trainData, testData)
      test.end<-length(fullData)
      ar.order=as.numeric(model$call$order[[2]])
      ma.order=as.numeric(model$call$order[[4]])
      if(as.numeric(model$call$order[[3]])!=0){
        fullData<-diff(fullData, lag=1, diff=as.numeric(model$call$order[[3]]))
        test.start=test.start-as.numeric(model$call$order[[3]])
        test.end=test.end-as.numeric(model$call$order[[3]])
      }else{
        fullData<-fullData
        test.start=test.start
        test.end=test.end
      }
      epsilon<-rnorm(n=test.end, mean=0, sd=sqrt(as.numeric(model[["sigma2"]])))
      if(ar.order+ma.order!=as.numeric(length(model$coef))){
        intercept<-as.numeric(model$coef[length(model$coef)])
        if(ar.order==0){
          ar.coef=NULL
          ma.coef<-as.numeric(model$coef[1:1+ma.order])
        }else{
          ar.coef<-as.numeric(model$coef[1:ar.order])
          if(ma.order!=0){
            ma.coef<-as.numeric(model$coef[1+ar.order:ar.order+ma.order])
          }else{
            ma.coef<-NULL
          }
        }
      }else{
        intercept<-0
        if(ar.order==0){
          ar.coef=NULL
          ma.coef<-as.numeric(model$coef[1:1+ma.order])
        }else{
          ar.coef<-as.numeric(model$coef[1:ar.order])
          if(ma.order!=0){
            ma.coef<-as.numeric(model$coef[1+ar.order:ar.order+ma.order])
          }else{
            ma.coef<-NULL
          }
        }
      }

      for(i in test.start:test.end){
        if(is.null(ar.coef)){
          ar.value=0
        }else{
          ar.part<-vector()
          for(j in 1:length(ar.coef)){
            ar.part[j]<-ar.coef[j]*fullData[i-j]
          }
          ar.value<-sum(ar.part)
        }
        if(is.null(ma.coef)){
          ma.value=0
        }else{
          ma.part<-vector()
          for(j in 1:length(ma.coef)){
            ma.part[j]<-ma.coef[j]*epsilon[i-j]
          }
          ma.value<-sum(ma.part)
        }
        forecast<-ar.value+ma.value+intercept
        forecasts<-append(forecasts, forecast)
      }
      forData<-zoo(forecasts, test.start:test.end)
      rmse=sqrt(sum((forData-fullData[test.start:test.end])^2)/length(forData))
      forData<-as.list(forData)
      rmse<-as.list(rmse)
      list<-c(forData, rmse)
      names(list)<-c("forecasts", "rmse")
      return(list)
    }
  }

  if(window=="recursive"){
    for (i in 1:n_windows) {
      ridge_cv <- cv.glmnet(x = X[1:(w_size + i - 1),], y = Y[1:(w_size + i - 1),],
                            type.measure = "mse",
                            nfold = 10,
                            alpha = 0,
                            standardize=standardize,
                            lambda=lambda.ridge)
      best_ridge_coef <- as.numeric(coef(ridge_cv, s = ridge_cv$lambda.1se))[-1]

      alasso.mod <- glmnet(x = X[1:(w_size + i - 1),], y = Y[1:(w_size + i - 1),],
                           alpha = 1,
                           penalty.factor = 1/abs(best_ridge_coef),
                           standardize=standardize)

      alasso_cv <- cv.glmnet(x = X[1:(w_size + i - 1),], y = Y[1:(w_size + i - 1),],
                             type.measure = "mse",
                             nfold = 10,
                             alpha = 1,
                             penalty.factor = 1/abs(best_ridge_coef),
                             keep = TRUE,
                             standardize=standardize,
                             lambda=lambda.lasso)

      best_alasso_coef <- coef(alasso_cv, s = alasso_cv$lambda.1se)
      best_alasso_coef = as.numeric(best_alasso_coef) [-1]
      nonzero_index<-which(best_alasso_coef!=0)
      covariates[i]<-list(nonzero_index)

      if (length(nonzero_index) !=0){
        nonzero_index <-as.vector(nonzero_index)
        xregs <- X[1:(w_size + i - 1), nonzero_index]
        newxregs <-t(matrix(X[w_size + i, nonzero_index]))

        Data_uni = ts(Y, frequency=12)
        AR_alasso <-Arima(Data_uni[1:(w_size + i - 1)], order=c(1,0,0), xreg =xregs)
        AR_alasso_predict <- predict(AR_alasso, newxreg = newxregs, n.ahead =1)$pred
        y_real =Data_uni[w_size + i]

        predicts[i] <- AR_alasso_predict
        trues[i] <- y_real
        e <- y_real - AR_alasso_predict
        errors[i] <- e
      } else {
        Data_uni = ts(Y, frequency=12)
        AR_alasso <-Arima(Data_uni[1:(w_size + i - 1)], order=c(1,0,0))
        AR_alasso_predict <- predict(AR_alasso, n.ahead =1)$pred
        y_real =Data_uni[w_size + i]

        predicts[i] <- AR_alasso_predict
        trues[i] <- y_real
        e <- y_real - AR_alasso_predict
        errors[i] <- e
      }
    }
  }else if(window=="rolling"){
    for (i in 1:n_windows) {
      ridge_cv <- cv.glmnet(x = X[i:(w_size + i - 1),], y = Y[i:(w_size + i - 1),],
                            type.measure = "mse",
                            nfold = 10,
                            alpha = 0,
                            standardize=standardize,
                            lambda=lambda.ridge)
      best_ridge_coef <- as.numeric(coef(ridge_cv, s = ridge_cv$lambda.1se))[-1]

      alasso.mod <- glmnet(x = X[i:(w_size + i - 1),], y = Y[i:(w_size + i - 1),],
                           alpha = 1,
                           penalty.factor = 1/abs(best_ridge_coef),
                           standardize=standardize)

      alasso_cv <- cv.glmnet(x = X[i:(w_size + i - 1),], y = Y[i:(w_size + i - 1),],
                             type.measure = "mse",
                             nfold = 10,
                             alpha = 1,
                             penalty.factor = 1/abs(best_ridge_coef),
                             keep = TRUE,
                             standardize=standardize,
                             lambda-lambda.lasso)

      best_alasso_coef <- coef(alasso_cv, s = alasso_cv$lambda.1se)
      best_alasso_coef = as.numeric(best_alasso_coef) [-1]
      nonzero_index<-which(best_alasso_coef!=0)
      covariates[i]<-list(nonzero_index)

      if (length(nonzero_index) !=0){
        nonzero_index <-as.vector(nonzero_index)
        xregs <- X[i:(w_size + i - 1), nonzero_index]
        newxregs <-t(matrix(X[w_size + i, nonzero_index]))

        Data_uni = ts(Y, frequency=12)
        AR_alasso <-Arima(Data_uni[i:(w_size + i - 1)], order=c(1,0,0), xreg =xregs)
        AR_alasso_predict <- predict(AR_alasso, newxreg = newxregs, n.ahead =1)$pred
        y_real =Data_uni[w_size + i]

        predicts[i] <- AR_alasso_predict
        trues[i] <- y_real
        e <- y_real - AR_alasso_predict
        errors[i] <- e
      } else {
        Data_uni = ts(Y, frequency=12)
        AR_alasso <-Arima(Data_uni[i:(w_size + i - 1)], order=c(1,0,0))
        AR_alasso_predict <- predict(AR_alasso, n.ahead =1)$pred
        y_real =Data_uni[w_size + i]

        predicts[i] <- AR_alasso_predict
        trues[i] <- y_real
        e <- y_real - AR_alasso_predict
        errors[i] <- e
      }
    }
  }else{
    ridge_cv <- cv.glmnet(x = X[1:w_size,], y = Y[1:w_size,],
                          type.measure = "mse",
                          nfold = 10,
                          alpha = 0,
                          standardize=standardize,
                          lambda=lambda.ridge)
    best_ridge_coef <- as.numeric(coef(ridge_cv, s = ridge_cv$lambda.1se))[-1]

    alasso.mod <- glmnet(x = X[1:w_size,], y = Y[1:w_size,],
                         alpha = 1,
                         penalty.factor = 1/abs(best_ridge_coef),
                         standardize=standardize)

    alasso_cv <- cv.glmnet(x = X[1:w_size,], y = Y[1:w_size,],
                           type.measure = "mse",
                           nfold = 10,
                           alpha = 1,
                           penalty.factor = 1/abs(best_ridge_coef),
                           keep = TRUE,
                           standardize=standardize,
                           lambda=lambda.lasso)

    best_alasso_coef <- coef(alasso_cv, s = alasso_cv$lambda.1se)
    best_alasso_coef = as.numeric(best_alasso_coef) [-1]
    nonzero_index<-which(best_alasso_coef!=0)
    covariates<-list(nonzero_index)

    if (length(nonzero_index) !=0){
      nonzero_index <-as.vector(nonzero_index)
      Data_uni = ts(Y, frequency=12)
      xregs <- X[1:w_size, nonzero_index]
      AR_alasso <-Arima(Data_uni[1:w_size], order=c(1,0,0), xreg =xregs)
      for(i in 1:n_windows){
        newxregs <-t(matrix(X[w_size + i, nonzero_index]))

        AR_alasso_predict <- predict(AR_alasso, newxreg = newxregs, n.ahead =1)$pred
        y_real =Data_uni[w_size + i]

        predicts[i] <- AR_alasso_predict
        trues[i] <- y_real
        e <- y_real - AR_alasso_predict
        errors[i] <- e
      }
    }else {
      Data_uni = ts(Y, frequency=12)
      AR_model <- Arima(Data_uni[1:w_size], order=c(1,0,0))
      predicts<-as.numeric(farima(model=AR_model, test=Data_uni[(w_size+1):length(Data)])$forecasts)
      trues<-Data_uni[(w_size+1):length(Data)]
      errors<-trues-predicts
    }
  }
  results<-Metrics(pred=predicts, true=trues)
  results$Variables<-covariates
  return(results)
}
















maeforecast.postnet<-function(data=NULL, w_size=NULL, window="recursive", h=0, pred=NULL, standardize=TRUE, alphas=c(0.2, 0.8, 0.02), y.index=1){
  if(is.null(data)|is.null(w_size)){
    stop("Have to provide values for data and w_size.")
  }
  if(window %in% c("recursive", "rolling", "fixed")==F){
    stop("Unsupported forecasting sheme. Has to be either 'recursive', 'rolling' or 'fixed'.")
  }
  if(is.null(pred)){
    pred=w_size
  }else{
    pred=pred
  }


  Data = data
  X = matrix(Data[1:(dim(Data)[1]-h),-y.index],nrow = (dim(Data)[1]-h))
  Y = matrix(Data[(1+h):dim(Data)[1],y.index],nrow = (dim(Data)[1]-h))
  w_size = w_size
  n_windows = dim(Y)[1] - w_size
  predicts<-c()
  trues<-c()
  errors<-c()
  covariates<-list()

  suppressMessages(require(glm2))
  suppressMessages(require(msaenet))
  suppressMessages(require(forecast))
  suppressMessages(require(glmnet))

  farima<-function(model=NULL, test=NULL){
    suppressMessages(require(stats))
    suppressMessages(require(forecast))
    suppressMessages(require(zoo))
    if(class(model)[1]!="ARIMA"){
      stop("The function only works with model fitted by Arima.")
    }else{
      forecasts<-vector()
      trainData<-as.numeric(model[["x"]])
      test.start<-length(trainData)+1
      testData<-as.numeric(test)
      fullData<-append(trainData, testData)
      test.end<-length(fullData)
      ar.order=as.numeric(model$call$order[[2]])
      ma.order=as.numeric(model$call$order[[4]])
      if(as.numeric(model$call$order[[3]])!=0){
        fullData<-diff(fullData, lag=1, diff=as.numeric(model$call$order[[3]]))
        test.start=test.start-as.numeric(model$call$order[[3]])
        test.end=test.end-as.numeric(model$call$order[[3]])
      }else{
        fullData<-fullData
        test.start=test.start
        test.end=test.end
      }
      epsilon<-rnorm(n=test.end, mean=0, sd=sqrt(as.numeric(model[["sigma2"]])))
      if(ar.order+ma.order!=as.numeric(length(model$coef))){
        intercept<-as.numeric(model$coef[length(model$coef)])
        if(ar.order==0){
          ar.coef=NULL
          ma.coef<-as.numeric(model$coef[1:1+ma.order])
        }else{
          ar.coef<-as.numeric(model$coef[1:ar.order])
          if(ma.order!=0){
            ma.coef<-as.numeric(model$coef[1+ar.order:ar.order+ma.order])
          }else{
            ma.coef<-NULL
          }
        }
      }else{
        intercept<-0
        if(ar.order==0){
          ar.coef=NULL
          ma.coef<-as.numeric(model$coef[1:1+ma.order])
        }else{
          ar.coef<-as.numeric(model$coef[1:ar.order])
          if(ma.order!=0){
            ma.coef<-as.numeric(model$coef[1+ar.order:ar.order+ma.order])
          }else{
            ma.coef<-NULL
          }
        }
      }

      for(i in test.start:test.end){
        if(is.null(ar.coef)){
          ar.value=0
        }else{
          ar.part<-vector()
          for(j in 1:length(ar.coef)){
            ar.part[j]<-ar.coef[j]*fullData[i-j]
          }
          ar.value<-sum(ar.part)
        }
        if(is.null(ma.coef)){
          ma.value=0
        }else{
          ma.part<-vector()
          for(j in 1:length(ma.coef)){
            ma.part[j]<-ma.coef[j]*epsilon[i-j]
          }
          ma.value<-sum(ma.part)
        }
        forecast<-ar.value+ma.value+intercept
        forecasts<-append(forecasts, forecast)
      }
      forData<-zoo(forecasts, test.start:test.end)
      rmse=sqrt(sum((forData-fullData[test.start:test.end])^2)/length(forData))
      forData<-as.list(forData)
      rmse<-as.list(rmse)
      list<-c(forData, rmse)
      names(list)<-c("forecasts", "rmse")
      return(list)
    }
  }

  SIS.gaussian<-function (X, Y, pred, scale = F){
    if (scale == T) {
      X <- scale(X)
    }
    p = dim(X)[2]
    IndicesSIS <- rep(0, p)
    beta <- rep(0, p)
    for (jj in 1:p) {
      beta[jj] <- abs(glm2(Y ~ X[, jj], family = gaussian)$coefficients[2])
    }
    IndicesSIS <- sort(beta, index = TRUE, decreasing = TRUE)$ix
    Xc <- X[, IndicesSIS[1:pred]]
    Xsis<-list()
    Xsis$Xs <- Xc
    Xsis$Index<-IndicesSIS
    return(Xsis)
  }

  if(window=="recursive"){
    for (i in 1:n_windows) {

      Xsis <- SIS.gaussian(X[1:(w_size + i - 1),],
                           matrix(Y[1:(w_size + i - 1),1]),
                           pred=pred,
                           scale=standardize)

      aenet.fit <- aenet(Xsis$Xs, matrix(Y[1:(w_size + i - 1),1]),
                         alphas=alphas, seed = 10)

      nonzero_index = which(!coef(aenet.fit) == 0)
      covariates[i]<-list(nonzero_index)

      if (length(nonzero_index) !=0){
        nonzero_index <-as.vector(nonzero_index)
        xregs <- X[1:(w_size + i - 1), Xsis$Index[1:pred][nonzero_index]]
        newxregs <-t(matrix(X[w_size + i,Xsis$Index[1:pred][nonzero_index]]))

        Data_uni = ts(Y, frequency=12)
        AR_net <-Arima(Data_uni[1:(w_size + i - 1)], order=c(1,0,0), xreg =xregs)
        AR_net_predict <- predict(AR_net, newxreg = newxregs, n.ahead =1)$pred
        y_real =Data_uni[w_size + i]

        predicts[i] <- AR_net_predict
        trues[i] <- y_real
        e <- y_real - AR_net_predict
        errors[i] <- e
      } else {
        Data_uni = ts(Y, frequency=12)
        AR_net <-Arima(Data_uni[1:(w_size + i - 1)], order=c(1,0,0))
        AR_net_predict <- predict(AR_net, n.ahead =1)$pred
        y_real =Data_uni[w_size + i]

        predicts[i] <- AR_net_predict
        trues[i] <- y_real
        e <- y_real - AR_net_predict
        errors[i] <- e
      }
    }
  }else if(window=="rolling"){
    for (i in 1:n_windows) {

      Xsis <- SIS.gaussian(X[i:(w_size + i - 1),],
                           matrix(Y[i:(w_size + i - 1),1]),
                           pred=pred,
                           scale=standardize)

      aenet.fit <- aenet(Xsis$Xs, matrix(Y[i:(w_size + i - 1),1]),
                         alphas=alphas, seed = 10)

      nonzero_index = which(!coef(aenet.fit) == 0)
      covariates[i]<-list(nonzero_index)

      if (length(nonzero_index) !=0){
        nonzero_index <-as.vector(nonzero_index)
        xregs <- X[i:(w_size + i - 1), Xsis$Index[1:pred][nonzero_index]]
        newxregs <-t(matrix(X[w_size + i,Xsis$Index[1:pred][nonzero_index]]))

        Data_uni = ts(Y, frequency=12)
        AR_net <-Arima(Data_uni[i:(w_size + i - 1)], order=c(1,0,0), xreg =xregs)
        AR_net_predict <- predict(AR_net, newxreg = newxregs, n.ahead =1)$pred
        y_real =Data_uni[w_size + i]

        predicts[i] <- AR_net_predict
        trues[i] <- y_real
        e <- y_real - AR_net_predict
        errors[i] <- e
      } else {
        Data_uni = ts(Y, frequency=12)
        AR_net <-Arima(Data_uni[i:(w_size + i - 1)], order=c(1,0,0))
        AR_net_predict <- predict(AR_net, n.ahead =1)$pred
        y_real =Data_uni[w_size + i]

        predicts[i] <- AR_net_predict
        trues[i] <- y_real
        e <- y_real - AR_net_predict
        errors[i] <- e
      }
    }
  }else{
    Xsis <- SIS.gaussian(X[1:w_size,],
                         matrix(Y[1:w_size,1]),
                         pred=pred,
                         scale=standardize)

    aenet.fit <- aenet(Xsis$Xs, matrix(Y[1:w_size,1]),
                       alphas=alphas, seed = 10)

    nonzero_index = which(!coef(aenet.fit) == 0)
    covariates<-list(nonzero_index)


    if (length(nonzero_index) !=0){
      nonzero_index <-as.vector(nonzero_index)
      xregs <- X[1:w_size, Xsis$Index[1:pred][nonzero_index]]
      Data_uni = ts(Y, frequency=12)
      AR_net <-Arima(Data_uni[1:w_size], order=c(1,0,0), xreg =xregs)
      for(i in 1:n_windows){
        newxregs <-t(matrix(X[w_size + i,Xsis$Index[1:pred][nonzero_index]]))
        AR_net_predict <- predict(AR_net, newxreg = newxregs, n.ahead =1)$pred
        y_real =Data_uni[w_size + i]

        predicts[i] <- AR_net_predict
        trues[i] <- y_real
        e <- y_real - AR_net_predict
        errors[i] <- e
      }
    }else {
      Data_uni = ts(Y, frequency=12)

      AR_model <- Arima(Data_uni[1:w_size], order=c(1,0,0))
      predicts<-as.numeric(farima(model=AR_model, test=Data_uni[(w_size+1):length(Data)])$forecasts)
      trues<-Data_uni[(w_size+1):length(Data)]
      errors<-trues-predicts
    }
  }
  results<-Metrics(pred=predicts, true=trues)
  results$Variables<-covariates
  return(results)
}








maeforecast.dfm2<-function(data=NULL, w_size=NULL, h=0, window="recursive", factor.num=3, method="two-step", clustor.type="partitional", y.index=1){
  if(is.null(data)|is.null(w_size)){
    stop("Have to provide values for data and w_size.")
  }
  if(window %in% c("recursive", "rolling", "fixed")==F){
    stop("Unsupported forecasting sheme. Has to be either 'recursive', 'rolling' or 'fixed'.")
  }

  Data = data
  X = matrix(Data[1:(dim(Data)[1]-h),-y.index],nrow = (dim(Data)[1]-h))
  Y = matrix(Data[(1+h):dim(Data)[1],y.index],nrow = (dim(Data)[1]-h))
  w_size = w_size
  n_windows = dim(Y)[1] - w_size
  predicts<-c()
  trues<-c()
  errors<-c()

  suppressMessages(require(nowcasting))
  suppressMessages(require(forecast))

  if(window=="recursive"){
    for(i in 1:n_windows){
      allfactors<-clust.factor(X[1:(w_size + i),], fac.num=factor.num, method=method, clustor.type=clustor.type)
      factors<-allfactors[1:(w_size + i - 1),]
      newfactors<-matrix(allfactors[w_size+i,], nrow=1)

      AR_dfm<-Arima(Y[1:(w_size+i-1), ], order=c(1,0,0), xreg=factors)
      AR_dfm_predict<-as.numeric(predict(AR_dfm, newxreg=newfactors, n.ahead=1)$pred)

      y_real<-Y[w_size+i,]
      e<-y_real-AR_dfm_predict

      trues[i] <- y_real
      predicts[i] <- AR_dfm_predict
      errors[i] <- e
    }
  }else if(window=="rolling"){
    for(i in 1:n_windows){
      allfactors<-clust.factor(X[i:(w_size + i),], fac.num=factor.num, method=method, clustor.type=clustor.type)
      factors<-allfactors[1:w_size,]
      newfactors<-matrix(allfactors[w_size+1,], nrow=1)

      AR_dfm<-Arima(Y[i:(w_size+i-1), ], order=c(1,0,0), xreg=factors)
      AR_dfm_predict<-as.numeric(predict(AR_dfm, newxreg=newfactors, n.ahead=1)$pred)

      y_real<-Y[w_size+i,]
      e<-y_real-AR_dfm_predict

      trues[i] <- y_real
      predicts[i] <- AR_dfm_predict
      errors[i] <- e
    }
  }else{
    allfactors<-clust.factor(X, fac.num=factor.num, method=method, clustor.type=clustor.type)
    factors<-allfactors[1:w_size,]
    AR_dfm<-Arima(Y[1:w_size, ], order=c(1,0,0), xreg=factors)
    for(i in 1:n_windows){
      newfactors<-matrix(allfactors[w_size+i,], nrow=1)
      AR_dfm_predict<-as.numeric(predict(AR_dfm, newxreg=newfactors, n.ahead=1)$pred)

      y_real<-Y[w_size+i,]
      e<-y_real-AR_dfm_predict

      trues[i] <- y_real
      predicts[i] <- AR_dfm_predict
      errors[i] <- e
    }
  }
  results<-Metrics(pred=predicts, true=trues)
  return(results)
}

















maeforecast.dfm<-function(data=NULL, w_size=NULL, h=0, window="recursive", factor.num=3, y.index=1){
  if(is.null(data)|is.null(w_size)){
    stop("Have to provide values for data and w_size.")
  }
  if(window %in% c("recursive", "rolling", "fixed")==F){
    stop("Unsupported forecasting sheme. Has to be either 'recursive', 'rolling' or 'fixed'.")
  }
  if("dynfactoR" %in% rownames(installed.packages())==F){
    stop("Package dynfactorR is not installed. To install, call library(devtools) and then call install_github('rbagd/dynfactoR').")
  }

  Data = data
  X = matrix(Data[1:(dim(Data)[1]-h),-y.index],nrow = (dim(Data)[1]-h))
  Y = matrix(Data[(1+h):dim(Data)[1],y.index],nrow = (dim(Data)[1]-h))
  w_size = w_size
  n_windows = dim(Y)[1] - w_size
  predicts<-c()
  trues<-c()
  errors<-c()

  suppressMessages(require(nowcasting))
  suppressMessages(require(forecast))
  suppressMessages(require(dynfactoR))

  if(window=="recursive"){
    for(i in 1:n_windows){
      fac.mod<-dynfactoR::dfm(X[1:(w_size + i),], r=factor.num, p=1, q=factor.num)
      factors<-fac.mod$twostep[1:(w_size + i - 1),]
      newfactors<-matrix(fac.mod$twostep[w_size+i,], nrow=1)

      AR_dfm<-Arima(Y[1:(w_size+i-1), ], order=c(1,0,0), xreg=factors)
      AR_dfm_predict<-as.numeric(predict(AR_dfm, newxreg=newfactors, n.ahead=1)$pred)

      y_real<-Y[w_size+i,]
      e<-y_real-AR_dfm_predict

      trues[i] <- y_real
      predicts[i] <- AR_dfm_predict
      errors[i] <- e
    }
  }else if(window=="rolling"){
    for(i in 1:n_windows){
      fac.mod<-dynfactoR::dfm(X[i:(w_size + i),], r=factor.num, p=1, q=factor.num)
      factors<-fac.mod$twostep[1:w_size,]
      newfactors<-matrix(fac.mod$twostep[w_size+1,], nrow=1)

      AR_dfm<-Arima(Y[i:(w_size+i-1), ], order=c(1,0,0), xreg=factors)
      AR_dfm_predict<-as.numeric(predict(AR_dfm, newxreg=newfactors, n.ahead=1)$pred)

      y_real<-Y[w_size+i,]
      e<-y_real-AR_dfm_predict

      trues[i] <- y_real
      predicts[i] <- AR_dfm_predict
      errors[i] <- e
    }
  }else{
    fac.mod<-dynfactoR::dfm(X, r=factor.num, p=1, q=factor.num)
    factors<-fac.mod$twostep[1:w_size,]
    AR_dfm<-Arima(Y[1:w_size, ], order=c(1,0,0), xreg=factors)
    for(i in 1:n_windows){
      newfactors<-matrix(fac.mod$twostep[w_size+i,], nrow=1)
      AR_dfm_predict<-as.numeric(predict(AR_dfm, newxreg=newfactors, n.ahead=1)$pred)

      y_real<-Y[w_size+i,]
      e<-y_real-AR_dfm_predict

      trues[i] <- y_real
      predicts[i] <- AR_dfm_predict
      errors[i] <- e
    }
  }
  results<-Metrics(pred=predicts, true=trues)
  return(results)
}









maeforecast.rf<-function(data=NULL, w_size=NULL, h=0, window="recursive", ntree=500, replace=TRUE, y.index=1){
  if(is.null(data)|is.null(w_size)){
    stop("Have to provide values for data and w_size.")
  }
  if(window %in% c("recursive", "rolling", "fixed")==F){
    stop("Unsupported forecasting sheme. Has to be either 'recursive', 'rolling' or 'fixed'.")
  }

  Data = data
  X = matrix(Data[1:(dim(Data)[1]-h),-y.index],nrow = (dim(Data)[1]-h))
  Y = matrix(Data[(1+h):dim(Data)[1],y.index],nrow = (dim(Data)[1]-h))
  w_size = w_size
  n_windows = dim(Y)[1] - w_size
  predicts<-c()
  trues<-c()

  suppressMessages(require(randomForest))

  if(window=="recursive"){
    for(i in 1:n_windows){
      X.train<-X[1:(w_size+i-1),]
      Y.train<-Y[1:(w_size+i-1),]
      X.test<-X[(w_size + i),]
      Y.test<-Y[(w_size + i),]
      rf.mod<-randomForest(x=X.train, y=Y.train)
      rf_predict<-predict(rf.mod, newdata=X.test)
      predicts[i]<-rf_predict
      trues[i]<-Y.test
    }
  }else if(window=="rolling"){
    for(i in 1:n_windows){
      X.train<-X[i:(w_size+i-1),]
      Y.train<-Y[i:(w_size+i-1),]
      X.test<-X[(w_size + i),]
      Y.test<-Y[(w_size + i),]
      rf.mod<-randomForest(x=X.train, y=Y.train)
      rf_predict<-predict(rf.mod, newdata=X.test)
      predicts[i]<-rf_predict
      trues[i]<-Y.test
    }
  }else{
    X.train<-X[i:w_size,]
    Y.train<-Y[i:w_size,]
    for(i in 1:n_windows){
      X.test<-X[(w_size + i),]
      Y.test<-Y[(w_size + i),]
      rf.mod<-randomForest(x=X.train, y=Y.train)
      rf_predict<-predict(rf.mod, newdata=X.test)
      predicts[i]<-rf_predict
      trues[i]<-Y.test
    }
  }
  results<-Metrics(pred=predicts, true=trues)
  return(results)
}








maeforecast.ar<-function(data=NULL, w_size=NULL, window="recursive", y.index=1){
  if(is.null(data)|is.null(w_size)){
    stop("Have to provide values for data and w_size.")
  }
  if(window %in% c("recursive", "rolling", "fixed")==F){
    stop("Unsupported forecasting sheme. Has to be either 'recursive', 'rolling' or 'fixed'.")
  }

  Data = ts(data[,y.index], frequency=12)
  w_size = w_size
  n_windows = length(Data) - w_size
  predicts<-c()
  trues<-c()
  errors<-c()

  farima<-function(model=NULL, test=NULL){
    suppressMessages(require(stats))
    suppressMessages(require(forecast))
    suppressMessages(require(zoo))
    if(class(model)[1]!="ARIMA"){
      stop("The function only works with model fitted by Arima.")
    }else{
      forecasts<-vector()
      trainData<-as.numeric(model[["x"]])
      test.start<-length(trainData)+1
      testData<-as.numeric(test)
      fullData<-append(trainData, testData)
      test.end<-length(fullData)
      ar.order=as.numeric(model$call$order[[2]])
      ma.order=as.numeric(model$call$order[[4]])
      if(as.numeric(model$call$order[[3]])!=0){
        fullData<-diff(fullData, lag=1, diff=as.numeric(model$call$order[[3]]))
        test.start=test.start-as.numeric(model$call$order[[3]])
        test.end=test.end-as.numeric(model$call$order[[3]])
      }else{
        fullData<-fullData
        test.start=test.start
        test.end=test.end
      }
      epsilon<-rnorm(n=test.end, mean=0, sd=sqrt(as.numeric(model[["sigma2"]])))
      if(ar.order+ma.order!=as.numeric(length(model$coef))){
        intercept<-as.numeric(model$coef[length(model$coef)])
        if(ar.order==0){
          ar.coef=NULL
          ma.coef<-as.numeric(model$coef[1:1+ma.order])
        }else{
          ar.coef<-as.numeric(model$coef[1:ar.order])
          if(ma.order!=0){
            ma.coef<-as.numeric(model$coef[1+ar.order:ar.order+ma.order])
          }else{
            ma.coef<-NULL
          }
        }
      }else{
        intercept<-0
        if(ar.order==0){
          ar.coef=NULL
          ma.coef<-as.numeric(model$coef[1:1+ma.order])
        }else{
          ar.coef<-as.numeric(model$coef[1:ar.order])
          if(ma.order!=0){
            ma.coef<-as.numeric(model$coef[1+ar.order:ar.order+ma.order])
          }else{
            ma.coef<-NULL
          }
        }
      }

      for(i in test.start:test.end){
        if(is.null(ar.coef)){
          ar.value=0
        }else{
          ar.part<-vector()
          for(j in 1:length(ar.coef)){
            ar.part[j]<-ar.coef[j]*fullData[i-j]
          }
          ar.value<-sum(ar.part)
        }
        if(is.null(ma.coef)){
          ma.value=0
        }else{
          ma.part<-vector()
          for(j in 1:length(ma.coef)){
            ma.part[j]<-ma.coef[j]*epsilon[i-j]
          }
          ma.value<-sum(ma.part)
        }
        forecast<-ar.value+ma.value+intercept
        forecasts<-append(forecasts, forecast)
      }
      forData<-zoo(forecasts, test.start:test.end)
      rmse=sqrt(sum((forData-fullData[test.start:test.end])^2)/length(forData))
      forData<-as.list(forData)
      rmse<-as.list(rmse)
      list<-c(forData, rmse)
      names(list)<-c("forecasts", "rmse")
      return(list)
    }
  }

  if(window=="recursive"){
    for(i in 1:n_windows){
      AR_model <- arima(Data[1:(w_size + i - 1)], order=c(1,0,0))
      AR_predict <- predict(AR_model, 1)
      y_real =Data[w_size + i]

      predicts[i] <- AR_predict$pred
      trues[i] <- y_real
      e<-y_real-AR_predict$pred
      errors[i] <- e
    }
  }else if(window=="rolling"){
    for(i in 1:n_windows){
      AR_model <- arima(Data[i:(w_size + i - 1)], order=c(1,0,0))
      AR_predict <- predict(AR_model, 1)
      y_real =Data[w_size + i]

      predicts[i] <- AR_predict$pred
      trues[i] <- y_real
      e<-y_real-AR_predict$pred
      errors[i] <- e
    }
  }else{
    AR_model <- Arima(Data[1:w_size], order=c(1,0,0))
    predicts<-as.numeric(farima(model=AR_model, test=Data[(w_size+1):length(Data)])$forecasts)
    trues<-Data[(w_size+1):length(Data)]
    errors<-trues-predicts
  }
  results<-Metrics(pred=predicts, true=trues)
  return(results)
}



maeforecast<-function(data=NULL, model="ar", w_size=NULL, window="recursive", y.index=1, ...){
  if(model %in% c("ar", "lasso", "postlasso", "ridge",
                  "alasso", "postalasso", "postnet",
                  "dfm", "dfm2", "rf")==FALSE){
    stop("Unsupported model type. Refer to help(maeforecast) for a list of supported models.")
  }
  FUN<-paste("maeforecast.", model, sep="")
  FUN<-match.fun(FUN)
  results<-FUN(data=data, w_size=w_size, window=window, y.index=y.index, ...)
  return(results)
}






