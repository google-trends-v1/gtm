maeforecast.lasso<-function(data=NULL, w_size=NULL, window="recursive"){
  if(is.null(data)|is.null(w_size)){
    return("Have to provide values for data and w_size")
  }
  if(window %in% c("recursive", "rolling", "fixed")==F){
    return("Unsupported forecasting sheme. Has to be either 'recursive', 'rolling' or 'fixed'.")
  }

  Data = data
  X = matrix(Data[,-1],nrow = dim(Data[,-1])[1])
  Y = matrix(Data[,1],nrow = dim(Data[,-1])[1])
  w_size = w_size
  n_windows = nrow(Data) - w_size
  predicts<-c()
  trues<-c()
  errors<-c()

  if(window=="recursive"){
    for(i in 1:n_windows){
      lasso.mod<-glmnet(x = X[1:(w_size + i - 1),], y = Y[1:(w_size + i - 1),],
                        alpha=1)

      lasso_cv<-cv.glmnet(x = X[1:(w_size + i - 1),], y = Y[1:(w_size + i - 1),],
                          type.measure="mse",
                          nfold=10,
                          alpha=1)

      lasso_predict<-predict(lasso.mod, s=lasso_cv$lambda.min, newx=matrix(X[(w_size+i),], ncol=ncol(X)))
      y_real<-Y[w_size+i,]
      e<-y_real-lasso_predict

      predicts[i]<-lasso_predict
      trues[i]<-y_real
      errors[i]<-e
    }
  }else if(window=="rolling"){
    for(i in 1:n_windows){
      lasso.mod<-glmnet(x = X[i:(w_size + i - 1),], y = Y[i:(w_size + i - 1),],
                        alpha=1)

      lasso_cv<-cv.glmnet(x = X[i:(w_size + i - 1),], y = Y[i:(w_size + i - 1),],
                          type.measure="mse",
                          nfold=10,
                          alpha=1)

      lasso_predict<-predict(lasso.mod, s=lasso_cv$lambda.min, newx=matrix(X[(w_size+i),], ncol=ncol(X)))
      y_real<-Y[w_size+i,]
      e<-y_real-lasso_predict

      predicts[i]<-lasso_predict
      trues[i]<-y_real
      errors[i]<-e
    }
  }else{
    lasso.mod<-glmnet(x = X[1:w_size,], y = Y[1:w_size,],
                      alpha=1)

    lasso_cv<-cv.glmnet(x = X[1:w_size,], y = Y[1:w_size,],
                        type.measure="mse",
                        nfold=10,
                        alpha=1)
    for(i in 1:n_windows){
      lasso_predict<-predict(lasso.mod, s=lasso_cv$lambda.min, newx=matrix(X[(w_size+i),], ncol=ncol(X)))
      y_real<-Y[w_size+i,]
      e<-y_real-lasso_predict

      predicts[i]<-lasso_predict
      trues[i]<-y_real
      errors[i]<-e
    }
  }

  mse<-mean(na.omit(errors)^2)

  forecasts<-data.frame(predicts, trues)
  colnames(forecasts) <-c('Forecasts','Realized')
  forecasts$Errors<-errors
  forecasts$Success <- ifelse(sign(forecasts$Forecasts) == sign(forecasts$Realized), 1, 0)
  success_ratio = sum(forecasts$Success)/nrow(forecasts)
  forecasts$Forecasted_Direction <- sign(forecasts$Forecasts)
  forecasts$Forecasted_Direction <- ifelse(sign(forecasts$Forecasted_Direction) == 1, 1, 0)
  forecasts$True_Direction <- sign(forecasts$Realized)
  forecasts$True_Direction <- ifelse(sign(forecasts$True_Direction) == 1, 1, 0)


  results<-list(Forecasts=forecasts, MSE=mse, SRatio=success_ratio)

  return(results)
}












maeforecast.arlasso<-function(data=NULL, w_size=NULL, window="recursive"){
  if(is.null(data)|is.null(w_size)){
    return("Have to provide values for data and w_size")
  }
  if(window %in% c("recursive", "rolling", "fixed")==F){
    return("Unsupported forecasting sheme. Has to be either 'recursive', 'rolling' or 'fixed'.")
  }

  Data = data
  X = matrix(Data[,-1],nrow = dim(Data[,-1])[1])
  Y = matrix(Data[,1],nrow = dim(Data[,-1])[1])
  w_size = w_size
  n_windows = nrow(Data) - w_size
  predicts<-c()
  trues<-c()
  errors<-c()

  if(window=="recursive"){
    for(i in 1:n_windows){
      lasso.mod<-glmnet(x = X[1:(w_size + i - 1),], y = Y[1:(w_size + i - 1),],
                        alpha=1)

      lasso_cv<-cv.glmnet(x = X[1:(w_size + i - 1),], y = Y[1:(w_size + i - 1),],
                          type.measure="mse",
                          nfold=10,
                          alpha=1)

      best_lasso_coef <- coef(lasso_cv, s = lasso_cv$lambda.min)
      best_lasso_coef = as.numeric(best_lasso_coef) [-1]
      nonzero_index<-which(best_lasso_coef!=0)

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
                        alpha=1)

      lasso_cv<-cv.glmnet(x = X[i:(w_size + i - 1),], y = Y[i:(w_size + i - 1),],
                          type.measure="mse",
                          nfold=10,
                          alpha=1)

      best_lasso_coef <- coef(lasso_cv, s = lasso_cv$lambda.min)
      best_lasso_coef = as.numeric(best_lasso_coef) [-1]
      nonzero_index<-which(best_lasso_coef!=0)

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
                      alpha=1)

    lasso_cv<-cv.glmnet(x = X[1:w_size,], y = Y[1:w_size,],
                        type.measure="mse",
                        nfold=10,
                        alpha=1)

    best_lasso_coef <- coef(lasso_cv, s = lasso_cv$lambda.min)
    best_lasso_coef = as.numeric(best_lasso_coef) [-1]
    nonzero_index<-which(best_lasso_coef!=0)

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
  mse<-mean(na.omit(errors)^2)

  forecasts<-data.frame(predicts, trues)
  colnames(forecasts) <-c('Forecasts','Realized')
  forecasts$Errors<-errors
  forecasts$Success <- ifelse(sign(forecasts$Forecasts) == sign(forecasts$Realized), 1, 0)
  success_ratio = sum(forecasts$Success)/nrow(forecasts)
  forecasts$Forecasted_Direction <- sign(forecasts$Forecasts)
  forecasts$Forecasted_Direction <- ifelse(sign(forecasts$Forecasted_Direction) == 1, 1, 0)
  forecasts$True_Direction <- sign(forecasts$Realized)
  forecasts$True_Direction <- ifelse(sign(forecasts$True_Direction) == 1, 1, 0)


  results<-list(Forecasts=forecasts, MSE=mse, SRatio=success_ratio)

  return(results)
}









maeforecast.ridge<-function(data=NULL, w_size=NULL, window="recursive"){
  if(is.null(data)|is.null(w_size)){
    return("Have to provide values for data and w_size")
  }
  if(window %in% c("recursive", "rolling", "fixed")==F){
    return("Unsupported forecasting sheme. Has to be either 'recursive', 'rolling' or 'fixed'.")
  }

  Data = data
  X = matrix(Data[,-1],nrow = dim(Data[,-1])[1])
  Y = matrix(Data[,1],nrow = dim(Data[,-1])[1])
  w_size = w_size
  n_windows = nrow(Data) - w_size
  predicts<-c()
  trues<-c()
  errors<-c()

  if(window=="recursive"){
    for(i in 1:n_windows){
      ridge.mod<-glmnet(x = X[1:(w_size + i - 1),], y = Y[1:(w_size + i - 1),],
                        alpha=0)

      ridge_cv<-cv.glmnet(x = X[1:(w_size + i - 1),], y = Y[1:(w_size + i - 1),],
                          type.measure="mse",
                          nfold=10,
                          alpha=0)


      ridge_predict<-predict(ridge.mod, s=ridge_cv$lambda.min, newx=matrix(X[(w_size+i),], ncol=ncol(X)))
      y_real<-Y[w_size+i,]
      e<-y_real-ridge_predict

      predicts[i] <- ridge_predict
      trues[i] <- y_real
      errors[i] <- e
    }
  }else if (window=="rolling"){
    for(i in 1:n_windows){
      ridge.mod<-glmnet(x = X[i:(w_size + i - 1),], y = Y[i:(w_size + i - 1),],
                        alpha=0)

      ridge_cv<-cv.glmnet(x = X[i:(w_size + i - 1),], y = Y[i:(w_size + i - 1),],
                          type.measure="mse",
                          nfold=10,
                          alpha=0)


      ridge_predict<-predict(ridge.mod, s=ridge_cv$lambda.min, newx=matrix(X[(w_size+i),], ncol=ncol(X)))
      y_real<-Y[w_size+i,]
      e<-y_real-ridge_predict

      predicts[i] <- ridge_predict
      trues[i] <- y_real
      errors[i] <- e
    }
  }else{
    ridge.mod<-glmnet(x = X[1:w_size,], y = Y[1:w_size,],
                      alpha=0)

    ridge_cv<-cv.glmnet(x = X[1:w_size,], y = Y[1:w_size,],
                        type.measure="mse",
                        nfold=10,
                        alpha=0)
    for(i in 1:n_windows){
      ridge_predict<-predict(ridge.mod, s=ridge_cv$lambda.min, newx=matrix(X[(w_size+i),], ncol=ncol(X)))
      y_real<-Y[w_size+i,]
      e<-y_real-ridge_predict

      predicts[i] <- ridge_predict
      trues[i] <- y_real
      errors[i] <- e
    }
  }
  mse<-mean(na.omit(errors)^2)

  forecasts<-data.frame(predicts, trues)
  colnames(forecasts) <-c('Forecasts','Realized')
  forecasts$Errors<-errors
  forecasts$Success <- ifelse(sign(forecasts$Forecasts) == sign(forecasts$Realized), 1, 0)
  success_ratio = sum(forecasts$Success)/nrow(forecasts)
  forecasts$Forecasted_Direction <- sign(forecasts$Forecasts)
  forecasts$Forecasted_Direction <- ifelse(sign(forecasts$Forecasted_Direction) == 1, 1, 0)
  forecasts$True_Direction <- sign(forecasts$Realized)
  forecasts$True_Direction <- ifelse(sign(forecasts$True_Direction) == 1, 1, 0)


  results<-list(Forecasts=forecasts, MSE=mse, SRatio=success_ratio)

  return(results)

}









maeforecast.alasso<-function(data=NULL, w_size=NULL, window="recursive"){
  if(is.null(data)|is.null(w_size)){
    return("Have to provide values for data and w_size")
  }
  if(window %in% c("recursive", "rolling", "fixed")==F){
    return("Unsupported forecasting sheme. Has to be either 'recursive', 'rolling' or 'fixed'.")
  }

  Data = data
  X = matrix(Data[,-1],nrow = dim(Data[,-1])[1])
  Y = matrix(Data[,1],nrow = dim(Data[,-1])[1])
  w_size = w_size
  n_windows = nrow(Data) - w_size
  predicts<-c()
  trues<-c()
  errors<-c()


  if(window=="recursive"){
    for (i in 1:n_windows) {
      ridge_cv <- cv.glmnet(x = X[1:(w_size + i - 1),], y = Y[1:(w_size + i - 1),],
                             type.measure = "mse",
                             nfold = 10,
                             alpha = 0)
      best_ridge_coef <- as.numeric(coef(ridge_cv, s = ridge_cv$lambda.1se))[-1]

      alasso.mod <- glmnet(x = X[1:(w_size + i - 1),], y = Y[1:(w_size + i - 1),],
                            alpha = 1,
                            penalty.factor = 1/abs(best_ridge_coef))

      alasso_cv <- cv.glmnet(x = X[1:(w_size + i - 1),], y = Y[1:(w_size + i - 1),],
                              type.measure = "mse",
                              nfold = 10,
                              alpha = 1,
                              penalty.factor = 1/abs(best_ridge_coef),
                              keep = TRUE)

      alasso_predict<-predict(alasso.mod, s=alasso_cv$lambda.min, newx=matrix(X[(w_size+i),], ncol=ncol(X)))
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
                            alpha = 0)
      best_ridge_coef <- as.numeric(coef(ridge_cv, s = ridge_cv$lambda.1se))[-1]

      alasso.mod <- glmnet(x = X[i:(w_size + i - 1),], y = Y[i:(w_size + i - 1),],
                           alpha = 1,
                           penalty.factor = 1/abs(best_ridge_coef))

      alasso_cv <- cv.glmnet(x = X[i:(w_size + i - 1),], y = Y[i:(w_size + i - 1),],
                             type.measure = "mse",
                             nfold = 10,
                             alpha = 1,
                             penalty.factor = 1/abs(best_ridge_coef),
                             keep = TRUE)

      alasso_predict<-predict(alasso.mod, s=alasso_cv$lambda.min, newx=matrix(X[(w_size+i),], ncol=ncol(X)))
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
                            alpha = 0)
      best_ridge_coef <- as.numeric(coef(ridge_cv, s = ridge_cv$lambda.1se))[-1]

      alasso.mod <- glmnet(x = X[1:w_size,], y = Y[1:w_size,],
                           alpha = 1,
                           penalty.factor = 1/abs(best_ridge_coef))

      alasso_cv <- cv.glmnet(x = X[1:w_size,], y = Y[1:w_size,],
                             type.measure = "mse",
                             nfold = 10,
                             alpha = 1,
                             penalty.factor = 1/abs(best_ridge_coef),
                             keep = TRUE)

      for(i in 1:n_windows){
        alasso_predict<-predict(alasso.mod, s=alasso_cv$lambda.min, newx=matrix(X[(w_size+i),], ncol=ncol(X)))
        y_real<-Y[w_size+i,]
        e<-y_real-alasso_predict

        predicts[i]<- alasso_predict
        trues[i]<- y_real
        errors[i]<- e
      }
    }
    mse<-mean(na.omit(errors)^2)

    forecasts<-data.frame(predicts, trues)
    colnames(forecasts) <-c('Forecasts','Realized')
    forecasts$Errors<-errors
    forecasts$Success <- ifelse(sign(forecasts$Forecasts) == sign(forecasts$Realized), 1, 0)
    success_ratio = sum(forecasts$Success)/nrow(forecasts)
    forecasts$Forecasted_Direction <- sign(forecasts$Forecasts)
    forecasts$Forecasted_Direction <- ifelse(sign(forecasts$Forecasted_Direction) == 1, 1, 0)
    forecasts$True_Direction <- sign(forecasts$Realized)
    forecasts$True_Direction <- ifelse(sign(forecasts$True_Direction) == 1, 1, 0)


    results<-list(Forecasts=forecasts, MSE=mse, SRatio=success_ratio)

    return(results)

}






















maeforecast.aralasso<-function(data=NULL, w_size=NULL, window="recursive"){
  if(is.null(data)|is.null(w_size)){
    return("Have to provide values for data and w_size")
  }
  if(window %in% c("recursive", "rolling", "fixed")==F){
    return("Unsupported forecasting sheme. Has to be either 'recursive', 'rolling' or 'fixed'.")
  }

  Data = data
  X = matrix(Data[,-1],nrow = dim(Data[,-1])[1])
  Y = matrix(Data[,1],nrow = dim(Data[,-1])[1])
  w_size = w_size
  n_windows = nrow(Data) - w_size
  predicts<-c()
  trues<-c()
  errors<-c()

  farima<-function(model=NULL, test=NULL){
    suppressMessages(require(stats))
    suppressMessages(require(forecast))
    suppressMessages(require(zoo))
    if(class(model)[1]!="ARIMA"){
      print("Error: the function only works with model fitted by Arima")
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
                            alpha = 0)
      best_ridge_coef <- as.numeric(coef(ridge_cv, s = ridge_cv$lambda.1se))[-1]

      alasso.mod <- glmnet(x = X[1:(w_size + i - 1),], y = Y[1:(w_size + i - 1),],
                           alpha = 1,
                           penalty.factor = 1/abs(best_ridge_coef))

      alasso_cv <- cv.glmnet(x = X[1:(w_size + i - 1),], y = Y[1:(w_size + i - 1),],
                             type.measure = "mse",
                             nfold = 10,
                             alpha = 1,
                             penalty.factor = 1/abs(best_ridge_coef),
                             keep = TRUE)


      best_alasso_coef <- coef(alasso_cv, s = alasso_cv$lambda.min)
      best_alasso_coef = as.numeric(best_alasso_coef) [-1]
      nonzero_index<-which(best_alasso_coef!=0)

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
                            alpha = 0)
      best_ridge_coef <- as.numeric(coef(ridge_cv, s = ridge_cv$lambda.1se))[-1]

      alasso.mod <- glmnet(x = X[i:(w_size + i - 1),], y = Y[i:(w_size + i - 1),],
                           alpha = 1,
                           penalty.factor = 1/abs(best_ridge_coef))

      alasso_cv <- cv.glmnet(x = X[i:(w_size + i - 1),], y = Y[i:(w_size + i - 1),],
                             type.measure = "mse",
                             nfold = 10,
                             alpha = 1,
                             penalty.factor = 1/abs(best_ridge_coef),
                             keep = TRUE)
      best_alasso_coef <- coef(alasso_cv, s = alasso_cv$lambda.min)
      best_alasso_coef = as.numeric(best_alasso_coef) [-1]
      nonzero_index<-which(best_alasso_coef!=0)

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
                          alpha = 0)
    best_ridge_coef <- as.numeric(coef(ridge_cv, s = ridge_cv$lambda.1se))[-1]

    alasso.mod <- glmnet(x = X[1:w_size,], y = Y[1:w_size,],
                         alpha = 1,
                         penalty.factor = 1/abs(best_ridge_coef))

    alasso_cv <- cv.glmnet(x = X[1:w_size,], y = Y[1:w_size,],
                           type.measure = "mse",
                           nfold = 10,
                           alpha = 1,
                           penalty.factor = 1/abs(best_ridge_coef),
                           keep = TRUE)

    best_alasso_coef <- coef(alasso_cv, s = alasso_cv$lambda.min)
    best_alasso_coef = as.numeric(best_alasso_coef) [-1]
    nonzero_index<-which(best_alasso_coef!=0)

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
  mse<-mean(na.omit(errors)^2)

  forecasts<-data.frame(predicts, trues)
  colnames(forecasts) <-c('Forecasts','Realized')
  forecasts$Errors<-errors
  forecasts$Success <- ifelse(sign(forecasts$Forecasts) == sign(forecasts$Realized), 1, 0)
  success_ratio = sum(forecasts$Success)/nrow(forecasts)
  forecasts$Forecasted_Direction <- sign(forecasts$Forecasts)
  forecasts$Forecasted_Direction <- ifelse(sign(forecasts$Forecasted_Direction) == 1, 1, 0)
  forecasts$True_Direction <- sign(forecasts$Realized)
  forecasts$True_Direction <- ifelse(sign(forecasts$True_Direction) == 1, 1, 0)


  results<-list(Forecasts=forecasts, MSE=mse, SRatio=success_ratio)

  return(results)

}
















maeforecast.arnet<-function(data=NULL, w_size=NULL, window="recursive", pred=NULL){
  if(is.null(data)|is.null(w_size)){
    return("Have to provide values for data and w_size")
  }
  if(window %in% c("recursive", "rolling", "fixed")==F){
    return("Unsupported forecasting sheme. Has to be either 'recursive', 'rolling' or 'fixed'.")
  }
  if(is.null(pred)){
    pred=w_size
  }else{
    pred=pred
  }

  Data = data
  X = matrix(Data[,-1],nrow = dim(Data[,-1])[1])
  Y = matrix(Data[,1],nrow = dim(Data[,-1])[1])
  w_size = w_size
  n_windows = nrow(Data) - w_size
  predicts<-c()
  trues<-c()
  errors<-c()

  farima<-function(model=NULL, test=NULL){
    suppressMessages(require(stats))
    suppressMessages(require(forecast))
    suppressMessages(require(zoo))
    if(class(model)[1]!="ARIMA"){
      print("Error: the function only works with model fitted by Arima")
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

  suppressMessages(require(glm2))
  suppressMessages(require(msaenet))
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
                           pred=pred)

      aenet.fit <- aenet(Xsis$Xs, matrix(Y[1:(w_size + i - 1),1]),
                         alphas = seq(0.2, 0.8, 0.02), seed = 10)

      nonzero_index = which(!coef(aenet.fit) == 0)

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
                           pred=pred)

      aenet.fit <- aenet(Xsis$Xs, matrix(Y[i:(w_size + i - 1),1]),
                         alphas = seq(0.2, 0.8, 0.02), seed = 10)

      nonzero_index = which(!coef(aenet.fit) == 0)

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
                         pred=pred)

    aenet.fit <- aenet(Xsis$Xs, matrix(Y[1:w_size,1]),
                       alphas = seq(0.2, 0.8, 0.02), seed = 10)

    nonzero_index = which(!coef(aenet.fit) == 0)


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


  mse<-mean(na.omit(errors)^2)

  forecasts<-data.frame(predicts, trues)
  colnames(forecasts) <-c('Forecasts','Realized')
  forecasts$Errors<-errors
  forecasts$Success <- ifelse(sign(forecasts$Forecasts) == sign(forecasts$Realized), 1, 0)
  success_ratio = sum(forecasts$Success)/nrow(forecasts)
  forecasts$Forecasted_Direction <- sign(forecasts$Forecasts)
  forecasts$Forecasted_Direction <- ifelse(sign(forecasts$Forecasted_Direction) == 1, 1, 0)
  forecasts$True_Direction <- sign(forecasts$Realized)
  forecasts$True_Direction <- ifelse(sign(forecasts$True_Direction) == 1, 1, 0)


  results<-list(Forecasts=forecasts, MSE=mse, SRatio=success_ratio)

  return(results)
}














maeforecast.dfm<-function(data=NULL, w_size=NULL, window="recursive", factor.num=3){
  if(is.null(data)|is.null(w_size)){
    return("Have to provide values for data and w_size")
  }
  if(window %in% c("recursive", "rolling", "fixed")==F){
    return("Unsupported forecasting sheme. Has to be either 'recursive', 'rolling' or 'fixed'.")
  }

  Data = data
  X = matrix(Data[,-1],nrow = dim(Data[,-1])[1])
  Y = matrix(Data[,1],nrow = dim(Data[,-1])[1])
  w_size = w_size
  n_windows = nrow(Data) - w_size
  predicts<-c()
  trues<-c()
  errors<-c()

  suppressMessages(require(nowcasting))
  suppressMessages(require(forecast))

  dfm<-function (X, r, p, q, max_iter = 100, threshold = 1e-04, rQ, rC){
    if (missing(rQ)) {
      rQ <- ""
    }
    if (missing(rC)) {
      rC <- ""
    }
    if (missing(q)) {
      q <- r
    }
    if (q > r) {
      stop("r must be larger than q.")
    }
    T <- dim(X)[1]
    N <- dim(X)[2]
    x <- apply(X, 2, function(z) {
      (z - mean(z, na.rm = TRUE))/sd(z, na.rm = TRUE)
    })
    Mx <- apply(X, 2, mean, na.rm = TRUE)
    Wx <- apply(X, 2, sd, na.rm = TRUE)
    W <- !is.na(x)
    A <- rbind(matrix(0, nrow = r, ncol = r * p), diag(1, nrow = r *
                                                         (p - 1), ncol = r * p))
    Q <- matrix(0, nrow = p * r, ncol = p * r)
    Q[1:r, 1:r] <- diag(1, r)
    eigen.decomp <- eigen(cov(x, use = "complete.obs"))
    v <- eigen.decomp$vectors[, 1:r]
    d <- eigen.decomp$values[1:r]
    chi <- x %*% v %*% t(v)
    d <- diag(1, r)
    F <- x %*% v
    F_pc <- F
    F <- na.omit(F)
    if (p > 0) {
      if (rQ == "identity") {
        fit <- VAR(F, p)
        A[1:r, 1:(r * p)] <- t(fit$A)
        Q[1:r, 1:r] <- diag(1, r)
      }
      else {
        fit <- VAR(F, p)
        A[1:r, 1:(r * p)] <- t(fit$A)
        H <- cov(fit$res)
        if (r > q) {
          q.decomp <- eigen(H)
          P <- q.decomp$vectors[, 1:q, drop = FALSE]
          M <- q.decomp$values[1:q]
          if (q == 1) {
            P <- P * P[1, ]
            Q[1:r, 1:r] <- P %*% t(P) * M
          }
          else {
            P <- P %*% diag(sign(P[1, ]))
            Q[1:r, 1:r] <- P %*% diag(M) %*% t(P)
          }
        }
        else {
          Q[1:r, 1:r] <- H
        }
      }
    }
    R <- diag(diag(cov(x - chi, use = "complete.obs")))
    Z <- fit$X
    initx <- Z[1, ]
    initV <- matrix(ginv(kronecker(A, A)) %*% as.numeric(Q),
                    ncol = r * p, nrow = r * p)
    C <- cbind(v, matrix(0, nrow = N, ncol = r * (p - 1)))
    previous_loglik <- -.Machine$double.xmax
    loglik <- 0
    num_iter <- 0
    LL <- c()
    converged <- 0
    kf_res <- K_filter(initx, initV, x, A, C, R, Q)
    ks_res <- K_smoother(A, kf_res$xitt, kf_res$xittm, kf_res$Ptt,
                         kf_res$Pttm, C, R, W)
    xsmooth <- ks_res$xitT
    Vsmooth <- ks_res$PtT
    Wsmooth <- ks_res$PtTm
    F_kal <- t(xsmooth[1:r, , drop = FALSE])
    if (rC == "upper" & (r > 1)) {
      dimC <- dim(C[, 1:r])
      rK <- rep(0, (r - 1) * r/2)
      irC <- which(matrix(upper.tri(C[, 1:r]) + 0) == 1)
      rH <- matrix(0, nrow = length(rK), ncol = prod(dimC))
      for (i in 1:length(rK)) {
        rH[i, irC[i]] <- 1
      }
    }
    while ((num_iter < max_iter) & !converged) {
      em_res <- Estep(t(x), A, C, Q, R, initx, initV, W)
      beta <- em_res$beta_t
      gamma <- em_res$gamma_t
      delta <- em_res$delta_t
      gamma1 <- em_res$gamma1_t
      gamma2 <- em_res$gamma2_t
      P1sum <- em_res$V1 + em_res$x1 %*% t(em_res$x1)
      x1sum <- em_res$x1
      loglik <- em_res$loglik_t
      num_iter <- num_iter + 1
      if (rC == "upper" & (r > 1)) {
        fp <- matrix(delta[, 1:r] %*% ginv(gamma[1:r, 1:r]))
        kronCR <- kronecker(ginv(gamma[1:r, 1:r]), R)
        sp <- kronCR %*% t(rH) %*% ginv(rH %*% kronCR %*%
                                          t(rH)) %*% (rK - rH %*% fp)
        C[, 1:r] <- matrix(fp + sp, nrow = dimC[1], ncol = dimC[2])
      }
      else {
        C[, 1:r] <- delta[, 1:r] %*% ginv(gamma[1:r, 1:r])
      }
      if (p > 0) {
        A_update <- beta[1:r, 1:(r * p), drop = FALSE] %*%
          solve(gamma1[1:(r * p), 1:(r * p)])
        A[1:r, 1:(r * p)] <- A_update
        if (rQ != "identity") {
          H <- (gamma2[1:r, 1:r] - A_update %*% t(beta[1:r,
                                                       1:(r * p), drop = FALSE]))/(T - 1)
          if (r > q) {
            h.decomp <- svd(H)
            P <- h.decomp$v[, 1:q, drop = FALSE]
            M <- h.decomp$d[1:q]
            if (q == 1) {
              P <- P * P[1, ]
              Q[1:r, 1:r] <- P %*% t(P) * M
            }
            else {
              P <- P %*% diag(sign(P[1, ]))
              Q[1:r, 1:r] <- P %*% diag(M) %*% t(P)
            }
          }
          else {
            Q[1:r, 1:r] <- H
          }
        }
      }
      xx <- as.matrix(na.omit(x))
      R <- (t(xx) %*% xx - C %*% t(delta))/T
      RR <- diag(R)
      RR[RR < 1e-07] <- 1e-07
      R <- diag(RR)
      R <- diag(diag(R))
      LL <- c(LL, loglik)
      initx <- x1sum
      initV <- P1sum - initx %*% t(initx)
      converged <- em_converged(loglik, previous_loglik, threshold = threshold)
      previous_loglik <- loglik
      if (num_iter < 25) {
        converged <- FALSE
      }
    }
    if (converged == TRUE) {
      cat("Converged after", num_iter, "iterations.\n")
    }
    else {
      cat("Maximum number of iterations reached.\n")
    }
    kf <- K_filter(initx, initV, x, A, C, R, Q)
    ks <- K_smoother(A, kf$xitt, kf$xittm, kf$Ptt, kf$Pttm, C,
                     R, W)
    xsmooth <- ks$xitT
    chi <- t(xsmooth) %*% t(C) %*% diag(Wx) + kronecker(matrix(1,
                                                               T, 1), t(Mx))
    F_hat <- t(xsmooth[1:r, , drop = FALSE])
    final_object <- list(pca = F_pc, qml = F_hat, twostep = F_kal,
                         A = A[1:r, ], C = C[, 1:r], Q = Q[1:q, 1:q], R = R, p = p,
                         data = x)
    class(final_object) <- c("dfm", "list")
    return(final_object)
  }

  if(window=="recursive"){
    for(i in 1:n_windows){
      suppressMessages(fac.mod<-dfm(X[1:(w_size + i),], r=factor.num, p=1, q=factor.num))
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
      suppressMessages(fac.mod<-dfm(X[i:(w_size + i),], r=factor.num, p=1, q=factor.num))
      factors<-fac.mod$twostep[i:(w_size + i - 1),]
      newfactors<-matrix(fac.mod$twostep[w_size+i,], nrow=1)

      AR_dfm<-Arima(Y[i:(w_size+i-1), ], order=c(1,0,0), xreg=factors)
      AR_dfm_predict<-as.numeric(predict(AR_dfm, newxreg=newfactors, n.ahead=1)$pred)

      y_real<-Y[w_size+i,]
      e<-y_real-AR_dfm_predict

      trues[i] <- y_real
      predicts[i] <- AR_dfm_predict
      errors[i] <- e
    }
  }else{
    suppressMessages(fac.mod<-dfm(X, r=factor.num, p=1, q=factor.num))
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
  mse<-mean(na.omit(errors)^2)

  forecasts<-data.frame(predicts, trues)
  colnames(forecasts) <-c('Forecasts','Realized')
  forecasts$Errors<-errors
  forecasts$Success <- ifelse(sign(forecasts$Forecasts) == sign(forecasts$Realized), 1, 0)
  success_ratio = sum(forecasts$Success)/nrow(forecasts)
  forecasts$Forecasted_Direction <- sign(forecasts$Forecasts)
  forecasts$Forecasted_Direction <- ifelse(sign(forecasts$Forecasted_Direction) == 1, 1, 0)
  forecasts$True_Direction <- sign(forecasts$Realized)
  forecasts$True_Direction <- ifelse(sign(forecasts$True_Direction) == 1, 1, 0)


  results<-list(Forecasts=forecasts, MSE=mse, SRatio=success_ratio)

  return(results)
}













maeforecast.ar<-function(data=NULL, w_size=NULL, window="recursive"){
  if(is.null(data)|is.null(w_size)){
    return("Have to provide values for data and w_size")
  }
  if(window %in% c("recursive", "rolling", "fixed")==F){
    return("Unsupported forecasting sheme. Has to be either 'recursive', 'rolling' or 'fixed'.")
  }

  Data = ts(data[,1], frequency=12)
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
      print("Error: the function only works with model fitted by Arima")
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
  mse<-mean(na.omit(errors)^2)

  forecasts<-data.frame(predicts, trues)
  colnames(forecasts) <-c('Forecasts','Realized')
  forecasts$Errors<-errors
  forecasts$Success <- ifelse(sign(forecasts$Forecasts) == sign(forecasts$Realized), 1, 0)
  success_ratio = sum(forecasts$Success)/nrow(forecasts)
  forecasts$Forecasted_Direction <- sign(forecasts$Forecasts)
  forecasts$Forecasted_Direction <- ifelse(sign(forecasts$Forecasted_Direction) == 1, 1, 0)
  forecasts$True_Direction <- sign(forecasts$Realized)
  forecasts$True_Direction <- ifelse(sign(forecasts$True_Direction) == 1, 1, 0)


  results<-list(Forecasts=forecasts, MSE=mse, SRatio=success_ratio)

  return(results)
}








maeforecast<-function(data=NULL, model="ar", w_size=NULL, window="recursive", pred=NULL, factor.num=NULL){
  if(model %in% c("lasso", "arlasso", "ridge", "alasso", "aralasso", "arnet", "dfm", "ar")==F){
    return("Unsupported model type.")
  }
  if(model=="ar"){
    results<-maeforecast.ar(data=data, w_size=w_size, window=window)
  }else if(model=="lasso"){
    results<-maeforecast.lasso(data=data, w_size=w_size, window=window)
  }else if(model=="arlasso"){
    results<-maeforecast.arlasso(data=data, w_size=w_size, window=window)
  }else if(model=="ridge"){
    results<-maeforecast.ridge(data=data, w_size=w_size, window=window)
  }else if(model=="alasso"){
    results<-maeforecast.alasso(data=data, w_size=w_size, window=window)
  }else if (model=="aralasso"){
    results<-maeforecast.aralasso(data=data, w_size=w_size, window=window)
  }else if(model=="arnet"){
    results<-maeforecast.arnet(data=data, w_size=w_size, window=window, pred=pred)
  }else{
    results<-maeforecast.dfm(data=data, w_size=w_size, window=window, factor.num=factor.num)
  }
  return(results)
}
