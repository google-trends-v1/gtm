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

  suppressMessages(require(dynfactoR))
  suppressMessages(require(nowcasting))
  suppressMessages(require(forecast))

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
