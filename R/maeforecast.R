maeforecast.lasso<-function(data=NULL, w_size=NULL, window="recursive", h=0, standardize=TRUE, lambda=NULL, y.index=1, t.select=NULL, t.update=FALSE){
  if(is.null(data)|is.null(w_size)){
    stop("Have to provide values for data and w_size.")
  }
  if(window %in% c("recursive", "rolling", "fixed")==F){
    stop("Unsupported forecasting sheme. Has to be either 'recursive', 'rolling' or 'fixed'.")
  }

  Data = as.matrix(data)

  if(!is.null(colnames(data))){
    varnames<-colnames(data)[-y.index]
  }else{
    varnames<-seq(from=1, to=dim(data)[2], by=1)[-y.index]
  }

  X = matrix(Data[1:(dim(Data)[1]-h),-y.index],nrow = (dim(Data)[1]-h))
  Y = matrix(Data[(1+h):dim(Data)[1],y.index],nrow = (dim(Data)[1]-h))

  w_size = w_size
  n_windows = dim(Y)[1] - w_size
  predicts<-c()
  trues<-c()
  errors<-c()
  covariates<-list()

  if(t.update==F){
    if(class(t.select)=="numeric"){
      ranks<-t.rank(Data[1:(w_size), ], y.index=y.index, h=h)
      X.index<-ranks[1:t.select,1]
      varnames.selected<-varnames[X.index]
    }else{
      X.index<-1:ncol(X)
      varnames.selected<-varnames
  }

  suppressMessages(require(glmnet))

  if(window=="recursive"){
    for(i in 1:n_windows){
      if(t.update){
        if(class(t.select)=="numeric"){
          ranks<-t.rank(Data[1:(w_size+i-1), ], y.index=y.index, h=h)
          X.index<-ranks[1:t.select,1]
          varnames.selected<-varnames[X.index]
        }else{
          X.index<-1:ncol(X)
          varnames.selected<-varnames
        }
      }

      lasso.mod<-glmnet(x = X[1:(w_size + i - 1), X.index], y = Y[1:(w_size + i - 1),],
                        alpha=1,
                        standardize=standardize)

      lasso_cv<-cv.glmnet(x = X[1:(w_size + i - 1), X.index], y = Y[1:(w_size + i - 1),],
                          type.measure="mse",
                          nfold=10,
                          alpha=1,
                          standardize=standardize,
                          lambda=lambda)


      best_lasso_coef <- coef(lasso_cv, s=lasso_cv$lambda.1se)
      best_lasso_coef = as.numeric(best_lasso_coef) [-1]
      nonzero_index<-which(best_lasso_coef!=0)
      covariates[i]<-list(varnames.selected[nonzero_index])

      lasso_predict<-predict(lasso.mod, s=lasso_cv$lambda.1se, newx=matrix(X[(w_size+i), X.index], ncol=length(X.index)))
      y_real<-Y[w_size+i,]
      e<-y_real-lasso_predict

      predicts[i]<-lasso_predict
      trues[i]<-y_real
      errors[i]<-e
    }
  }else if(window=="rolling"){
    for(i in 1:n_windows){

      if(t.update){
        if(class(t.select)=="numeric"){
          ranks<-t.rank(Data[i:(w_size+i-1), ], y.index=y.index, h=h)
          X.index<-ranks[1:t.select,1]
          varnames.selected<-varnames[X.index]
        }else{
          X.index<-1:ncol(X)
          varnames.selected=varnames
        }
      }


      lasso.mod<-glmnet(x = X[i:(w_size + i - 1), X.index], y = Y[i:(w_size + i - 1),],
                        alpha=1, standardize=standardize)

      lasso_cv<-cv.glmnet(x = X[i:(w_size + i - 1), X.index], y = Y[i:(w_size + i - 1),],
                          type.measure="mse",
                          nfold=10,
                          alpha=1,
                          standardize=standardize,
                          lambda=lambda)

      best_lasso_coef <- coef(lasso_cv, s = lasso_cv$lambda.1se)
      best_lasso_coef = as.numeric(best_lasso_coef) [-1]
      nonzero_index<-which(best_lasso_coef!=0)
      covariates[i]<-list(varnames.selected[nonzero_index])

      lasso_predict<-predict(lasso.mod, s=lasso_cv$lambda.1se, newx=matrix(X[(w_size+i), X.index], ncol=length(X.index)))
      y_real<-Y[w_size+i,]
      e<-y_real-lasso_predict

      predicts[i]<-lasso_predict
      trues[i]<-y_real
      errors[i]<-e
    }
  }else{

    if(t.update){
      if(class(t.select)=="numeric"){
        ranks<-t.rank(Data[1:w_size, ], y.index=y.index, h=h)
        X.index<-ranks[1:t.select,1]
        varnames.selected<-varnames[X.index]
      }else{
        X.index<-1:ncol(X)
        varnames.selected<-varnames
      }
    }


    lasso.mod<-glmnet(x = X[1:w_size, X.index], y = Y[1:w_size,],
                      alpha=1, standardize=standardize)

    lasso_cv<-cv.glmnet(x = X[1:w_size, X.index], y = Y[1:w_size,],
                        type.measure="mse",
                        nfold=10,
                        alpha=1,
                        standardize=standardize,
                        lambda=lambda)

    best_lasso_coef <- coef(lasso_cv, s = lasso_cv$lambda.1se)
    best_lasso_coef = as.numeric(best_lasso_coef) [-1]
    nonzero_index<-which(best_lasso_coef!=0)
    covariates<-list(varnames.selected[nonzero_index])

    for(i in 1:n_windows){
      lasso_predict<-predict(lasso.mod, s=lasso_cv$lambda.1se, newx=matrix(X[(w_size+i), X.index], ncol=length(X.index)))
      y_real<-Y[w_size+i,]
      e<-y_real-lasso_predict

      predicts[i]<-lasso_predict
      trues[i]<-y_real
      errors[i]<-e
    }
  }
  results<-Metrics(pred=predicts, true=trues, h=h)
  results$Data<-Data
  results$Model<-list(Model="Lasso", Window=window, Size=w_size, Horizon=h, Preselection=t.select, Update=t.update, Index=y.index, Standardize=standardize, Lambda=lambda)
  results$Variables<-covariates
  class(results)<-"Maeforecast"
  return(results)
  }
}













maeforecast.postlasso<-function(data=NULL, w_size=NULL, window="recursive", h=0, standardize=TRUE, lambda=NULL, y.index=1, t.select=NULL, t.update=FALSE){
  if(is.null(data)|is.null(w_size)){
    stop("Have to provide values for data and w_size.")
  }
  if(window %in% c("recursive", "rolling", "fixed")==F){
    stop("Unsupported forecasting sheme. Has to be either 'recursive', 'rolling' or 'fixed'.")
  }

  Data = as.matrix(data)
  if(!is.null(colnames(data))){
    varnames<-colnames(data)[-y.index]
  }else{
    varnames<-seq(from=1, to=dim(data)[2], by=1)[-y.index]
  }

  X = matrix(Data[1:(dim(Data)[1]-h),-y.index],nrow = (dim(Data)[1]-h))
  Y = matrix(Data[(1+h):dim(Data)[1],y.index],nrow = (dim(Data)[1]-h))

  w_size = w_size
  n_windows = dim(Y)[1] - w_size
  predicts<-c()
  trues<-c()
  errors<-c()
  covariates<-list()

  if(t.update==F){
    if(class(t.select)=="numeric"){
      ranks<-t.rank(Data[1:(w_size), ], y.index=y.index, h=h)
      X.index<-ranks[1:t.select,1]
      varnames.selected<-varnames[X.index]
    }else{
      X.index<-1:ncol(X)
      varnames.selected<-varnames
    }
  }

  suppressMessages(require(glmnet))
  suppressMessages(require(forecast))

  if(window=="recursive"){
    for(i in 1:n_windows){

      if(t.update){
        if(class(t.select)=="numeric"){
          ranks<-t.rank(Data[1:(w_size+i-1), ], y.index=y.index, h=h)
          X.index<-ranks[1:t.select,1]
          varnames.selected<-varnames[X.index]
          }else{
            X.index<-1:ncol(X)
            varnames.selected<-varnames
          }
        }

      lasso.mod<-glmnet(x = X[1:(w_size + i - 1), X.index], y = Y[1:(w_size + i - 1),],
                        alpha=1, standardize=standardize)

      lasso_cv<-cv.glmnet(x = X[1:(w_size + i - 1), X.index], y = Y[1:(w_size + i - 1),],
                          type.measure="mse",
                          nfold=10,
                          alpha=1,
                          standardize=standardize,
                          lambda=lambda)


      best_lasso_coef <- coef(lasso_cv, s = lasso_cv$lambda.1se)
      best_lasso_coef = as.numeric(best_lasso_coef) [-1]
      nonzero_index<-which(best_lasso_coef!=0)
      covariates[i]<-list(varnames.selected[nonzero_index])

      if (length(nonzero_index) !=0){
        nonzero_index <-as.vector(nonzero_index)
        xregs <- X[1:(w_size + i - 1), X.index][, nonzero_index]
        newxregs <-matrix(matrix(X[(w_size + i), X.index], ncol=length(X.index), nrow=1)[, nonzero_index], nrow=1)

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



      if(t.update){
        if(class(t.select)=="numeric"){

          ranks<-t.rank(Data[i:(w_size+i-1), ], y.index=y.index, h=h)

          X.index<-ranks[1:t.select,1]

          varnames.selected<-varnames[X.index]

          }else{

            X.index<-1:ncol(X)

            varnames.selected<-varnames

          }

        }

      lasso.mod<-glmnet(x = X[i:(w_size + i - 1), X.index], y = Y[i:(w_size + i - 1),],
                        alpha=1, standardize=standardize)

      lasso_cv<-cv.glmnet(x = X[i:(w_size + i - 1), X.index], y = Y[i:(w_size + i - 1),],
                          type.measure="mse",
                          nfold=10,
                          alpha=1,
                          standardize=standardize,
                          lambda=lambda)

      best_lasso_coef <- coef(lasso_cv, s = lasso_cv$lambda.1se)
      best_lasso_coef = as.numeric(best_lasso_coef) [-1]
      nonzero_index<-which(best_lasso_coef!=0)
      covariates[i]<-list(varnames.selected[nonzero_index])

      if (length(nonzero_index) !=0){
        nonzero_index <-as.vector(nonzero_index)
        xregs <- X[i:(w_size + i - 1), X.index][, nonzero_index]
        newxregs <-matrix(matrix(X[(w_size + i), X.index], ncol=length(X.index), nrow=1)[, nonzero_index], nrow=1)

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

    if(t.update){
      if(class(t.select)=="numeric"){

        ranks<-t.rank(Data[1:w_size, ], y.index=y.index, h=h)

        X.index<-ranks[1:t.select,1]

        varnames.selected<-varnames[X.index]

         }else{

           X.index<-1:ncol(X)

           varnames.selected<-varnames

         }
      }

    lasso.mod<-glmnet(x = X[1:w_size, X.index], y = Y[1:w_size,],
                      alpha=1, standardize=standardize)

    lasso_cv<-cv.glmnet(x = X[1:w_size, X.index], y = Y[1:w_size,],
                        type.measure="mse",
                        nfold=10,
                        alpha=1,
                        standardize=standardize,
                        lambda=lambda)

    best_lasso_coef <- coef(lasso_cv, s = lasso_cv$lambda.1se)
    best_lasso_coef = as.numeric(best_lasso_coef) [-1]
    nonzero_index<-which(best_lasso_coef!=0)
    covariates<-list(varnames.selected[nonzero_index])

    if (length(nonzero_index) !=0){
      nonzero_index <-as.vector(nonzero_index)
      Data_uni = ts(Y, frequency=12)
      xregs <- X[1:w_size, X.index][, nonzero_index]
      AR_lasso <-Arima(Data_uni[1:w_size], order=c(1,0,0), xreg =xregs)
      for(i in 1:n_windows){
        newxregs <-matrix(matrix(X[(w_size + i), X.index], ncol=length(X.index), nrow=1)[, nonzero_index], nrow=1)

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
  results<-Metrics(pred=predicts, true=trues, h=h)
  results$Data<-Data
  results$Model<-list(Model="Post Lasso", Window=window, Size=w_size, Horizon=h, Preselection=t.select, Update=t.update, Index=y.index, Standardize=standardize, Lambda=lambda)
  results$Variables<-covariates
  class(results)<-"Maeforecast"
  return(results)
}








maeforecast.ridge<-function(data=NULL, w_size=NULL, window="recursive", h=0, standardize=TRUE, lambda=NULL, y.index=1, t.select=NULL, t.update=FALSE){
  if(is.null(data)|is.null(w_size)){
    stop("Have to provide values for data and w_size.")
  }
  if(window %in% c("recursive", "rolling", "fixed")==F){
    stop("Unsupported forecasting sheme. Has to be either 'recursive', 'rolling' or 'fixed'.")
  }

  Data = as.matrix(data)
  if(!is.null(colnames(data))){
    varnames<-colnames(data)[-y.index]
  }else{
    varnames<-seq(from=1, to=dim(data)[2], by=1)[-y.index]
  }

  X = matrix(Data[1:(dim(Data)[1]-h),-y.index],nrow = (dim(Data)[1]-h))
  Y = matrix(Data[(1+h):dim(Data)[1],y.index],nrow = (dim(Data)[1]-h))

  w_size = w_size
  n_windows = dim(Y)[1] - w_size
  predicts<-c()
  trues<-c()
  errors<-c()

  if(t.update==F){
    if(class(t.select)=="numeric"){
      ranks<-t.rank(Data[1:(w_size), ], y.index=y.index, h=h)
      X.index<-ranks[1:t.select,1]
      varnames.selected<-varnames[X.index]
    }else{
      X.index<-1:ncol(X)
      varnames.selected<-varnames
    }
  }

  suppressMessages(require(glmnet))

  if(window=="recursive"){
    for(i in 1:n_windows){

      if(t.update){
        if(class(t.select)=="numeric"){

          ranks<-t.rank(Data[1:(w_size+i-1), ], y.index=y.index, h=h)

          X.index<-ranks[1:t.select,1]

          varnames.selected<-varnames[X.index]

          }else{

            X.index<-1:ncol(X)

            varnames.selected<-varnames
          }

        }

      ridge.mod<-glmnet(x = X[1:(w_size + i - 1), X.index], y = Y[1:(w_size + i - 1),],
                        alpha=0, standardize=standardize)

      ridge_cv<-cv.glmnet(x = X[1:(w_size + i - 1), X.index], y = Y[1:(w_size + i - 1),],
                          type.measure="mse",
                          nfold=10,
                          alpha=0,
                          standardize=standardize,
                          lambda=lambda)

      ridge_predict<-predict(ridge.mod, s=ridge_cv$lambda.1se, newx=matrix(X[(w_size+i), X.index], ncol=length(X.index)))
      y_real<-Y[w_size+i,]
      e<-y_real-ridge_predict

      predicts[i] <- ridge_predict
      trues[i] <- y_real
      errors[i] <- e
    }
  }else if (window=="rolling"){
    for(i in 1:n_windows){

      if(t.update){
        if(class(t.select)=="numeric"){
          ranks<-t.rank(Data[i:(w_size+i-1), ], y.index=y.index, h=h)

          X.index<-ranks[1:t.select,1]

          varnames.selected<-varnames[X.index]
      }else{

        X.index<-1:ncol(X)

        varnames.selected<-varnames
      }

        }

      ridge.mod<-glmnet(x = X[i:(w_size + i - 1), X.index], y = Y[i:(w_size + i - 1),],
                        alpha=0, standardize=standardize)

      ridge_cv<-cv.glmnet(x = X[i:(w_size + i - 1), X.index], y = Y[i:(w_size + i - 1),],
                          type.measure="mse",
                          nfold=10,
                          alpha=0,
                          standardize=standardize,
                          lambda=lambda)

      ridge_predict<-predict(ridge.mod, s=ridge_cv$lambda.1se, newx=matrix(X[(w_size+i), X.index], ncol=length(X.index)))
      y_real<-Y[w_size+i,]
      e<-y_real-ridge_predict

      predicts[i] <- ridge_predict
      trues[i] <- y_real
      errors[i] <- e
    }
  }else{

    if(t.update){
      if(class(t.select)=="numeric"){

        ranks<-t.rank(Data[1:w_size, ], y.index=y.index, h=h)

        X.index<-ranks[1:t.select,1]

        varnames.selected<-varnames[X.index]

        }else{

          X.index<-1:ncol(X)

          varnames.selected<-varnames
        }
      }

    ridge.mod<-glmnet(x = X[1:w_size, X.index], y = Y[1:w_size,],
                      alpha=0, standardize=standardize)

    ridge_cv<-cv.glmnet(x = X[1:w_size, X.index], y = Y[1:w_size,],
                        type.measure="mse",
                        nfold=10,
                        alpha=0,
                        standardize=standardize,
                        lambda=lambda)

    for(i in 1:n_windows){
      ridge_predict<-predict(ridge.mod, s=ridge_cv$lambda.1se, newx=matrix(X[(w_size+i), X.index], ncol=length(X.index)))
      y_real<-Y[w_size+i,]
      e<-y_real-ridge_predict

      predicts[i] <- ridge_predict
      trues[i] <- y_real
      errors[i] <- e
    }
  }
  results<-Metrics(pred=predicts, true=trues, h=h)
  results$Data<-Data
  results$Model<-list(Model="Ridge", Window=window, Size=w_size, Horizon=h, Preselection=t.select, Update=t.update, Index=y.index, Standardize=standardize, Lambda=lambda)
  class(results)<-"Maeforecast"
  return(results)
}









maeforecast.alasso<-function(data=NULL, w_size=NULL, window="recursive", h=0, standardize=TRUE, lambda.ridge=NULL, lambda.lasso=NULL, y.index=1, t.select=NULL, t.update=FALSE){
  if(is.null(data)|is.null(w_size)){
    stop("Have to provide values for data and w_size.")
  }
  if(window %in% c("recursive", "rolling", "fixed")==F){
    stop("Unsupported forecasting sheme. Has to be either 'recursive', 'rolling' or 'fixed'.")
  }

  Data = as.matrix(data)
  if(!is.null(colnames(data))){
    varnames<-colnames(data)[-y.index]
  }else{
    varnames<-seq(from=1, to=dim(data)[2], by=1)[-y.index]
  }

  X = matrix(Data[1:(dim(Data)[1]-h),-y.index],nrow = (dim(Data)[1]-h))
  Y = matrix(Data[(1+h):dim(Data)[1],y.index],nrow = (dim(Data)[1]-h))

  w_size = w_size
  n_windows = dim(Y)[1] - w_size
  predicts<-c()
  trues<-c()
  errors<-c()
  covariates<-list()

  if(t.update==F){
    if(class(t.select)=="numeric"){
      ranks<-t.rank(Data[1:(w_size), ], y.index=y.index, h=h)
      X.index<-ranks[1:t.select,1]
      varnames.selected<-varnames[X.index]
    }else{
      X.index<-1:ncol(X)
      varnames.selected<-varnames
    }
  }

  suppressMessages(require(glmnet))

  if(window=="recursive"){
    for (i in 1:n_windows) {

      {
        if(class(t.select)=="numeric"){

           ranks<-t.rank(Data[1:(w_size+i-1), ], y.index=y.index, h=h)

           X.index<-ranks[1:t.select,1]

           varnames.selected<-varnames[X.index]

           }else{

             X.index<-1:ncol(X)

             varnames.selected<-varnames

           }
        }

      ridge_cv <- cv.glmnet(x = X[1:(w_size + i - 1), X.index], y = Y[1:(w_size + i - 1),],
                             type.measure = "mse",
                             nfold = 10,
                             alpha = 0,
                            standardize=standardize,
                            lambda=lambda.ridge)
      best_ridge_coef <- as.numeric(coef(ridge_cv, s = ridge_cv$lambda.1se))[-1]

      alasso.mod <- glmnet(x = X[1:(w_size + i - 1), X.index], y = Y[1:(w_size + i - 1),],
                            alpha = 1,
                            penalty.factor = 1/abs(best_ridge_coef),
                           standardize=standardize)

      alasso_cv <- cv.glmnet(x = X[1:(w_size + i - 1), X.index], y = Y[1:(w_size + i - 1),],
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
      covariates[i]<-list(varnames.selected[nonzero_index])

      alasso_predict<-predict(alasso.mod, s=alasso_cv$lambda.1se, newx=matrix(X[(w_size+i), X.index], ncol=length(X.index)))
      y_real<-Y[w_size+i,]
      e<-y_real-alasso_predict

      predicts[i]<- alasso_predict
      trues[i]<- y_real
      errors[i]<- e
    }
  }else if(window=="rolling"){
    for (i in 1:n_windows) {

      if(t.update){
        if(class(t.select)=="numeric"){

          ranks<-t.rank(Data[i:(w_size+i-1), ], y.index=y.index, h=h)

          X.index<-ranks[1:t.select,1]

          varnames.selected<-varnames[X.index]

          }else{

            X.index<-1:ncol(X)

            varnames.selected<-varnames

          }
        }

      ridge_cv <- cv.glmnet(x = X[i:(w_size + i - 1), X.index], y = Y[i:(w_size + i - 1),],
                            type.measure = "mse",
                            nfold = 10,
                            alpha = 0,
                            standardize=standardize,
                            lambda=lambda.ridge)
      best_ridge_coef <- as.numeric(coef(ridge_cv, s = ridge_cv$lambda.1se))[-1]

      alasso.mod <- glmnet(x = X[i:(w_size + i - 1), X.index], y = Y[i:(w_size + i - 1),],
                           alpha = 1,
                           penalty.factor = 1/abs(best_ridge_coef),
                           standardize=standardize)

      alasso_cv <- cv.glmnet(x = X[i:(w_size + i - 1), X.index], y = Y[i:(w_size + i - 1),],
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
      covariates[i]<-list(varnames.selected[nonzero_index])

      alasso_predict<-predict(alasso.mod, s=alasso_cv$lambda.1se, newx=matrix(X[(w_size+i), X.index], ncol=length(X.index)))
      y_real<-Y[w_size+i,]
      e<-y_real-alasso_predict

      predicts[i]<- alasso_predict
      trues[i]<- y_real
      errors[i]<- e
    }
  }else{

    if(t.update){
      if(class(t.select)=="numeric"){

        ranks<-t.rank(Data[1:w_size, ], y.index=y.index, h=h)

        X.index<-ranks[1:t.select,1]

        varnames.selected<-varnames[X.index]

         }else{

           X.index<-1:ncol(X)

           varnames.selected<-varnames

         }
      }

      ridge_cv <- cv.glmnet(x = X[1:w_size, X.index], y = Y[1:w_size,],
                            type.measure = "mse",
                            nfold = 10,
                            alpha = 0,
                            standardize=standardize,
                            lambda=lambda.ridge)
      best_ridge_coef <- as.numeric(coef(ridge_cv, s = ridge_cv$lambda.1se))[-1]

      alasso.mod <- glmnet(x = X[1:w_size, X.index], y = Y[1:w_size,],
                           alpha = 1,
                           penalty.factor = 1/abs(best_ridge_coef),
                           standardize=standardize)

      alasso_cv <- cv.glmnet(x = X[1:w_size, X.index], y = Y[1:w_size,],
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
      covariates<-list(varnames.selected[nonzero_index])

      for(i in 1:n_windows){
        alasso_predict<-predict(alasso.mod, s=alasso_cv$lambda.1se, newx=matrix(X[(w_size+i), X.index], ncol=length(X.index)))
        y_real<-Y[w_size+i,]
        e<-y_real-alasso_predict

        predicts[i]<- alasso_predict
        trues[i]<- y_real
        errors[i]<- e
      }
    }
  results<-Metrics(pred=predicts, true=trues, h=h)
  results$Data<-Data
  results$Model<-list(Model="Adaptive Lasso", Window=window, Size=w_size, Horizon=h, Preselection=t.select, Update=t.update, Index=y.index, Standardize=standardize, RidgeLambda=lambda.ridge, LassoLambda=lambda.lasso)
  results$Variables<-covariates
  class(results)<-"Maeforecast"
  return(results)

}






















maeforecast.postalasso<-function(data=NULL, w_size=NULL, window="recursive", h=0, standardize=TRUE, lambda.ridge=NULL, lambda.lasso=NULL, y.index=1, t.select=NULL, t.update=FALSE){
  if(is.null(data)|is.null(w_size)){
    stop("Have to provide values for data and w_size.")
  }
  if(window %in% c("recursive", "rolling", "fixed")==F){
    stop("Unsupported forecasting sheme. Has to be either 'recursive', 'rolling' or 'fixed'.")
  }

  Data = as.matrix(data)
  if(!is.null(colnames(data))){
    varnames<-colnames(data)[-y.index]
  }else{
    varnames<-seq(from=1, to=dim(data)[2], by=1)[-y.index]
  }

  X = matrix(Data[1:(dim(Data)[1]-h),-y.index],nrow = (dim(Data)[1]-h))
  Y = matrix(Data[(1+h):dim(Data)[1],y.index],nrow = (dim(Data)[1]-h))

  w_size = w_size
  n_windows = dim(Y)[1] - w_size
  predicts<-c()
  trues<-c()
  errors<-c()
  covariates<-list()

  if(t.update==F){
    if(class(t.select)=="numeric"){
      ranks<-t.rank(Data[1:(w_size), ], y.index=y.index, h=h)
      X.index<-ranks[1:t.select,1]
      varnames.selected<-varnames[X.index]
    }else{
      X.index<-1:ncol(X)
      varnames.selected<-varnames
    }
  }

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

      if(t.update){
        if(class(t.select)=="numeric"){

          ranks<-t.rank(Data[1:(w_size+i-1), ], y.index=y.index, h=h)

           X.index<-ranks[1:t.select,1]

           varnames.selected<-varnames[X.index]

           }else{

             X.index<-1:ncol(X)

             varnames.selected<-varnames

           }
        }

      ridge_cv <- cv.glmnet(x = X[1:(w_size + i - 1), X.index], y = Y[1:(w_size + i - 1),],
                            type.measure = "mse",
                            nfold = 10,
                            alpha = 0,
                            standardize=standardize,
                            lambda=lambda.ridge)
      best_ridge_coef <- as.numeric(coef(ridge_cv, s = ridge_cv$lambda.1se))[-1]

      alasso.mod <- glmnet(x = X[1:(w_size + i - 1), X.index], y = Y[1:(w_size + i - 1),],
                           alpha = 1,
                           penalty.factor = 1/abs(best_ridge_coef),
                           standardize=standardize)

      alasso_cv <- cv.glmnet(x = X[1:(w_size + i - 1), X.index], y = Y[1:(w_size + i - 1),],
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
      covariates[i]<-list(varnames.selected[nonzero_index])

      if (length(nonzero_index) !=0){
        nonzero_index <-as.vector(nonzero_index)
        xregs <- X[1:(w_size + i - 1), X.index][, nonzero_index]
        newxregs <-matrix(matrix(X[(w_size + i), X.index], ncol=length(X.index), nrow=1)[, nonzero_index], nrow=1)

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

      if(t.update){
        if(class(t.select)=="numeric"){

          ranks<-t.rank(Data[i:(w_size+i-1), ], y.index=y.index, h=h)

          X.index<-ranks[1:t.select,1]

          varnames.selected<-varnames[X.index]

          }else{

            X.index<-1:ncol(X)

            varnames.selected<-varnames

          }
        }


      ridge_cv <- cv.glmnet(x = X[i:(w_size + i - 1), X.index], y = Y[i:(w_size + i - 1),],
                            type.measure = "mse",
                            nfold = 10,
                            alpha = 0,
                            standardize=standardize,
                            lambda=lambda.ridge)
      best_ridge_coef <- as.numeric(coef(ridge_cv, s = ridge_cv$lambda.1se))[-1]

      alasso.mod <- glmnet(x = X[i:(w_size + i - 1), X.index], y = Y[i:(w_size + i - 1),],
                           alpha = 1,
                           penalty.factor = 1/abs(best_ridge_coef),
                           standardize=standardize)

      alasso_cv <- cv.glmnet(x = X[i:(w_size + i - 1), X.index], y = Y[i:(w_size + i - 1),],
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
      covariates[i]<-list(varnames.selected[nonzero_index])

      if (length(nonzero_index) !=0){
        nonzero_index <-as.vector(nonzero_index)
        xregs <- X[i:(w_size + i - 1), X.index][, nonzero_index]
        newxregs <-matrix(matrix(X[(w_size + i), X.index], ncol=length(X.index), nrow=1)[, nonzero_index], nrow=1)

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

    if(t.update){
      if(class(t.select)=="numeric"){

        ranks<-t.rank(Data[1:w_size, ], y.index=y.index, h=h)

        X.index<-ranks[1:t.select,1]

        varnames.selected<-varnames[X.index]

        }else{

          X.index<-1:ncol(X)

          varnames.selected<-varnames
        }
      }

    ridge_cv <- cv.glmnet(x = X[1:w_size, X.index], y = Y[1:w_size,],
                          type.measure = "mse",
                          nfold = 10,
                          alpha = 0,
                          standardize=standardize,
                          lambda=lambda.ridge)
    best_ridge_coef <- as.numeric(coef(ridge_cv, s = ridge_cv$lambda.1se))[-1]

    alasso.mod <- glmnet(x = X[1:w_size, X.index], y = Y[1:w_size,],
                         alpha = 1,
                         penalty.factor = 1/abs(best_ridge_coef),
                         standardize=standardize)

    alasso_cv <- cv.glmnet(x = X[1:w_size, X.index], y = Y[1:w_size,],
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
    covariates<-list(varnames.selected[nonzero_index])

    if (length(nonzero_index) !=0){
      nonzero_index <-as.vector(nonzero_index)
      Data_uni = ts(Y, frequency=12)
      xregs <- X[1:w_size, X.index][, nonzero_index]
      AR_alasso <-Arima(Data_uni[1:w_size], order=c(1,0,0), xreg =xregs)
      for(i in 1:n_windows){
        newxregs <-matrix(matrix(X[(w_size + i), X.index], ncol=length(X.index), nrow=1)[, nonzero_index], nrow=1)

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
  results<-Metrics(pred=predicts, true=trues, h=h)
  results$Data<-Data
  results$Model<-list(Model="Post Adaptive Lasso", Window=window, Size=w_size, Horizon=h, Preselection=t.select, Update=t.update, Index=y.index, Standardize=standardize, RidgeLambda=lambda.ridge, LassoLambda=lambda.lasso)
  results$Variables<-covariates
  class(results)<-"Maeforecast"
  return(results)
}
















maeforecast.postnet<-function(data=NULL, w_size=NULL, window="recursive", h=0, pred=NULL, standardize=TRUE, alphas=c(0.2, 0.8, 0.02), y.index=1, t.select=NULL, t.update=FALSE){
  if(is.null(data)|is.null(w_size)){
    stop("Have to provide values for data and w_size.")
  }
  if(window %in% c("recursive", "rolling", "fixed")==F){
    stop("Unsupported forecasting sheme. Has to be either 'recursive', 'rolling' or 'fixed'.")
  }


  Data = as.matrix(data)
  if(!is.null(colnames(data))){
    varnames<-colnames(data)[-y.index]
  }else{
    varnames<-seq(from=1, to=dim(data)[2], by=1)[-y.index]
  }

  X = matrix(Data[1:(dim(Data)[1]-h),-y.index],nrow = (dim(Data)[1]-h))
  Y = matrix(Data[(1+h):dim(Data)[1],y.index],nrow = (dim(Data)[1]-h))

  if(is.null(pred)){
    if(!is.null(t.select)){
      pred=t.select
    }else{
      pred=w_size
    }
  }else{
    pred=pred
  }

  w_size = w_size
  n_windows = dim(Y)[1] - w_size
  predicts<-c()
  trues<-c()
  errors<-c()
  covariates<-list()

  if(t.update==F){
    if(class(t.select)=="numeric"){
      ranks<-t.rank(Data[1:(w_size), ], y.index=y.index, h=h)
      X.index<-ranks[1:t.select,1]
      varnames.selected<-varnames[X.index]
    }else{
      X.index<-1:ncol(X)
      varnames.selected<-varnames
    }
  }

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

      if(t.update){
        if(class(t.select)=="numeric"){

          ranks<-t.rank(Data[1:(w_size+i-1), ], y.index=y.index, h=h)

          X.index<-ranks[1:t.select,1]

          varnames.selected<-varnames[X.index]

          }else{

            X.index<-1:ncol(X)

            varnames.selected<-varnames

          }
        }

      Xsis <- SIS.gaussian(X[1:(w_size + i - 1), X.index],
                           matrix(Y[1:(w_size + i - 1),1]),
                           pred=pred,
                           scale=standardize)

      aenet.fit <- aenet(Xsis$Xs, matrix(Y[1:(w_size + i - 1),1]),
                         alphas=alphas, seed = 10)

      nonzero_index = which(!coef(aenet.fit) == 0)
      covariates[i]<-list(varnames.selected[nonzero_index])

      if (length(nonzero_index) !=0){
        nonzero_index <-as.vector(nonzero_index)
        xregs <- X[1:(w_size + i - 1), X.index][, Xsis$Index[1:pred][nonzero_index]]
        newxregs <-matrix(matrix(X[(w_size + i), X.index], ncol=length(X.index), nrow=1)[, Xsis$Index[1:pred][nonzero_index]], nrow=1)

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

      if(t.update){
        if(class(t.select)=="numeric"){

          ranks<-t.rank(Data[i:(w_size+i-1), ], y.index=y.index, h=h)

          X.index<-ranks[1:t.select,1]

          varnames.selected<-varnames[X.index]

          }else{

            X.index<-1:ncol(X)

            varnames.selected<-varnames

          }
        }

      Xsis <- SIS.gaussian(X[i:(w_size + i - 1), X.index],
                           matrix(Y[i:(w_size + i - 1),1]),
                           pred=pred,
                           scale=standardize)

      aenet.fit <- aenet(Xsis$Xs, matrix(Y[i:(w_size + i - 1),1]),
                         alphas=alphas, seed = 10)

      nonzero_index = which(!coef(aenet.fit) == 0)
      covariates[i]<-list(varnames.selected[nonzero_index])

      if (length(nonzero_index) !=0){
        nonzero_index <-as.vector(nonzero_index)
        xregs <- X[i:(w_size + i - 1), X.index][, Xsis$Index[1:pred][nonzero_index]]
        newxregs <-matrix(matrix(X[(w_size + i), X.index], ncol=length(X.index), nrow=1)[, Xsis$Index[1:pred][nonzero_index]], nrow=1)

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

    if(t.update){
      if(class(t.select)=="numeric"){

        ranks<-t.rank(Data[1:w_size, ], y.index=y.index, h=h)

        X.index<-ranks[1:t.select,1]

        varnames.selected<-varnames[X.index]

        }else{

          X.index<-1:ncol(X)

          varnames.selected<-varnames

        }
      }

    Xsis <- SIS.gaussian(X[1:w_size, X.index],
                         matrix(Y[1:w_size,1]),
                         pred=pred,
                         scale=standardize)

    aenet.fit <- aenet(Xsis$Xs, matrix(Y[1:w_size,1]),
                       alphas=alphas, seed = 10)

    nonzero_index = which(!coef(aenet.fit) == 0)
    covariates<-list(varnames.selected[nonzero_index])


    if (length(nonzero_index) !=0){
      nonzero_index <-as.vector(nonzero_index)
      xregs <- X[1:w_size, X.index][, Xsis$Index[1:pred][nonzero_index]]
      Data_uni = ts(Y, frequency=12)
      AR_net <-Arima(Data_uni[1:w_size], order=c(1,0,0), xreg =xregs)
      for(i in 1:n_windows){
        newxregs <-matrix(matrix(X[(w_size + i), X.index], ncol=length(X.index), nrow=1)[, Xsis$Index[1:pred][nonzero_index]], nrow=1)
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
  results<-Metrics(pred=predicts, true=trues, h=h)
  results$Data<-Data
  results$Model<-list(Model="Post Adaptive ElasticNet", Window=window, Size=w_size, Horizon=h, Preselection=t.select, Update=t.update, Index=y.index, Standardize=standardize, Alphas=alphas, Preditors=pred)
  results$Variables<-covariates
  class(results)<-"Maeforecast"
  return(results)
}








maeforecast.dfm2<-function(data=NULL, w_size=NULL, h=0, window="recursive", factor.num=3, method="two-step", clustor.type="partitional", y.index=1, t.select=NULL, t.update=FALSE){
  if(is.null(data)|is.null(w_size)){
    stop("Have to provide values for data and w_size.")
  }
  if(window %in% c("recursive", "rolling", "fixed")==F){
    stop("Unsupported forecasting sheme. Has to be either 'recursive', 'rolling' or 'fixed'.")
  }

  Data = as.matrix(data)
  if(!is.null(colnames(data))){
    varnames<-colnames(data)[-y.index]
  }else{
    varnames<-seq(from=1, to=dim(data)[2], by=1)[-y.index]
  }

  X = matrix(Data[1:(dim(Data)[1]-h),-y.index],nrow = (dim(Data)[1]-h))
  Y = matrix(Data[(1+h):dim(Data)[1],y.index],nrow = (dim(Data)[1]-h))

  w_size = w_size
  n_windows = dim(Y)[1] - w_size
  predicts<-c()
  trues<-c()
  errors<-c()

  if(t.update==F){
    if(class(t.select)=="numeric"){
      ranks<-t.rank(Data[1:(w_size), ], y.index=y.index, h=h)
      X.index<-ranks[1:t.select,1]
      varnames.selected<-varnames[X.index]
    }else{
      X.index<-1:ncol(X)
      varnames.selected<-varnames
    }
  }

  suppressMessages(require(nowcasting))
  suppressMessages(require(forecast))

  if(window=="recursive"){
    for(i in 1:n_windows){

      if(t.update){
        if(class(t.select)=="numeric"){

          ranks<-t.rank(Data[1:(w_size+i-1), ], y.index=y.index, h=h)

          X.index<-ranks[1:t.select,1]

          varnames.selected<-varnames[X.index]

          }else{

            X.index<-1:ncol(X)

            varnames.selected<-varnames

          }
        }

      allfactors<-clust.factor(X[1:(w_size + i), X.index], fac.num=factor.num, method=method, clustor.type=clustor.type)
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

      if(t.update){
        if(class(t.select)=="numeric"){

          ranks<-t.rank(Data[i:(w_size+i-1), ], y.index=y.index, h=h)

          X.index<-ranks[1:t.select,1]

          varnames.selected<-varnames[X.index]

          }else{

            X.index<-1:ncol(X)

            varnames.selected<-varnames

          }
        }

      allfactors<-clust.factor(X[i:(w_size + i), X.index], fac.num=factor.num, method=method, clustor.type=clustor.type)
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

    if(t.update){
      if(class(t.select)=="numeric"){

        ranks<-t.rank(Data[1:w_size, ], y.index=y.index, h=h)

        X.index<-ranks[1:t.select,1]

        varnames.selected<-varnames[X.index]

        }else{

          X.index<-1:ncol(X)

          varnames.selected<-varnames

        }
      }

    allfactors<-clust.factor(X[, X.index], fac.num=factor.num, method=method, clustor.type=clustor.type)
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
  results<-Metrics(pred=predicts, true=trues, h=h)
  results$Data<-Data
  results$Model<-list(Model="Dynamic Factor Model 2", Window=window, Size=w_size, Horizon=h, Preselection=t.select, Update=t.update, Index=y.index, Factors=factor.num,  Method=method, Clustor=clustor.type)
  class(results)<-"Maeforecast"
  return(results)
}

















maeforecast.dfm<-function(data=NULL, w_size=NULL, h=0, window="recursive", factor.num=3, y.index=1, t.select=NULL, t.update=FALSE){
  if(is.null(data)|is.null(w_size)){
    stop("Have to provide values for data and w_size.")
  }
  if(window %in% c("recursive", "rolling", "fixed")==F){
    stop("Unsupported forecasting sheme. Has to be either 'recursive', 'rolling' or 'fixed'.")
  }
  if("dynfactoR" %in% rownames(installed.packages())==F){
    stop("Package dynfactorR is not installed. To install, call library(devtools) and then call install_github('rbagd/dynfactoR').")
  }

  Data = as.matrix(data)
  if(!is.null(colnames(data))){
    varnames<-colnames(data)[-y.index]
  }else{
    varnames<-seq(from=1, to=dim(data)[2], by=1)[-y.index]
  }

  X = matrix(Data[1:(dim(Data)[1]-h),-y.index],nrow = (dim(Data)[1]-h))
  Y = matrix(Data[(1+h):dim(Data)[1],y.index],nrow = (dim(Data)[1]-h))

  w_size = w_size
  n_windows = dim(Y)[1] - w_size
  predicts<-c()
  trues<-c()
  errors<-c()

  if(t.update==F){
    if(class(t.select)=="numeric"){
      ranks<-t.rank(Data[1:(w_size), ], y.index=y.index, h=h)
      X.index<-ranks[1:t.select,1]
      varnames.selected<-varnames[X.index]
    }else{
      X.index<-1:ncol(X)
      varnames.selected<-varnames
    }
  }

  suppressMessages(require(nowcasting))
  suppressMessages(require(forecast))
  suppressMessages(require(dynfactoR))

  if(window=="recursive"){
    for(i in 1:n_windows){

      if(t.update){
        if(class(t.select)=="numeric"){

          ranks<-t.rank(Data[1:(w_size+i-1), ], y.index=y.index, h=h)

          X.index<-ranks[1:t.select,1]

          varnames.selected<-varnames[X.index]

          }else{

            X.index<-1:ncol(X)

            varnames.selected<-varnames

          }
        }

      fac.mod<-dynfactoR::dfm(X[1:(w_size + i), X.index], r=factor.num, p=1, q=factor.num)
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

      if(t.update){
        if(class(t.select)=="numeric"){

          ranks<-t.rank(Data[i:(w_size+i-1), ], y.index=y.index, h=h)

          X.index<-ranks[1:t.select,1]

          varnames.selected<-varnames[X.index]

          }else{

            X.index<-1:ncol(X)

            varnames.selected<-varnames

          }
        }

      fac.mod<-dynfactoR::dfm(X[i:(w_size + i), X.index], r=factor.num, p=1, q=factor.num)
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

    if(t.update){
      if(class(t.select)=="numeric"){

        ranks<-t.rank(Data[1:w_size, ], y.index=y.index, h=h)

        X.index<-ranks[1:t.select,1]

        varnames.selected<-varnames[X.index]

        }else{

          X.index<-1:ncol(X)

          varnames.selected<-varnames

        }
      }

    fac.mod<-dynfactoR::dfm(X[, X.index], r=factor.num, p=1, q=factor.num)
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
  results<-Metrics(pred=predicts, true=trues, h=h)
  results$Data<-Data
  results$Model<-list(Model="Dynamic Factor Model 1", Window=window, Size=w_size, Horizon=h, Preselection=t.select, Update=t.update, Index=y.index, Factors=factor.num, Method="two-step", Clustor="None")
  class(results)<-"Maeforecast"
  return(results)
}









maeforecast.rf<-function(data=NULL, w_size=NULL, h=0, window="recursive", ntree=500, replace=TRUE, y.index=1, t.select=NULL, t.update=FALSE){
  if(is.null(data)|is.null(w_size)){
    stop("Have to provide values for data and w_size.")
  }
  if(window %in% c("recursive", "rolling", "fixed")==F){
    stop("Unsupported forecasting sheme. Has to be either 'recursive', 'rolling' or 'fixed'.")
  }

  Data = as.matrix(data)
  if(!is.null(colnames(data))){
    varnames<-colnames(data)[-y.index]
  }else{
    varnames<-seq(from=1, to=dim(data)[2], by=1)[-y.index]
  }

  X = matrix(Data[1:(dim(Data)[1]-h),-y.index],nrow = (dim(Data)[1]-h))
  Y = matrix(Data[(1+h):dim(Data)[1],y.index],nrow = (dim(Data)[1]-h))

  w_size = w_size
  n_windows = dim(Y)[1] - w_size
  predicts<-c()
  trues<-c()

  if(t.update==F){
    if(class(t.select)=="numeric"){
      ranks<-t.rank(Data[1:(w_size), ], y.index=y.index, h=h)
      X.index<-ranks[1:t.select,1]
      varnames.selected<-varnames[X.index]
    }else{
      X.index<-1:ncol(X)
      varnames.selected<-varnames
    }
  }

  suppressMessages(require(randomForest))

  if(window=="recursive"){
    for(i in 1:n_windows){

      if(t.update){
        if(class(t.select)=="numeric"){

          ranks<-t.rank(Data[1:(w_size+i-1), ], y.index=y.index, h=h)

          X.index<-ranks[1:t.select,1]

          varnames.selected<-varnames[X.index]

          }else{

            X.index<-1:ncol(X)

            varnames.selected<-varnames
          }

        }

      X.train<-X[1:(w_size+i-1), X.index]
      Y.train<-Y[1:(w_size+i-1),]
      X.test<-X[(w_size + i), X.index]
      Y.test<-Y[(w_size + i),]
      rf.mod<-randomForest(x=X.train, y=Y.train)
      rf_predict<-predict(rf.mod, newdata=X.test)
      predicts[i]<-rf_predict
      trues[i]<-Y.test
    }
  }else if(window=="rolling"){
    for(i in 1:n_windows){

      if(t.update){
        if(class(t.select)=="numeric"){

          ranks<-t.rank(Data[i:(w_size+i-1), ], y.index=y.index, h=h)

          X.index<-ranks[1:t.select,1]

          varnames.selected<-varnames[X.index]

          }else{

            X.index<-1:ncol(X)

            varnames.selected<-varnames

          }
        }

      X.train<-X[i:(w_size+i-1), X.index]
      Y.train<-Y[i:(w_size+i-1),]
      X.test<-X[(w_size + i), X.index]
      Y.test<-Y[(w_size + i),]
      rf.mod<-randomForest(x=X.train, y=Y.train)
      rf_predict<-predict(rf.mod, newdata=X.test)
      predicts[i]<-rf_predict
      trues[i]<-Y.test
    }
  }else{

    if(t.update){
      if(class(t.select)=="numeric"){

        ranks<-t.rank(Data[1:w_size, ], y.index=y.index, h=h)

        X.index<-ranks[1:t.select,1]

        varnames.selected<-varnames[X.index]

        }else{

          X.index<-1:ncol(X)

          varnames.selected<-varnames

        }
      }

    X.train<-X[i:w_size, X.index]
    Y.train<-Y[i:w_size,]
    for(i in 1:n_windows){
      X.test<-X[(w_size + i), X.index]
      Y.test<-Y[(w_size + i),]
      rf.mod<-randomForest(x=X.train, y=Y.train)
      rf_predict<-predict(rf.mod, newdata=X.test)
      predicts[i]<-rf_predict
      trues[i]<-Y.test
    }
  }
  results<-Metrics(pred=predicts, true=trues, h=h)
  results$Data<-Data
  results$Model<-list(Model="Random Forest", Window=window, Size=w_size, Horizon=h, Preselection=t.select, Update=t.update, Index=y.index, Trees=ntree, Replace=replace)
  class(results)<-"Maeforecast"
  return(results)
}








maeforecast.ar<-function(data=NULL, w_size=NULL, window="recursive", y.index=1, h=0){
  if(is.null(data)|is.null(w_size)){
    stop("Have to provide values for data and w_size.")
  }
  if(window %in% c("recursive", "rolling", "fixed")==F){
    stop("Unsupported forecasting sheme. Has to be either 'recursive', 'rolling' or 'fixed'.")
  }

  Data<-as.matrix(data)
  Data_uni = ts(Data[,y.index], frequency=12)

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




  if(h==0){
    w_size = w_size
    n_windows = length(Data_uni) - w_size
    predicts<-c()
    trues<-c()
    errors<-c()

    if(window=="recursive"){
      for(i in 1:n_windows){
        AR_model <- arima(Data_uni[1:(w_size + i - 1)], order=c(1,0,0))
        AR_predict <- predict(AR_model, 1)
        y_real =Data_uni[w_size + i]

        predicts[i] <- AR_predict$pred
        trues[i] <- y_real
        e<-y_real-AR_predict$pred
        errors[i] <- e
      }
    }else if(window=="rolling"){
      for(i in 1:n_windows){
        AR_model <- arima(Data_uni[i:(w_size + i - 1)], order=c(1,0,0))
        AR_predict <- predict(AR_model, 1)
        y_real =Data_uni[w_size + i]

        predicts[i] <- AR_predict$pred
        trues[i] <- y_real
        e<-y_real-AR_predict$pred
        errors[i] <- e
      }
    }else{
      AR_model <- Arima(Data_uni[1:w_size], order=c(1,0,0))
      predicts<-as.numeric(farima(model=AR_model, test=Data_uni[(w_size+1):length(Data_uni)])$forecasts)
      trues<-Data_uni[(w_size+1):length(Data_uni)]
      errors<-trues-predicts
    }
    results<-Metrics(pred=predicts, true=trues)
    results$Data<-Data
    results$Model<-list(Model="AR", Window=window, Size=w_size, Horizon=h, Index=y.index)
    class(results)<-"Maeforecast"
    return(results)
  }else{
    Y<-as.numeric(Data_uni[(1+h):length(Data_uni)])
    X<-as.numeric(Data_uni[1:(length(Data_uni)-h)])
    w_size = w_size
    n_windows = length(Y) - w_size
    predicts<-c()
    trues<-c()
    errors<-c()

    if(window=="recursive"){
      for(i in 1:n_windows){
        AR_model<-lm(Y[1:(w_size+i-1)]~X[1:(w_size+i-1)])
        AR_predict<-X[(w_size+i)]*as.numeric(coef(AR_model))[2]+as.numeric(coef(AR_model))[1]
        y_real<-Y[(w_size+i)]

        predicts[i] <- AR_predict
        trues[i] <- y_real
        e<-y_real-AR_predict
        errors[i] <- e
      }
    }else if(window=="rolling"){
      for(i in 1:n_windows){
        AR_model<-lm(Y[i:(w_size+i-1)]~X[i:(w_size+i-1)])
        AR_predict<-X[(w_size+i)]*as.numeric(coef(AR_model))[2]+as.numeric(coef(AR_model))[1]
        y_real<-Y[(w_size+i)]

        predicts[i] <- AR_predict
        trues[i] <- y_real
        e<-y_real-AR_predict
        errors[i] <- e
      }
    }else{
      AR_model<-lm(Y[1:w_size]~X[1:w_size])
      for(i in 1:n_windows){
        AR_predict<-X[(w_size+i)]*as.numeric(coef(AR_model))[2]+as.numeric(coef(AR_model))[1]
        y_real<-Y[(w_size+i)]

        predicts[i] <- AR_predict
        trues[i] <- y_real
        e<-y_real-AR_predict
        errors[i] <- e
      }
    }
    results<-Metrics(pred=predicts, true=trues, h=h)
    results$Data<-Data
    results$Model<-list(Model="AR", Window=window, Size=w_size, Horizon=h, Index=y.index)
    class(results)<-"Maeforecast"
    return(results)
  }
}








maeforecast.arimax<-function(data=NULL, w_size=NULL, window="recursive", y.index=1, h=0){
  if(is.null(data)|is.null(w_size)){
    stop("Have to provide values for data and w_size.")
  }
  if(window %in% c("recursive", "rolling", "fixed")==F){
    stop("Unsupported forecasting sheme. Has to be either 'recursive', 'rolling' or 'fixed'.")
  }

  Data = as.matrix(data)

  X = matrix(Data[1:(dim(Data)[1]-h),-y.index],nrow = (dim(Data)[1]-h))
  Y = matrix(Data[(1+h):dim(Data)[1],y.index],nrow = (dim(Data)[1]-h))

  if(window=="recursive"){
    for(i in 1:n_windows){
        xregs <- X[1:(w_size + i - 1), ]
        newxregs <-matrix(X[(w_size + i), ], ncol=ncol(X), nrow=1)

        Data_uni = ts(Y, frequency=12)
        AR_lasso <-Arima(Data_uni[1:(w_size + i - 1)], order=c(1,0,0), xreg =xregs)
        AR_lasso_predict <- predict(AR_lasso, newxreg = newxregs, n.ahead =1)$pred
        y_real =Data_uni[w_size + i]

        predicts[i] <- AR_lasso_predict
        trues[i] <- y_real
        e <- y_real - AR_lasso_predict
        errors[i] <- e
    }
  }else if(window=="rolling"){
    for(i in 1:n_windows){
        xregs <- X[i:(w_size + i - 1), ]
        newxregs <-matrix(X[(w_size + i), ], ncol=ncol(X), nrow=1)

        Data_uni = ts(Y, frequency=12)
        AR_lasso <-Arima(Data_uni[i:(w_size + i - 1)], order=c(1,0,0), xreg =xregs)
        AR_lasso_predict <- predict(AR_lasso, newxreg = newxregs, n.ahead =1)$pred
        y_real =Data_uni[w_size + i]

        predicts[i] <- AR_lasso_predict
        trues[i] <- y_real
        e <- y_real - AR_lasso_predict
        errors[i] <- e
    }
  }else{
      Data_uni = ts(Y, frequency=12)
      xregs <- X[1:w_size, ]
      AR_lasso <-Arima(Data_uni[1:w_size], order=c(1,0,0), xreg =xregs)
      for(i in 1:n_windows){
        newxregs <-matrix(X[(w_size + i), ], ncol=lncol(X), nrow=1)

        AR_lasso_predict <- predict(AR_lasso, newxreg = newxregs, n.ahead =1)$pred
        y_real =Data_uni[w_size + i]

        predicts[i] <- AR_lasso_predict
        trues[i] <- y_real
        e <- y_real - AR_lasso_predict
        errors[i] <- e
      }
  }
}



maeforecast.rw<-function(data=NULL, w_size=NULL, window="recursive", y.index=1, h=0){
  if(is.null(data)|is.null(w_size)){
    stop("Have to provide values for data and w_size.")
  }
  if(window %in% c("recursive", "rolling", "fixed")==F){
    stop("Unsupported forecasting sheme. Has to be either 'recursive', 'rolling' or 'fixed'.")
  }

  RW<-function(pred=NULL, true=NULL, h=0){
    if(is.null(pred)|is.null(true)){
      stop("Arguments 'pred' and 'true' cannot be ommitted.")
    }
    forecasts<-data.frame(Forecasts=pred, Realized=true)
    forecasts$Errors<-forecasts$Realized-forecasts$Forecasts
    mse<-mean(na.omit(forecasts$Errors)^2)
    if(h==0|h==1){
      forecasts$Success <- ifelse(sign(forecasts$Forecasts) == sign(forecasts$Realized), 1, 0)
      success_ratio = sum(forecasts$Success)/nrow(forecasts)
      forecasts$Forecasted_Direction <- sign(forecasts$Forecasts)
      forecasts$Forecasted_Direction <- ifelse(sign(forecasts$Forecasted_Direction) == 1, 1, 0)
      forecasts$True_Direction <- sign(forecasts$Realized)
      forecasts$True_Direction <- ifelse(sign(forecasts$True_Direction) == 1, 1, 0)
    }else{
      foredir<-vector()
      truedir<-vector()
      foredir[1:(h-1)]<-NA
      truedir[1:(h-1)]<-NA
      for(i in h:length(forecasts$Forecasts)){
        foredir[i]<-forecasts$Forecasts[(i+1-h)]
        truedir[i]<-forecasts$Realized[(i+1-h)]
      }
      forecasts$Success<-ifelse(sign(foredir)==sign(truedir), 1, 0)
      success_ratio = sum(na.omit(forecasts$Success))/nrow(na.omit(forecasts))
      forecasts$Forecasted_Direction <- ifelse(sign(foredir) == 1, 1, 0)
      forecasts$True_Direction <- ifelse(sign(truedir) == 1, 1, 0)
    }
    results<-list(Forecasts=forecasts, MSE=mse, SRatio=success_ratio)
    return(results)
  }

  Data = as.numeric(data[,y.index])
  #predicts<-Data[w_size:(length(Data)-h-1)]
  trues<-Data[(w_size+1+h):length(Data)]
  predicts<-rep(0, length(trues))
  results<-RW(pred=predicts, true=trues, h=h)
  results$Data<-Data
  results$Model<-list(Model="Random Walk", Window=window, Size=w_size, Horizon=h, Index=y.index)
  class(results)<-"Maeforecast"
  return(results)
}



maeforecast<-function(data=NULL, model="ar", w_size=NULL, window="recursive", y.index=1, h=0, ...){
  if(model %in% c("ar", "lasso", "postlasso", "ridge",
                  "alasso", "postalasso", "postnet",
                  "dfm", "dfm2", "rf", "rw")==FALSE){
    stop("Unsupported model type. Refer to help(maeforecast) for a list of supported models.")
  }
  FUN<-paste("maeforecast.", model, sep="")
  FUN<-match.fun(FUN)
  results<-FUN(data=data, w_size=w_size, window=window, y.index=y.index, h=h, ...)
  return(results)
}


summary.Maeforecast <- function(x, digits=7){
  stopifnot(inherits(x, "Maeforecast"))
  cat("\t\n",
      sprintf("Model: %s\n", x$Model$Model),
      sprintf("Forecasting Window: %s\n", x$Model$Window),
      sprintf("Window Size: %s\n", x$Model$Size),
      sprintf("Forecasting Horizon: %s\n", x$Model$Horizon),
      sprintf("Preselection Number: %s\n", ifelse(is.null(x$Model$Preselection), "N/A", x$Model$Preselection)),
      sprintf("MSE: %s\n", base::round(x$MSE, digits=digits)),
      sprintf("Sucess Ratio: %s", base::round(x$SRatio, digits=digits)))
}





