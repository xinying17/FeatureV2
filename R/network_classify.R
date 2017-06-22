# main function
network_classify <- function(L='label',data_train,data_test,nf,p,corr,f_type,s,nc){

  nc = round(max(1,min(nc,nf/5)))
  newdata <- network_features(L='label',data_train,data_test,nf,p,corr,f_type,s,nc)

  # test
  require('e1071', quietly = TRUE)
  library(e1071)
  if(f_type<4){
    wts <- nrow(newdata$new_train) / table(newdata$train_label)
    model1 <- svm(newdata$new_train,newdata$train_label,type="C-classification",class.weights = wts)
    prediction <- predict(model1, newdata$new_test)
  }
  else{
    data_trainm <- data_train[,colnames(data_train)!=L]
    data_testm <- data_test[,colnames(data_test)!=L]
    # feature selection
    if(nf>0) {
      nf = min(round(10*nf),ncol(data_trainm))
      # rank feature by ttest
      indx <- rankfeature(L,data_train,classes,nf)
      train_label <- data_train[,colnames(data_train)==L]
      data_trainm <- data_trainm[,indx]
      test_label <- data_test[,colnames(data_test)==L]
      data_testm <- data_testm[,indx]
    }

    Data_train = cbind(data_trainm,newdata$new_train)
    Data_test = cbind(data_testm,newdata$new_test)

    wts <- nrow(newdata$new_train) / table(newdata$train_label)
    model1 <- svm(Data_train,newdata$train_label,type="C-classification",class.weights = wts)
    prediction <- predict(model1, Data_test)
  }

  return(list(pred = prediction, acc = sum(prediction==newdata$test_label)/nrow(data_test)))
}

rankfeature <- function(L,data_train,classes,nf){

  ind <- data_train[,colnames(data_train)==L]==classes[1]
  data_train <- data_train[ ,colnames(data_train)!=L] # remove labels
  z <- NULL
  for(i in 1:ncol(data_train)){
    y = data_train[,i]
    y1 = y[ind]
    y2 = y[!ind]
    result <- t.test(y1,y2)
    z[[i]] <- result$statistic
  }

  indx <- c(1:ncol(data_train))
  indx <- indx[order(-z)]
  return(indx[1:nf])
}

