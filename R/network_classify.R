# main function
network_classify <- function(L='label',data_train,data_test,nf,p,corr,f_type,s,nc,kern){

  nc = round(max(1,min(nc,nf/5)))
  newdata <- network_features(L='label',data_train,data_test,nf,p,corr,f_type,s,nc)

  # test
  require('e1071', quietly = TRUE)
  library(e1071)
  if(f_type<4){
    wts <- nrow(newdata$new_train) / table(newdata$train_label)
    model1 <- svm(newdata$new_train,newdata$train_label,type="C-classification",class.weights = wts,kernel = kern)
    prediction <- predict(model1, newdata$new_test)
  }
  else{
    data_trainm <- data_train[,colnames(data_train)!=L]
    data_testm <- data_test[,colnames(data_test)!=L]
    # feature selection
    if(nf>0) {
      nf = round(min(ncol(data_train)-1,10*nf))
      # rank feature by ttest
      classes <- unique(data_train$label)
      indx <- rankfeature(L,data_train,classes,nf)
      train_label <- data_train[,colnames(data_train)==L]
      data_trainm <- data_trainm[,indx]
      test_label <- data_test[,colnames(data_test)==L]
      data_testm <- data_testm[,indx]
    }

    Data_train = cbind(data_trainm,newdata$new_train)
    Data_test = cbind(data_testm,newdata$new_test)

    wts <- nrow(newdata$new_train) / table(newdata$train_label)
    model1 <- svm(Data_train,newdata$train_label,type="C-classification",class.weights = wts,kernel = kern)
    prediction <- predict(model1, Data_test)
  }

  return(list(pred = prediction, acc = sum(prediction==newdata$test_label)/nrow(data_test)))
}


