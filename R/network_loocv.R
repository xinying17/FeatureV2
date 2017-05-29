network_loocv_new <- function(data_matrix,L,nf,p,corr,f_type,s,nc){

  names(data_matrix)[colnames(data_matrix)==L] <- paste("label")
  pred <- NULL

  classes <- unique(data_matrix$label)
  for(i in 1:nrow(data_matrix)){
    data_train <- data_matrix[-i,]  # training data
    data_test <- data_matrix[i,] # test data

    result <- network_classify(L='label',data_train,data_test,nf,p,corr,f_type,s,nc)
    pred[[i]] <- result$pred
  }
  predx <- as.numeric(as.factor(pred))
  testx <- as.numeric(as.factor(data_matrix$label))
  accuracy = sum(predx==testx)/nrow(data_matrix)
  return(accuracy)
}
