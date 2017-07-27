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
