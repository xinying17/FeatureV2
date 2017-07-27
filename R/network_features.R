network_features <- function(L='label',data_train,data_test,nf,p,corr,f_type,s,nc)
{
  classes <- unique(data_train$label)

  names(data_train)[colnames(data_train)==L] <- paste("label")
  names(data_test)[colnames(data_test)==L] <- paste("label")

  data_trainm <- data_train[,colnames(data_train)!=L]
  data_testm <- data_test[,colnames(data_test)!=L]

  train_label <- data_train$label
  test_label <- data_test$label

  # feature selection
  if(nf>0) {
    nf = round(min(ncol(data_train),nf))

    # rank feature by ttest
    indx <- rankfeature(L,data_train,classes,nf)
    train_label <- data_train[,colnames(data_train)==L]
    data_trainm <- data_trainm[,indx]
    test_label <- data_test[,colnames(data_test)==L]
    data_testm <- data_testm[,indx]

  }

  # build network for each class
  train_nets <- structure(list(types = character(),
                               featureIDX = list(),
                               nets = list()))
  if(f_type<3){
    aa=1
    for(t in classes){
      class_train_data <- data_trainm[train_label==t,]
      nets <- network_build(class_train_data, p, corr)
      train_nets$types[[aa]] <- t
      train_nets$featureIDX[[aa]] <- colnames(data_train)
      train_nets$nets[[aa]] <- nets
      aa=aa+1
    }

    # new features
    new_train = NULL
    new_test = NULL
    if(f_type==1){ # new features with different power of laplacian matrix
      for(b in 1:length(train_nets$types)){
        nets <- train_nets$nets[[b]]

        r <- eigen(nets$laplacian)
        V <- r$vectors
        lam <- r$values
        lam[lam<0] = 0
        Lmbd = diag(lam ** abs(s))
        if(s<0){
          Lmbd = ginv(Lmbd)
        }
        newL = V %*% Lmbd %*% solve(V)
        new_train <- cbind(new_train,as.matrix(data_trainm) %*% newL)
        new_test <- cbind(new_test,as.matrix(data_testm) %*% newL)
      }
    }

    if(f_type==2){ # single network intergration value
      for(b in 1:length(train_nets$types)){
        nets <- train_nets$nets[[b]]

        r <- eigen(nets$laplacian)
        V <- r$vectors
        lam <- r$values
        lam[lam<0] = 0
        Lmbd = diag(lam ** abs(s))
        newL = V %*% Lmbd %*% solve(V)
        lap_fun <- function(x) {x %*% newL %*% x}
        new_train <- cbind(new_train,apply(as.matrix(data_trainm),1,lap_fun))
        new_test <- cbind(new_test,apply(as.matrix(data_testm),1,lap_fun))
      }
    }
  }

  if(f_type>2){ # subnetwork integration value
    # build network for each class
    aa = 1
    for(t in classes){
      class_train <- data_trainm[train_label==t,]
      clusters <- hclust(dist(t(as.matrix(class_train))),method = "ward.D")
      clusterCut <- cutree(clusters, nc)
      for(i in 1:nc){
        x = data.frame(class_train[,clusterCut==i])
        if(ncol(x)>2){
          nets <- network_build(as.matrix(x), p, corr)
          train_nets$types[[aa]] <- t
          train_nets$featureIDX[[aa]] <- colnames(x)
          train_nets$nets[[aa]] <- nets
          aa = aa+1
        }
      }

    }
    new_train <- matrix(nrow = nrow(data_train),ncol = length(train_nets$types))
    new_test <- matrix(nrow = nrow(data_test),ncol = length(train_nets$types))

    # new train data
    for(b in 1:length(train_nets$types)){
      nets <- train_nets$nets[[b]]
      smooth_value <- smoothness(Lap = nets$laplacian,
                                 data_trainm[,train_nets$featureIDX[[b]]],s)
      new_train[,b] <- smooth_value
    }

    # new test data
    for(b in 1:length(train_nets$types)){
      nets <- train_nets$nets[[b]]
      smooth_value <- smoothness(nets$laplacian,
                                 data_testm[,train_nets$featureIDX[[b]]],s)
      new_test[,b] <- smooth_value
    }
  }


  # remove na and inf
  new_data <- rbind(new_train,new_test)
  new_data <- data.frame(scale(new_data))

  is.na(new_data) <- sapply(new_data, is.infinite)
  is.na(new_data) <- sapply(new_data, is.nan)
  ind_na <- colSums(is.na(new_data))==0

  xx = seq(from=1,to=nrow(new_train),by=1)
  new_train <- new_data[xx,]
  new_test <- new_data[-xx,]


  return(list(new_train = new_train, new_test = new_test, train_label = train_label, test_label = test_label))

}

smoothness <- function(Lap,f,s){
  # f is the function vector
  # s is parameter on Laplacian function
  f <- as.matrix(f)
  require('expm', quietly = TRUE)
  require('matrixcalc',quietly = TRUE)
  # if(is.positive.definite(Lap,tol=0)){
  #lap_fun <- function(x) {x %*% expm(s*logm(as.matrix(Lap))) %*% x}
  r <- eigen(Lap)
  V <- r$vectors
  lam <- r$values
  lam[lam<0] = 0
  Lmbd = diag(lam ** s)
  newL = V %*% Lmbd %*% solve(V)
  lap_fun <- function(x) {x %*% newL %*% x}

  smooth_value <- apply(f, 1, lap_fun)
  return(smooth_value)
}


