####load libraries
library(ggplot2)
library(Rfast)
###########################################################
#####Utilities used in both DCOL-PCA & DCOL-CCA############
###########################################################
#Ordered distance of rows of a|x
scol.matrix.order <- function(a,x){ # x is the vector, a is the matrix, find ordered distance of rows.of.a|x
  library(Rfast)
  ##when a is a vector
  if (is.null(nrow(a))|isTRUE(nrow(a) == 1)){
    a<-as.vector(a)
    a<-a[order(x)]
    d<-a[2:length(a)]-a[1:(length(a)-1)]
    dd<-sum(d^2,na.rm=T)
  }
  else{
    a<-a[,order(x)]
    d<-a[,2:ncol(a)]-a[,1:(ncol(a)-1)]
    dd <- rowSums(d^2)
  }
  return(dd)
}

###compute the DCOL matrix (in parallel)
findDCOL <-function(a, b, nNodes = 16){
  ##when a and b are two vectors
  if (is.null(nrow(a))|isTRUE(nrow(a) == 1)){
    return(scol.matrix.order(a, b)) #find ordered distance of rows.of.a|x
  }
  
  #parallel
  library(parallel)
  # print('parSapply')
  if (nNodes > 1){ #parallel computing
    if (nNodes > detectCores()-2){
      cat('Set too many nodes')
      nodes <- max(1,detectCores()-2)
    }
    cl <- makeCluster(nNodes)
    clusterExport(cl, c("scol.matrix.order"))
    dcolab <- parSapply(cl, 1:nrow(b), function(i) scol.matrix.order(a, b[i,]))
    dcolba <- t(parSapply(cl, 1:nrow(a), function(j) scol.matrix.order(b, a[j,])))
    stopCluster(cl)
  }
  else{
    print('Number of nodes:1')
    dcolab <- sapply(1:nrow(b), function(i) scol.matrix.order(a, b[i,]))
    dcolba <- t(sapply(1:nrow(a), function(j) scol.matrix.order(b, a[j,])))
  }
  #Retain the smaller dcol value to ensure the dcol matrix is symmetric
  dcol <- dcolba
  dcol.diff <- dcolab-dcolba
  sel <- which(dcol.diff < 0)
  dcol[sel]<- dcolab[sel]
  return(dcol)
}

#Rcpp
library(Rcpp)
sourceCpp("path/to/your/matrix_multiplication.cpp")
getCov <- function(DCOLMatrix, X, Y){  ##X and Y, row: sample col: feature
  ###two vectors
  if (is.null(nrow(DCOLMatrix)) | isTRUE(nrow(DCOLMatrix) == 1)){ 
    print('X and Y are two vectors')
    n <- length(X)
    var_Y <- var(Y)
    print(var_Y)
    value <- sqrt(max(0,1-(DCOLMatrix)/(2*(n-2)*var_Y)))
    return(value)
  }else{
    CovMatrix <- matrix(NA, nrow = dim(DCOLMatrix)[1],ncol = dim(DCOLMatrix)[2])
    (n <- dim(X)[1])
    Var_list <- apply(Y, 2, var)  ##note that this one is the sample variance 
    scale.matrix <- diag(1/(2*(n-2)*Var_list))
    CovMatrix <- 1- eigenMapMatMult(DCOLMatrix, scale.matrix)
    CovMatrix[CovMatrix < 0] <- 0
    CovMatrix <- sqrt(CovMatrix)
    return(CovMatrix)
  }
}

##visualize a matrix
image.real <- function(mat, title = NA, label = NA) { 
  mat <- t(mat)[,nrow(mat):1]
  image(mat, axes = FALSE, main = title)
  axis(1, at = seq(0, 1, length = nrow(mat)), labels = label)
  axis(2, at = seq(0, 1, length = ncol(mat)), labels = label)
  box() 
}

#Version1
Dcol_PCA0 <- function(X, image = 0, k = 4, labels = NA, Scale = TRUE, nNodes = 1){ 
  ##cols of X: features
  ##k: visualize dimension
  X.o <- X  ##used in the final step
  if (Scale == TRUE) {X <- scale(X)}
  # print(apply(X, 2, sd)[1:min(15,dim(X)[2])])
  DcolMatrix <- findDCOL(t(X), t(X), nNodes) #modified here
  cov_D <- getCov(DcolMatrix, X, X) #DCOL-correlation matrix
  # image.real(cov_D, paste('Dcol Covariance,', 'dim = ', dim(cov_D)[1]))
  
  #Eigendecomposition
  library(RSpectra) 
  nPC <- min(100, dim(cov_D)[1])
  cat('nPC:', nPC,'\n')
  eigenDcomp <- eigs_sym(cov_D, nPC, which = "LM")  # "LM" is the default
  par(mfrow = c(2,1))
  plot(eigenDcomp$values, main = paste('top ', nPC, 'eigenvalues'))
  plot(cumsum(eigenDcomp$values/sum(eigenDcomp$values)), main = 'Proportion of variances', 
       xlab = 'PC index', ylab = 'proportions')
  X.proj <- eigenMapMatMult(as.matrix(X.o), eigenDcomp$vectors)
  df <- data.frame(X.proj[,1:min(k, dim(X)[2])])
  if (image == 1){
    if (k <= 2){
      plot(df, col = as.factor(labels), pch = 16, main = 'Embeddings generated by DCOL-PCA' )
    }
    else{
      library(GGally)
      print(ggpairs(as.data.frame(df), mapping = aes(color = labels), upper  = 'blank')+ggtitle('Embeddings generated by DCOL-PCA'))
    }
  }
  return(list(cov_D= cov_D, Evectors = eigenDcomp$vectors, data.r = df))
}

#alternative version of DCOL-PCA
Dcol_PCA <- function(X, image = 0, k = 4, labels = NA, Scale = TRUE, nNodes = 1){ ##cols of X: features
  ##k: vitualize dimension
  X.o <- X  ##used in the final step
  X <- t(X) #cell-cell similarity
  if (Scale == TRUE) {X <- scale(X)}
  DcolMatrix <- findDCOL(t(X), t(X), nNodes)
  cov_D <- getCov(DcolMatrix, X, X)
  library(RSpectra) #for eigendecomposition
  nPC <- min(100, dim(cov_D)[1])
  cat('nPC:', nPC)
  eigenDcomp <- eigs_sym(cov_D, nPC, which = "LM")  # "LM" is the default
  par(mfrow = c(2,1))
  plot(eigenDcomp$values, main = 'Eigenvalues')
  plot(cumsum(eigenDcomp$values/sum(eigenDcomp$values)), main = 'Proportion of variances', 
       xlab = 'PC index', ylab = 'proportions')
  PCs <- eigenDcomp$vectors  ##difference
  ###change axis
  df <- PCs[,1:min(k, dim(cov_D)[2])]
  library(GGally)
  if (image == 1){
    # df <- data.frame(X.o%*%eigenDcomp$vectors[,1:min(k, dim(X)[2])])
    if (k <= 2){
      plot(df,col = as.factor(labels), pch = 16, main = 'Embeddings generated by DCOL-PCA' )
    }
    else{
      cols <- rainbow(length(table(labels)))
      print(ggpairs(as.data.frame(df), mapping = aes(color = labels), upper  = 'blank')+ggtitle('Embeddings generated by DCOL-PCA'))
    }
  }
  return(list(cov_D= cov_D, Evectors = eigenDcomp$vectors, data.r = df))
}

#########Visualize results############
plot.result<- function(reduced.data, group.info, k = 2) {
  df.emd <- as.data.frame(reduced.data[, 1:k])
  colnames(df.emd) <- paste0('factor', 1:k)
  group.info <- as.factor(group.info)
  
  # Calculate ARI
  kmeans.result <- kmeans(reduced.data[, 1:k], centers = length(table(group.info)))
  ari <- round(mclust::adjustedRandIndex(kmeans.result$cluster, group.info), 3)
  
  if (k <= 2) {
    p <- ggplot(df.emd, aes(x = factor1, y = factor2, color = group.info)) +
      geom_point(size = 2) +
      labs(title = paste('ARI =', ari), x = "Factor1", y = "Factor2") +
      scale_color_discrete(name = "Group") +
      theme_bw() +
      theme(plot.title = element_text(size = 12))
     return(p) 
  } else {
    # Use pairs to create a pairwise plot
    pairs(df.emd[, 1:k], col = group.info, pch = 16, cex = 0.8, lower.panel = NULL,
          main = paste('ARI =', ari))
    return(NULL)
  }
}




