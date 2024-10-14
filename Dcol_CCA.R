####load libraries
library(ggplot2)
library(iCluster)
library(SNFtool)
library(multiway) #for SCA == jSVD
# library(RegularizedSCA)
library(CCA)
library(pls)
library(ggplot2)
library(mvtnorm)
library(mixOmics)
library(MLmetrics)
library(BiocParallel)

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
    # dd<-apply(abs(d),1,sum,na.rm=T)
    # dd<-apply(d^2,1,sum,na.rm=T)
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
  print('parSapply')
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
    print('Number of node:1')
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
sourceCpp("code/matrix_multiplication.cpp")
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

##Visualize a matrix
image.real <- function(mat, title = NA, label = NA) { 
  mat <- t(mat)[,nrow(mat):1]
  image(mat, axes = FALSE, main = title)
  axis(1, at = seq(0, 1, length = nrow(mat)), labels = label)
  axis(2, at = seq(0, 1, length = ncol(mat)), labels = label)
  box() 
}

#############################################
#################CCA related#################
#############################################

#general eigen decomposition
geigen_m <- function (Amat, Bmat, Cmat){
  Bdim <- dim(Bmat)
  Cdim <- dim(Cmat)
  if (Bdim[1] != Bdim[2]) 
    stop("BMAT is not square")
  if (Cdim[1] != Cdim[2]) 
    stop("CMAT is not square")
  p <- Bdim[1]
  q <- Cdim[1]
  s <- min(c(p, q))
  if (max(abs(Bmat - t(Bmat)))/max(abs(Bmat)) > 1e-10) 
    stop("BMAT not symmetric.")
  if (max(abs(Cmat - t(Cmat)))/max(abs(Cmat)) > 1e-10) 
    stop("CMAT not symmetric.")
  Bmat <- (Bmat + t(Bmat))/2
  Cmat <- (Cmat + t(Cmat))/2
  Bfac <- chol(Bmat)
  Cfac <- chol(Cmat)
  Bfacinv <- solve(Bfac)
  Cfacinv <- solve(Cfac)
  Dmat <- eigenMapMatMult(t(Bfacinv), eigenMapMatMult(Amat, Cfacinv))
  if (p >= q) {
    result <- svd(Dmat)
    values <- result$d
    Lmat <- eigenMapMatMult(Bfacinv, result$u)
    Mmat <- eigenMapMatMult(Cfacinv, result$v)
  }
  else {
    result <- svd(t(Dmat))
    values <- result$d
    Lmat <- eigenMapMatMult(Bfacinv, result$v)
    Mmat <- eigenMapMatMult(Cfacinv, result$u)
  }
  geigenlist <- list(values, Lmat, Mmat)
  names(geigenlist) <- c("values", "Lmat", "Mmat")
  return(geigenlist)
}

Dcol_CCA <- function(X, Y, Scale = TRUE, nNodes = 2, k = 5) { #X:p*k, Y:q*k
  ##row: samples col: features
  library(Rfast)
  X.o <- X
  Y.o <- Y
  if (Scale == TRUE) {
    X <- scale(X)
    Y <- scale(Y)
  }
  Xnames <- dimnames(X)[[2]] #feature names of X 
  Ynames <- dimnames(Y)[[2]] #feature names of Y
  ind.names <- dimnames(X)[[1]]
  Cxx <- cova(X, large = TRUE) 
  Cxx <- Cxx + 0.05*diag(diag(Cxx))   ###regularization
  cat('Cxx done!\n')
  Cyy <- cova(Y, large = TRUE)
  Cyy <- Cyy + 0.05*diag(diag(Cyy))
  cat('Cyy done!\n')
  Cxy <- findDCOL(t(X), t(Y), nNodes)   
  Cxy <- getCov(Cxy, X, Y)
  cat('Cxy done!\n')
  res <- geigen_m(Cxy, Cxx, Cyy)
  cat('eigen-decomposition done!')
  names(res) <- c("cor", "xcoef", "ycoef")
  data.r.X <- eigenMapMatMult(X.o, res$xcoef)[,1:k] #lower-dimensional embeddings of X
  data.r.Y <- eigenMapMatMult(Y.o, res$ycoef)[,1:k] #lower-dimensional embeddings of Y
  return(list(Cxx = Cxx, Cyy = Cyy, Cxy = Cxy, cor = res$cor, 
              names = list(Xnames = Xnames, Ynames = Ynames, ind.names = ind.names), 
              xcoef = res$xcoef, ycoef = res$ycoef, data.r.X = data.r.X, data.r.Y = data.r.Y))
}

#classical CCA + Dcol_CCA
CCA_all <- function(X, Y, Scale = TRUE, nNodes = 16) {                #X:p*k, Y:q*k
  library(Rfast)
  #for(i in 1:nrow(X)) X[i,]<-(X[i,]-mean(X[i,],na.rm=T))/sd(X[i,],na.rm=T)
  #for(i in 1:nrow(Y)) Y[i,]<-(Y[i,]-mean(Y[i,],na.rm=T))/sd(Y[i,],na.rm=T)
  X.o <- X
  Y.o <- Y
  if (Scale == TRUE) {
    X <- scale(X)
    Y <- scale(Y)
  }
  print(apply(X, 2, sd)[1:10]+apply(Y, 2, sd)[1:10])
  Xnames <- dimnames(X)[[2]]
  Ynames <- dimnames(Y)[[2]]
  ind.names <- dimnames(X)[[1]]
  
  Cxx <- cora(X, large = TRUE)
  Cxx <- Cxx + 0.05*diag(diag(Cxx))  
  Cyy <- cora(Y, large = TRUE)
  Cyy <- Cyy + 0.05*diag(diag(Cyy))
  
  #CCA
  Cxy <- var(X, Y)
  res <- geigen_m(Cxy, Cxx, Cyy)
  names(res) <- c("cor", "xcoef", "ycoef")
  # scores <- comput((X), (Y), res)
  # res[['scores']] <- scores
  data.r.X <- X.o%*%res$xcoef
  data.r.Y <- Y.o%*%res$ycoef
  res[['data.r.X']] <- data.r.X
  res[['data.r.Y']] <- data.r.Y
  cat('CCA done!\n')
  
  #Dcol CCA
  dcol.Cxy <- findDCOL(t(X), t(Y), nNodes) 
  dcol.Cxy <- getCov(dcol.Cxy, X, Y)
  
  if (dim(dcol.Cxy)[2] < 15){
    print(dcol.Cxy)
  }else{
    print(dcol.Cxy[1:10, 1:10])
    
  }
  dcol.res <- geigen(dcol.Cxy, Cxx, Cyy)
  #dcol.res <- geigen(dcol.res.s <- geigen(sign(Cxy)*dcol.Cxy, Cxx, Cyy)
  names(dcol.res) <- c("dcol.cor", "dcol.xcoef", "dcol.ycoef")
  # dcol.scores <- comput(X, Y, dcol.res)
  # dcol.res[['dcol.scores']] <- dcol.scores
  
  dcol.data.r.X <- eigenMapMatMult(X.o, dcol.res$dcol.xcoef)
  dcol.data.r.Y <- eigenMapMatMult(Y.o, dcol.res$dcol.ycoef)
  # dcol.data.r.X <- X.o%*%dcol.res$dcol.xcoef
  # dcol.data.r.Y <- Y.o%*%dcol.res$dcol.ycoef
  # 
  dcol.res[['dcol.data.r.X']] <- dcol.data.r.X
  dcol.res[['dcol.data.r.Y']] <- dcol.data.r.Y
  cat('Dcol CCA done!\n')
  
  #######save
  return(list(res = res, 
              dcol.res = dcol.res, Cxy = Cxy, dcol.Cxy = dcol.Cxy))
}

###################################
############plot###################
###################################
#input list from CCA_All 
plot.each.cca <- function(cca.all.res, member.info, dataname, Ari = TRUE, dim = 5 ){
  #k is the index of the method being plotted
  if (names(cca.all.res)[1] == 'cor'){
    method <- 'CCA'
  }else{
    method <- 'Dcol_CCA'
  }
  if (Ari == TRUE){
    kmeans.result.x <- kmeans(cca.all.res[[4]][,1:dim], centers = length(table(member.info)))
    ari.x <- round(mclust::adjustedRandIndex(kmeans.result.x$cluster, member.info), 3)
    print(paste(method,'(on projected x)', '|', dataname, ':', ari.x, '| dim =', dim))
    
    kmeans.result.y <- kmeans(cca.all.res[[5]][,1:dim], centers = length(table(member.info)))
    ari.y <- round(mclust::adjustedRandIndex(kmeans.result.y$cluster, member.info), 3)
    print(paste(method,'(on projected y)', '|', dataname, ':', ari.y, '| dim =', dim))
    library(GGally)
    print(ggpairs(as.data.frame(cca.all.res[[4]][,1:dim]), mapping = aes(color = member.info), upper  = 'blank')+
            ggtitle(paste(method,'(on projected x)', '|', dataname, ':', ari.x, '| dim =', dim))+theme(legend.position = 'bottom'))

    print(ggpairs(as.data.frame(cca.all.res[[5]][,1:dim]), mapping = aes(color = member.info), upper  = 'blank')+
        ggtitle(paste(method,'(on projected y)', '|', dataname, ':', ari.y, '| dim =', dim))+theme(legend.position = 'bottom'))
return(c(ari.x, ari.y))
  }
}

