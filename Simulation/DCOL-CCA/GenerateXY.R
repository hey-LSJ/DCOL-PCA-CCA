#Define simulation settings
nfac=3
n.simu=1
nx=c(1000) # dimension of X
ny=c(1000) # dimension of Y
Kx=c(20, 40) # number of effective variables in X
Ky=c(20, 40) # number of effective variables in Y per each effective variable in X
links=c(1) # types of link function, 0 is linear only, 1 is half linear half mix, 2 is all mixed
sn=c(3.16) # signal to noise
N=c(200, 1000, 2000, 5000) # sample size
rho=c(0.1)

combos <- expand.grid(nx, ny, Kx, Ky, links, sn, N, rho)  
colnames(combos)<-c("nx", "ny", 'kx', 'ky', 'links', 'sn', 'N',"rho")
to.delete<-which(combos[,3]*combos[,4]>combos[,2])
if(length(to.delete)>0) combos<-combos[-to.delete,]
combos2<-new("list")
for(i in 1:nrow(combos)) combos2[[i]]<-combos[i,]
combos2

#Here, vec is the setting vector
GenerateXY.f <- function(vec, seed = 1, nfac = 3){ #seed = 123
  set.seed(seed)
  library(mvtnorm)
  vec<-as.numeric(vec)
  
  nx<-vec[1] # dimension of X
  ny<-vec[2]  # dimension ofc Y
  kx<-vec[3] # number of effective variables in X
  ky<-vec[4] # number of effective variables in Y per each effective variable in X
  links<-vec[5]# types of link function, 0 is linear only, 1 is half linear half mix, 2 is all mixed
  sn<-vec[6] 
  N<-vec[7]
  rho<-vec[8] 
  
  
  sig<-matrix(rho, ncol=nx, nrow=nx)
  diag(sig) <- 1
  X<-mvtnorm::rmvnorm(N, mean=rep(0,nx), sigma=sig)  ###zero mean
  
  sig <- matrix(rho, ncol=ny, nrow=ny)
  diag(sig)<-1
  Y<-mvtnorm::rmvnorm(N, mean=rep(0,ny), sigma=sig)
  
  X0 <- X
  ##link functions, you can add your own functions here
  link.functions <- function(x,type){
    if (type == 1) y <- runif(1,-2,2)*x
    if (type == 2) y <- abs(x)
    if (type == 3) y <- (2*x^2) 
    if (type == 4) y <- sin(1*(x-0.5)*pi)
    if (type == 5) y <- 1*(x>quantile(x, 0.25) & x<quantile(x, 0.75))
    if (type == 6) y <- cos(0.75*x)^4
    return(y)
  }
  
  for(i in 1:(kx*ky)){
    alphas <- runif(nfac, min=1, max=2)*sample(c(-1,1), nfac, replace=T)
    sample_idx <- sample(1:kx,nfac)
    if (links == 1){#mix of linear and nonlinear
      probabilities <- c(0.15, rep((1-0.15)/5, 5))  #0.15 --> 0.12
      sample_fun <- sample(1:6, nfac, prob = probabilities)
    }else{ ###pure nonlinear
      sample_fun <- sample(2:6, nfac)
      cat('i=',i, 'sample_idx:', sample_idx, '|', 'sample_fun:', sample_fun, '\n')
    }
    
    idx <- 1:nfac
    Z <- sapply(idx, function(idx){alphas[idx]*link.functions(X[,sample_idx[idx]], sample_fun[idx])})
    Y[,i] = apply(Z, 1, sum)
    link.type[[i]] <- sample_fun
  }
  
  ###add noise-m and standardization
  for(i in 1:ncol(X)) {
    X[,i] <-  X[,i] + rnorm(nrow(X), mean=0, sd= sd(X[,i])/sn)
    X[,i]<- (X[,i]-mean(X[,i]))/sd(X[,i]) }
  for(i in 1:ncol(Y)) {
    Y[,i] <- Y[,i] + rnorm(nrow(Y), mean=0, sd= sd(Y[,i])/sn)
    Y[,i]<- (Y[,i]-mean(Y[,i]))/sd(Y[,i]) }

  
  result <- list()
  result$X <- X
  result$Y <- Y
  return(result)  #columns are featrues
}

#To simulate X and Y, run the following code:
simulated.data <- GenerateXY.f(combos2[[1]])
X <- simulated.data$X
Y <- simulated.data$Y
