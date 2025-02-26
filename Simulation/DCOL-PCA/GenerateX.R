#Define simulation settings
n.simu=1
n.genes=c(500, 1000, 2000) 
n.samples = c(500, 1000, 2000) 
n.grps = c(5, 10, 15)
n.fun.types = c(6)
epsilon= c(.2)

combos<-expand.grid(n.samples, n.genes,  n.grps , n.fun.types, epsilon)  

colnames(combos)<-c('n.samples','n.genes', 'n.grps', 'n.fun.types', 'epsilon')
combos$'aver.grp.size' <- floor(combos$n.genes/combos$n.grps)
combos <- combos[order(combos$n.genes),]

combos2<-new("list")
for(i in 1:nrow(combos)) combos2[[i]]<-combos[i,]
combos2

GenerateX.f <- function(seed = 1, n.genes=100, n.samples=100, n.grps=10, aver.grp.size=10, n.fun.types=4, epsilon=0.1, n.depend=0)
{
    set.seed(seed*123)
    link<-function(x, type)
    {
        if(type == 1) return(1*(((x-0.5) %% 2)-1>0))
        if(type == 2) return(cos(2*x))
        if(type == 3) return(x^2)
        if(type == 4) return(abs(x))
        if(type == 5) return(x^3)
        if(type == 6) return(cosh(0.2*x)) 
    }
    
    a<-matrix((runif(n.genes*n.samples)-0.5)*4,ncol=n.samples)
    curr.count<-0
    g<-new("list")
    for(i in 1:n.grps)
    {
        this.size <- aver.grp.size
        if(this.size < 2) this.size<-2
        
        this.mat<-matrix(0, nrow=this.size, ncol=n.samples)
        this.mat[1,]<-(runif(n.samples)-0.5)*4
        for(j in 2:this.size)
        {
            if(n.depend==0)
            {
                this.basis<-c(1, rep(0,j-2))
            }else{
                if(j-1 <= n.depend)
                {
                    this.basis<-rep(1, j-1)
                }else{
                    this.basis<-sample(c(rep(1, n.depend), rep(0,j-1-n.depend)), j-1, replace=F)
                }
                
            }
            if(sum(this.basis) > 0)
            {
                x <- rnorm(n.samples, mean = 0, sd = 0.01)  
                for(k in which(this.basis == 1))
                {
                    x<-x+link(this.mat[k,], sample(n.fun.types,1))*runif(1,min=-1,max=1)
                }
                this.mat[j,]<-x
                this.mat[j,]<-(this.mat[j,]-mean(this.mat[j,]))/sd(this.mat[j,])
            }else{
                this.mat[j,]<-(runif(n.samples)-0.5)*4
            }
        }
        if(n.depend == 0)
        {
            type <- sample(n.fun.types,1)
            # cat('type', type, '\n')
            this.mat[1,]<-link(this.mat[1,], type)
            this.mat[1,]<-(this.mat[1,]-mean(this.mat[1,]))/sd(this.mat[1,])
        }
        
        if(curr.count+this.size <= n.genes)
        {
            a[(curr.count+1):(curr.count+this.size),]<-this.mat
            g[[length(g)+1]]<-(curr.count+1):(curr.count+this.size)
        }
        curr.count<-curr.count + this.size
    }
    a<-a+matrix(rnorm(n.genes*n.samples, sd=epsilon),ncol=n.samples)
    g2<-rep(0, nrow(a))
    for(i in 1:length(g)) g2[g[[i]]]<-i
    
    dcol.pca.simu <-new("list")
    dcol.pca.simu$data <- a
    dcol.pca.simu$grps <- g2
    return(dcol.pca.simu)
}

#To simulate X, run the following code:
current_setting <- combos2[[1]]
simulated.data <-  GenerateX.f(seed = 1, 
                   n.genes = current_setting$n.genes, 
                   n.samples = current_setting$n.samples, 
                   n.grps = current_setting$n.grps, 
                   aver.grp.size = current_setting$aver.grp.size, 
                   n.fun.types = current_setting$n.fun.types, 
                   epsilon = current_setting$epsilon)
X <- simulated.data$data
groups <- simulated.data$grps
dim(X)
