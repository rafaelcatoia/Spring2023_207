############# ---- Conditional distributions for the posterior

### -- Sigma_U
draw.Sigma_u <- function(mu_u,u,P=30){
  #mu_u vector k
  #u matrix nxn
  N=dim(u)[1]
  K=dim(u)[2]
  
  nu <- (N +1) * K + P + 1
  
  matrix_S <- diag(1,K) + mu_u%*%t(mu_u)
  ## initializing loop
  for (i in 1:N){
    matrix_S <- matrix_S + (u[i,]-mu_u)%*%t(u[i,]-mu_u)
  }
  return(LaplacesDemon::rinvwishart(nu,matrix_S))
}

### -- Sigma_V
draw.Sigma_v <- function(mu_v,v,P=30){
  #mu_u vector k
  #u matrix nxn
  K=dim(v)[2]
  nu <- (P +1) * K + P + 1
  
  matrix_S <- diag(1,K) + mu_v%*%t(mu_v)
  ## initializing loop
  for (i in 1:P){
    matrix_S <- matrix_S + (v[i,]-mu_v)%*%t(v[i,]-mu_v)
  }
  return(LaplacesDemon::rinvwishart(nu,matrix_S))
}

### -- sigma2
draw.sigma2 <- function(u,v,X=X){
  N = dim(u)[1]
  P = dim(v)[1]
  shape <- (N*P)/2 + 2
  scale <- 0
  for(i in 1:N){
    for(j in 1:P){
      scale <- scale + (as.numeric(X[i,j]-t(u[i,])%*%v[j,])^2)/2
    }
  }
  return(LaplacesDemon::rinvgamma(n = 1,shape,scale))
}

### mu_u
draw.mu_u <- function(Sigma_u,u){
  N = dim(u)[1]
  S_post = Sigma_u/(N+1)
  bar_u <- colMeans(u)
  mean_post = bar_u*N/(N+1)
  return(mvtnorm::rmvnorm(1,mean_post,S_post)[1,])
}


### mu_v
draw.mu_v <- function(Sigma_v,v){
  P = dim(v)[1]
  S_post = Sigma_v/(P+1)
  bar_v <- colMeans(v)
  mean_post = bar_v*P/(P+1)
  return(mvtnorm::rmvnorm(1,mean_post,S_post)[1,])
}

### u
draw.u <- function(Sigma_u,sigma2,v,mu_u,X = X,K=2){
  ## every row is a mvnorm
  N <- nrow(X)
  P <- dim(v)[1]
  matrix_U <- matrix(NA,nrow = N,ncol=K)
  ## initializing
  
  for (i in 1:N){
    S <- Sigma_u
    mu <- mu_u
    
    for (j in 1:P){
      mu0 =  mu
      Lamb0 =  S
      inv_Sig <- v[j,]%*%t(v[j,])/sigma2
      
      ybar <- rep(0,K)
      ybar[1] <- X[i,j]/v[j,1]
      
      ## updating parameters
      inv_Lamb0 <- solve(Lamb0)
      S <- solve( inv_Lamb0 + inv_Sig)
      mu <- S%*%(inv_Lamb0%*%mu0+inv_Sig%*%ybar)
      mu <- mu[,1]
    }
    matrix_U[i,] <- mvtnorm::rmvnorm(1,mu,S)
  }
  return(matrix_U)
}


### v
draw.v <- function(Sigma_v,sigma2,u,mu_v,X = X,K=2){
  ## every row is a mvnorm
  N <- dim(u)[1]
  P <- ncol(X)
  matrix_V <- matrix(NA,nrow = P,ncol=K)
  ## initializing
  
  for (j in 1:P){
    S <- Sigma_v
    mu <- mu_v
    
    for (i in 1:N){
      mu0 =  mu
      Lamb0 =  S
      inv_Sig <- u[i,]%*%t(u[i,])/sigma2
      
      ybar <- rep(0,K)
      ybar[1] <- X[i,j]/u[i,1]
      
      ## updating parameters
      inv_Lamb0 <- solve(Lamb0)
      S <- solve(inv_Lamb0 + inv_Sig)
      mu <- S%*%(inv_Lamb0%*%mu0+inv_Sig%*%ybar)
      mu <- mu[,1]
    }
    matrix_V[j,] <- mvtnorm::rmvnorm(1,mu,S)
  }
  return(matrix_V)
}