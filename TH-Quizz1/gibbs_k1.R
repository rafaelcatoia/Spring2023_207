library(dplyr) ; library(ggplot2)
mydir = "/Users/rafaelcatoia/MyDrive/20_UCSC/000_Courses/207-BayesianModels/TH-Quizz1/"
mydir_fig = "/Users/rafaelcatoia/MyDrive/20_UCSC/000_Courses/207-BayesianModels/TH-Quizz1/Figures/"
setwd(mydir)
df_ratings <- data.table::fread("SimMovieRating.csv") %>% data.frame()
df_covariates <- data.table::fread("SimMovieCovariates.csv") %>% data.frame()

df_stack <- df_ratings %>% 
  mutate(User = factor(1:nrow(df_ratings))) %>%
  tidyr::pivot_longer(!User, names_to = "MovieID",
                      values_to = "Rating") %>%
  left_join(df_covariates %>%
              transmute(MovieID=MovieID,
                        Genre=ifelse(Action==1,'Action','Romance'),
                        Reeves=ifelse(Reeves==1,'Yes','No')))


### functions
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



### gibbs
X <- df_ratings- mean(rowMeans(df_ratings))
N <- dim(X)[1]
P <- dim(X)[2] 
K=1

list.Sigma_u  <- list()
list.mu_u <- list()
list.Sigma_v <- list()
list.mu_v <- list()
list.sigma2 <- list()
list.u <- list()
list.v <- list()

set.seed(57)
nGibbs = 100
burnin = 100
jump = 10
nit = nGibbs*jump + burnin
## initial values
list.sigma2[[1]] <- 1
list.Sigma_u[[1]] <- diag(1,K)
list.Sigma_v[[1]] <- diag(1,K)
list.mu_u[[1]] <- rep(0,K)
list.mu_v[[1]] <- rep(0,K)
list.u[[1]] <- matrix(0.1,nrow=N,ncol=K)
list.v[[1]] <- matrix(0.1,nrow=P,ncol=K)


for(it in 2:nit){
  ## initial values
  list.sigma2[[it]] <- draw.sigma2(
    u = list.u[[it-1]],v = list.v[[it-1]],X = X)
  
  list.Sigma_u[[it]] <- draw.Sigma_u(
    mu_u = list.mu_u[[it-1]],u = list.u[[it-1]])
  
  list.Sigma_v[[it]] <- draw.Sigma_v(
    mu_v = list.mu_v[[it-1]],v = list.v[[it-1]])
  
  list.mu_u[[it]] <- draw.mu_u(
    Sigma_u = list.Sigma_u[[it]],u = list.u[[it-1]])
  
  list.mu_v[[it]] <- draw.mu_v(
    Sigma_v = list.Sigma_v[[it]],v = list.u[[it-1]])
  
  list.u[[it]] <- draw.u(
    Sigma_u = list.Sigma_u[[it]],
    sigma2 = list.sigma2[[it]],
    mu_u = list.mu_u[[it]],
    v = list.v[[it-1]],
    X = X,K = K)
  
  list.v[[it]] <- draw.v(
    Sigma_v = list.Sigma_v[[it]],
    sigma2 = list.sigma2[[it]],
    mu_v = list.mu_v[[it]],
    u = list.u[[it]],
    X = X,K = K)
  if(it%%50==0){
    cat(paste('seed=57, iteration ',it,'------------------------------ \n'))
  }
}

saveRDS(list.sigma2,file = 'k1_list.sigma2')
saveRDS(list.Sigma_u,file = 'k1_list.Sigma_u')
saveRDS(list.Sigma_v,file = 'k1_list.Sigma_v')
saveRDS(list.mu_u,file = 'k1_list.mu_u')
saveRDS(list.mu_v,file = 'k1_list.mu_v')
saveRDS(list.u,file = 'k1_list.u')
saveRDS(list.v,file = 'k1_list.v')

## cleaning the chain::
chain_positions <- seq(burnin+1,(nGibbs*jump + burnin),jump)
chain_sigma2<- list.sigma2[chain_positions] 
chain_Sigma_u <- list.Sigma_u[chain_positions]
chain_Sigma_v <- list.Sigma_v[chain_positions]
chain_mu_u <- list.mu_u[chain_positions] 
chain_mu_v <- list.mu_v[chain_positions] 
chain_u <- list.u[chain_positions]
chain_v <- list.v[chain_positions]


##sigma2
par(mfrow=c(3,1))
plot(unlist(chain_sigma2),type='l')
acf(unlist(chain_sigma2))
pacf(unlist(chain_sigma2))

##chain_Sigma_u
plot(unlist(chain_Sigma_u),type='l')
acf(unlist(chain_Sigma_u))
pacf(unlist(chain_Sigma_u))

##chain_Sigma_v
plot(unlist(chain_Sigma_v),type='l')
acf(unlist(chain_Sigma_v))
pacf(unlist(chain_Sigma_v))


##chain_mu_u
plot(unlist(chain_mu_u),type='l')
acf(unlist(chain_mu_u))
pacf(unlist(chain_mu_u))

##chain_mu_v
plot(unlist(chain_mu_v),type='l')
acf(unlist(chain_mu_v))
pacf(unlist(chain_mu_v))

## chain_u
chain_u_matrix <- plyr::laply(chain_u,function(x){x})
colMeans(chain_u_matrix)

## chain_v
chain_v_matrix <- plyr::laply(chain_v,function(x){x})
par(mfrow=c(1,1))
plot(colMeans(chain_v_matrix),colMeans(X))


############## u1
plot(chain_u_matrix[,1],unlist(chain_mu_u))
plot(chain_u_matrix[,2],unlist(chain_mu_u))
plot(chain_u_matrix[,3],unlist(chain_mu_u))


############## v1
plot(chain_v_matrix[,1],unlist(chain_mu_v))
plot(chain_v_matrix[,2],unlist(chain_mu_v))
plot(chain_v_matrix[,3],unlist(chain_mu_v))


############## u'v and sigma - no deu
par(mfrow=c(5,10))
for(i in 1:5){
  for(j in 1:10){
    plot(sqrt(unlist(chain_sigma2)),
         chain_u_matrix[,i] * chain_v_matrix[,j],
         xlab = '',ylab = '',pch=16)
  }
}
par(mfrow=c(1,1))


### natalismo
plot(chain_u_matrix[,1]*chain_v_matrix[,1],type = 'l')
acf(chain_u_matrix[,1]*chain_v_matrix[,1])

plot(X[1,1]-chain_u_matrix[,1]*chain_v_matrix[,1],unlist(chain_sigma2))
##### 


reeves_check <- data.frame(
  mean_post_v = colMeans(chain_v_matrix),
  reeves=ifelse(df_covariates$Reeves==1,'Y','N'),
  genre=ifelse(df_covariates$Action==1,'Action','Romance'))

reeves_check