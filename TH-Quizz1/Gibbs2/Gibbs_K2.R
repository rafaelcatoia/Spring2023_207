library(dplyr) ; library(ggplot2)
mydir = "/Users/rafaelcatoia/MyDrive/20_UCSC/000_Courses/207-BayesianModels/TH-Quizz1/"
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

######## --- Generators -------------------------
source(file = 'Generators_Functions.R')
source(file = 'Checking_Convergence.R')

######## --- New Dir -------------------------
mydir = "/Users/rafaelcatoia/MyDrive/20_UCSC/000_Courses/207-BayesianModels/TH-Quizz1/Gibbs2/"
setwd(mydir)



######## --- Gibbs -------------------------
mydir = "/Users/rafaelcatoia/MyDrive/20_UCSC/000_Courses/207-BayesianModels/TH-Quizz1/"

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
nGibbs = 500
burnin = 10000
jump = 200
nit = nGibbs*jump + burnin
## initial values
list.sigma2[[1]] <- 1/4
list.Sigma_u[[1]] <- diag(1,K)
list.Sigma_v[[1]] <- diag(1,K)
list.mu_u[[1]] <- rep(0,K)
list.mu_v[[1]] <- rep(0,K)
list.u[[1]] <- matrix(1,nrow=N,ncol=K)
list.v[[1]] <- matrix(1,nrow=P,ncol=K)


for(it in 2:nit){
  ## initial values
  list.u[[it]] <- draw.u(
    Sigma_u = list.Sigma_u[[it-1]],
    sigma2 = list.sigma2[[it-1]],
    mu_u = list.mu_u[[it-1]],
    v = list.v[[it-1]],
    X = X,
    K = K)
  
  list.v[[it]] <- draw.v(
    Sigma_v = list.Sigma_v[[it-1]],
    sigma2 = list.sigma2[[it-1]],
    mu_v = list.mu_v[[it-1]],
    u = list.u[[it]],
    X = X,
    K = K)
  
  
  list.sigma2[[it]] <- draw.sigma2(
    u = list.u[[it]],v = list.v[[it]],X = X)
  
  list.Sigma_u[[it]] <- draw.Sigma_u(
    mu_u = list.mu_u[[it-1]],u = list.u[[it]])
  
  list.Sigma_v[[it]] <- draw.Sigma_v(
    mu_v = list.mu_v[[it-1]],v = list.v[[it]])
  
  list.mu_u[[it]] <- draw.mu_u(
    Sigma_u = list.Sigma_u[[it]],u = list.u[[it]])
  
  list.mu_v[[it]] <- draw.mu_v(
    Sigma_v = list.Sigma_v[[it]],v = list.v[[it]])
  
  if(it%%50==0){
    cat(paste('seed=57, iteration ',it,'------------------------------ \n'))
  }
}

### Saving chain ---------------------------
# saveRDS(list.sigma2,file = 'k1_list.sigma2')
# saveRDS(list.Sigma_u,file = 'k1_list.Sigma_u')
# saveRDS(list.Sigma_v,file = 'k1_list.Sigma_v')
# saveRDS(list.mu_u,file = 'k1_list.mu_u')
# saveRDS(list.mu_v,file = 'k1_list.mu_v')
# saveRDS(list.u,file = 'k1_list.u')
# saveRDS(list.v,file = 'k1_list.v')
###



###### ---- Checking convergence 

chain_positions <- seq(burnin+1,(nGibbs*jump + burnin),jump)
chain_sigma2<- list.sigma2[chain_positions] 
chain_Sigma_u <- list.Sigma_u[chain_positions]
chain_Sigma_v <- list.Sigma_v[chain_positions]
chain_mu_u <- list.mu_u[chain_positions] 
chain_mu_v <- list.mu_v[chain_positions] 
chain_u <- list.u[chain_positions]
chain_v <- list.v[chain_positions]

par(mfrow=c(1,2))
plot(unlist(list.sigma2),type='l')
acf(unlist(list.sigma2))

plot(unlist(chain_sigma2),type='l')
acf(unlist(chain_sigma2))
par(mfrow=c(1,1))

verifying_convergence_sigma(list.Sigma_u)
verifying_convergence_sigma(chain_Sigma_u)

verifying_convergence_sigma(list.Sigma_v)
verifying_convergence_sigma(chain_Sigma_v)

verifying_convergence_mu(list.mu_u)
verifying_convergence_mu(chain_mu_u)

verifying_convergence_mu(list.mu_v)
verifying_convergence_mu(chain_mu_v)

verifying_convergence_latent(list.v)
verifying_convergence_latent(chain_v)

verifying_convergence_latent(list.u)
verifying_convergence_latent(chain_u)

