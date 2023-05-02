## check_convergence ------
nGibbs = 500
burnin = 10000
jump = 200
source("/Users/rafaelcatoia/MyDrive/20_UCSC/000_Courses/207-BayesianModels/TH-Quizz1/Checking_Convergence.R")

###### --------- K=1
list.sigma2  <- readRDS(file = 'k1_list.sigma2')
list.Sigma_u <- readRDS(file = 'k1_list.Sigma_u')
list.Sigma_v <- readRDS(file = 'k1_list.Sigma_v')
list.mu_u    <- readRDS(file = 'k1_list.mu_u')
list.mu_v    <- readRDS(file = 'k1_list.mu_v')
list.u       <- readRDS(file = 'k1_list.u')
list.v       <- readRDS(file = 'k1_list.v')

### cleaning the chain
chain_positions <- seq(burnin+1,(nGibbs*jump + burnin),jump)
chain_k1_sigma2 <- list.sigma2[chain_positions] 
chain_k1_Sigma_u <- list.Sigma_u[chain_positions]
chain_k1_Sigma_v <- list.Sigma_v[chain_positions]
chain_k1_mu_u <- list.mu_u[chain_positions] 
chain_k1_mu_v <- list.mu_v[chain_positions] 
chain_k1_u <- list.u[chain_positions]
chain_k1_v <- list.v[chain_positions]

###### --------- K=2
list.sigma2  <- readRDS(file = 'list.sigma2')
list.Sigma_u <- readRDS(file = 'list.Sigma_u')
list.Sigma_v <- readRDS(file = 'list.Sigma_v')
list.mu_u    <- readRDS(file = 'list.mu_u')
list.mu_v    <- readRDS(file = 'list.mu_v')
list.u       <- readRDS(file = 'list.u')
list.v       <- readRDS(file = 'list.v')

### cleaning the chain
chain_positions <- seq(burnin+1,(nGibbs*jump + burnin),jump)
chain_k2_sigma2 <- list.sigma2[chain_positions] 
chain_k2_Sigma_u <- list.Sigma_u[chain_positions]
chain_k2_Sigma_v <- list.Sigma_v[chain_positions]
chain_k2_mu_u <- list.mu_u[chain_positions] 
chain_k2_mu_v <- list.mu_v[chain_positions] 
chain_k2_u <- list.u[chain_positions]
chain_k2_v <- list.v[chain_positions]


check_chain <- function(chain,xlab='',ylab=''){
  par(mfrow=c(1,2))
  plot(unlist(chain),type='l',xlab=xlab,ylab=ylab)
  acf(unlist(chain))
  par(mfrow=c(1,1))
}



## Sigma2
check_chain(chain_k1_sigma2)
check_chain(chain_k2_sigma2)

## Sigma_u
check_chain(chain_k1_Sigma_u)
verifying_convergence_sigma(chain_k2_Sigma_u)

## Sigma_v
check_chain(chain_k1_Sigma_v)
verifying_convergence_sigma(chain_k2_Sigma_v)

## mu_u
check_chain(chain_k1_mu_u)
verifying_convergence_mu(chain_k2_mu_u)

## mu_v
check_chain(chain_k1_mu_v)
verifying_convergence_mu(chain_k2_mu_v)

##u
chain_u_plyr <- plyr::laply(chain_k1_u,function(el){el})
check_chain(chain_u_plyr[,1])
verifying_convergence_latent(chain_k2_u)

##v
chain_v_plyr <- plyr::laply(chain_k1_v,function(el){el})
check_chain(chain_v_plyr[,1])
verifying_convergence_latent(chain_k2_v)






