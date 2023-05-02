## check_convergence ------
# nGibbs = 500
# burnin = 10000
# jump = 200

list.sigma2  <- readRDS(file = 'k1_list.sigma2')
list.Sigma_u <- readRDS(file = 'k1_list.Sigma_u')
list.Sigma_v <- readRDS(file = 'k1_list.Sigma_v')
list.mu_u    <- readRDS(file = 'k1_list.mu_u')
list.mu_v    <- readRDS(file = 'k1_list.mu_v')
list.u       <- readRDS(file = 'k1_list.u')
list.v       <- readRDS(file = 'k1_list.v')

### sigma2 list.sigma2
chain_positions <- seq(burnin+1,(nGibbs*jump + burnin),jump)
chain_sigma2<- list.sigma2[chain_positions] 
chain_Sigma_u <- list.Sigma_u[chain_positions]
chain_Sigma_v <- list.Sigma_v[chain_positions]
chain_mu_u <- list.mu_u[chain_positions] 
chain_mu_v <- list.mu_v[chain_positions] 
chain_u <- list.u[chain_positions]
chain_v <- list.v[chain_positions]