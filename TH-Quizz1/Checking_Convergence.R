## For k=2 -------------------------------------------

verifying_convergence_sigma <- function(aux_list){
  chain <-plyr::laply(aux_list, function(el){el})
  par(mfrow=c(1,2))
  plot(chain[,1,1],type = 'l')
  acf(chain[,1,1])
  par(mfrow=c(1,1))
}

verifying_convergence_mu <- function(aux_list){
  chain <-plyr::laply(aux_list, function(el){el})
  par(mfrow=c(1,2))
  plot(chain[,1],type = 'l')
  acf(chain[,1])
  par(mfrow=c(1,1))
}


verifying_convergence_latent <- function(aux_list){
  chain <-plyr::laply(aux_list, function(el){el})
  par(mfrow=c(1,2))
  plot(chain[,1,1],type = 'l')
  acf(chain[,1,1])
  par(mfrow=c(1,1))
}



