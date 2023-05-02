##----- Figures
## Loading chains::
nGibbs = 500
burnin = 10000
jump = 200
source("/Users/rafaelcatoia/MyDrive/20_UCSC/000_Courses/207-BayesianModels/TH-Quizz1/Checking_Convergence.R/")

###### --------- K=1
mydir = "/Users/rafaelcatoia/MyDrive/20_UCSC/000_Courses/207-BayesianModels/TH-Quizz1/Gibbs2/"
setwd(mydir)
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

rm(list.sigma2,list.Sigma_u,list.Sigma_v,list.mu_u,list.mu_v,list.u,list.v)

#############################################################################
# Generating things that were asked
#############################################################################

## loading data::
library(dplyr) ; library(ggplot2)

df_ratings <- data.table::fread("SimMovieRating.csv") %>% data.frame()
df_covariates <- data.table::fread("SimMovieCovariates.csv") %>% data.frame()


### Vizualizing posterior means for V
chain_plyr <- plyr::laply(chain_k1_v,function(el){el})
df_covariates$ExpectedPostVK1 = colMeans(chain_plyr)
chain_plyr <- plyr::laply(chain_k2_v,function(el){el})
expectedPostK2 <- apply(chain_plyr,MARGIN = c(2,3),FUN = mean)
df_covariates$expectedPostVK21 <- expectedPostK2[,1] 
df_covariates$expectedPostVK22 <- expectedPostK2[,2] 

df_stack <- df_ratings %>% 
  mutate(User = factor(1:nrow(df_ratings))) %>%
  tidyr::pivot_longer(!User, names_to = "MovieID",
                      values_to = "Rating") %>%
  left_join(df_covariates %>%
              transmute(MovieID=MovieID,
                        Movie=factor(1:nrow(df_covariates)),
                        Genre=ifelse(Action==1,'Action','Romance'),
                        Reeves=ifelse(Reeves==1,'Yes','No'),
                        ExpecPostV=ExpectedPostVK1,
                        ExpecPostVk1=expectedPostVK21,
                        ExpecPostVk2=expectedPostVK22
                        ))

p1 = df_stack %>% select(Movie,Genre,Reeves,
                         ExpecPostV,ExpecPostVk1,ExpecPostVk2) %>% distinct() %>% 
  mutate(Movie=as.numeric(Movie)) %>% 
  ggplot(aes(x=Movie,y=ExpecPostV,color=Genre,shape=Reeves))+
  geom_point(alpha=0.5,size=5)+theme_bw(base_size = 16)+theme(legend.position = 'bottom')+
  ggtitle('k=1')

p2 = df_stack %>% select(Movie,Genre,Reeves,
                         ExpecPostV,ExpecPostVk1,ExpecPostVk2) %>% distinct() %>% 
  ggplot(aes(x=ExpecPostVk1,y=ExpecPostVk2,color=Genre,shape=Reeves))+
  geom_point(alpha=0.5,size=5)+theme_bw(base_size = 16)+theme(legend.position = 'bottom')+
  ggtitle('k=2')
ggsave(filename = 'v_inference.pdf',path = 'Figures/',
        ggpubr::ggarrange(p1,p2,common.legend = TRUE, legend="bottom",nrow = 1),device = 'pdf',width = 8,height=4,dpi=1000)


######## ----------------------
######## fited vs observed ----------------------
# For k = 2 ----
v_post <- plyr::laply(chain_k1_v,function(el){el})
u_post <- plyr::laply(chain_k1_u,function(el){el})


mean(u_post[,1]*v_post[,3])

X_hat <- matrix(NA,20,30)

for(i in 1:dim(u_post)[2]){
  for(j in 1:dim(v_post)[2]){
    X_hat[i,j] <- mean(u_post[,i]*v_post[,j])
  }
}

X_hat=X_hat+mean(rowMeans(df_ratings))
colnames(X_hat) <- colnames(df_ratings)
X_hat_stacked <- X_hat %>% data.frame() %>% 
  mutate(User = factor(1:nrow(X_hat))) %>%
  tidyr::pivot_longer(!User, names_to = "MovieID",
                      values_to = "FitedX_k1")
df_stack = df_stack %>% left_join(X_hat_stacked)

# For k = 2 ----

v_post <- plyr::laply(chain_k2_v,function(el){el})
u_post <- plyr::laply(chain_k2_u,function(el){el})
X_hat <- matrix(NA,20,30)
vec_prod <- numeric(length = dim(v_post)[1])
for(i in 1:dim(u_post)[2]){
  for(j in 1:dim(v_post)[2]){
    vec_prod <- numeric(length = dim(v_post)[1])
    for(b in 1:dim(v_post)[1]){
      vec_prod[b]<- u_post[b,i,]%*%v_post[b,j,]
    }
    X_hat[i,j] <- mean(vec_prod)
  }
}

X_hat=X_hat+mean(rowMeans(df_ratings))

colnames(X_hat) <- colnames(df_ratings)
X_hat_stacked <- X_hat %>% data.frame() %>% 
  mutate(User = factor(1:nrow(X_hat))) %>%
  tidyr::pivot_longer(!User, names_to = "MovieID",
                      values_to = "FitedX_k2")


df_stack = df_stack %>% left_join(X_hat_stacked)

df_stack %>% select(Rating,FitedX_k1,Genre,Reeves) %>% distinct() %>% 
  ggplot(aes(x=Rating,y=FitedX_k1,color=Genre,shape=Reeves))+
  geom_point(alpha=0.5,size=4)+theme_bw(base_size = 16)+theme(legend.position = 'bottom')+
  geom_abline(slope=1, intercept = 0)+
  facet_wrap(.~Genre)+
  ggtitle('k=2')

df_stack %>% select(Rating,FitedX_k1,Genre,Reeves) %>% distinct() %>% 
  group_by(Reeves) %>% summarise(n())

