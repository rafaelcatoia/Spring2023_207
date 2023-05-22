### -------------------- Script THQ-2 -----
library(nlme)
data(BodyWeight)
data <- data.frame(BodyWeight)

### Packages ---------------------------------------
library(ggplot2) ; library(dplyr)
### ------------------------------------------------


### Explanatory Data Analysis
data %>% mutate(Rat=as.integer(Rat)) %>% group_by(Rat,Time,Diet) %>% 
  summarise(Count = n()) %>% arrange(Rat)

## 
df_rats = data %>% mutate(Rat = as.integer(Rat))

df_rats = df_rats %>% 
  left_join(
    df_rats %>%
      filter(Time==1) %>%
      mutate(OriginWeight = weight) %>%
      select(-weight,-Time)
    ) %>% 
  mutate(deltaWeight = weight - OriginWeight)

df_rats %>% ggplot(aes(x=Time,y=weight,color=Diet,shape=Diet,linetype=Diet))+
  geom_point() + geom_line()+
  facet_wrap(Rat~.,ncol = 2,scales='free')+
  theme_minimal()
  
df_rats %>% ggplot(aes(x=Time,y=deltaWeight,color=Diet,shape=Diet,linetype=Diet))+
  geom_point() + geom_line()+
  facet_wrap(Rat~.,ncol = 2)+
  theme_minimal() +
  theme(legend.position = 'bottom')

### ----- Part 1 
df_rats %>% ggplot(aes(x=Time,y=weight,color=Diet,shape=Diet,linetype=Diet))+
  geom_point() + #geom_line()+
  facet_grid(~Diet)+
  theme_minimal() +
  theme(legend.position = 'bottom')

df_rats$Time %>% unique() %>% length()


#### for D=1 
df_rats_D1 = df_rats %>% filter(Diet==1) %>% select(-Diet)
df_rats_D2 = df_rats %>% filter(Diet==2) %>% select(-Diet)
df_rats_D3 = df_rats %>% filter(Diet==3) %>% select(-Diet)

X_D1 = df_rats_D1$Time
X_D1 = bind_cols(Intercept = 1,Time = X_D1)
Y_D1 = df_rats_D1$weight

X_D2 = df_rats_D2$Time
X_D2 = bind_cols(Intercept = 1,Time = X_D2)
Y_D2 = df_rats_D2$weight

X_D3 = df_rats_D3$Time
X_D3 = bind_cols(Intercept = 1,Time = X_D3)
Y_D3 = df_rats_D3$weight
#### -------------------------------------------------------

#### ----- First 

P1_draw_posterior <- function(B=10000,X,Y){
  # X = X_D2
  # Y = Y_D2
  # X = Xfull
  # Y = Yfull
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  n=nrow(X) ; k = ncol(X)
  XtX_inv <- solve(t(X)%*%X)
  beta_hat <- XtX_inv%*%t(X)%*%(Y)
  s2 <- (t(Y - X%*%beta_hat)%*%
    (Y - X%*%beta_hat) ) /(n-k)

  vet_sigma <- LaplacesDemon::rinvchisq(B,df = n-k,scale = s2)
  beta <- matrix(NA,nrow=B,ncol=ncol(X))
  for(i in 1:B){
    varcov<-as.matrix(XtX_inv*vet_sigma[i])
    beta[i,] <- mvtnorm::rmvnorm(1,beta_hat[,1],varcov)
  }
  output <- bind_cols(beta,vet_sigma)
   colnames(output) = c('Beta0','Beta1','Sig2')
  return(output)
}

set.seed(1234)
K1 = P1_draw_posterior(X = X_D1,Y=Y_D1,B = 1000)
K2 = P1_draw_posterior(X = X_D2,Y=Y_D2,B = 1000)
K3 = P1_draw_posterior(X = X_D3,Y=Y_D3,B = 1000)


p1_post <- bind_rows(
  K1 %>% mutate(Diet=1),
  K2 %>% mutate(Diet=2),
  K3 %>% mutate(Diet=3)) 


p1 = p1_post %>% mutate(Diet=factor(Diet)) %>% 
  tidyr::pivot_longer(
    cols=c('Beta0','Beta1','Sig2'),
    names_to = 'Parameter',values_to = 'Posterior') %>% 
  filter(Parameter =='Beta0') %>% 
  ggplot(aes(Posterior,color=Diet))+
  geom_density()+
  #facet_wrap(Diet~Parameter,scales = 'free_y',ncol = 2) +
  theme_minimal() +
  theme(legend.position = 'bottom')+
  ggtitle('Beta0')

p2 = p1_post %>% mutate(Diet=factor(Diet)) %>% 
  tidyr::pivot_longer(
    cols=c('Beta0','Beta1','Sig2'),
    names_to = 'Parameter',values_to = 'Posterior') %>% 
  filter(Parameter =='Beta1') %>% 
  ggplot(aes(Posterior,color=Diet))+
  geom_density()+
  #facet_wrap(Diet~Parameter,scales = 'free_y',ncol = 2) +
  theme_minimal() +
  theme(legend.position = 'bottom')+
  ggtitle('Beta1')

p3 = p1_post %>% mutate(Diet=factor(Diet)) %>% 
  tidyr::pivot_longer(
    cols=c('Beta0','Beta1','Sig2'),
    names_to = 'Parameter',values_to = 'Posterior') %>% 
  filter(Parameter =='Sig2') %>% 
  ggplot(aes(Posterior,color=Diet))+
  geom_density()+
  #facet_wrap(Diet~Parameter,scales = 'free_y',ncol = 2) +
  theme_minimal() +
  theme(legend.position = 'bottom')+
  ggtitle('Sig2')

ggpubr::ggarrange(p1,p2,p3,ncol=3,common.legend = T,legend = 'bottom') +
  theme(legend.position = 'bottom')



#expected posterior
expecPost_K1 <- colMeans(K1)
expecPost_K2 <- colMeans(K2)
expecPost_K3 <- colMeans(K3)


df_rats %>% mutate(Diet=as.factor(Diet)) %>% 
  ggplot(aes(y=weight,x=Time,shape=Diet,color=Diet))+
  geom_point()+
  geom_abline(slope = expecPost_K1[2],
              intercept = expecPost_K1[1],
              color='firebrick',linetype=1)+
  geom_abline(slope = expecPost_K2[2],
              intercept = expecPost_K2[1] ,
              color='forestgreen',linetype=1)+
  geom_abline(slope = expecPost_K3[2],
              intercept = expecPost_K3[1] ,
              color='darkblue',linetype=1)+
  theme_minimal()+theme(legend.position = 'bottom')

##### ------------------------------- Part 2 of the problem

X_D1 = df_rats_D1$Time
X_D1 = bind_cols(Intercept = 1,Time = X_D1)
Y_D1 = df_rats_D1$weight

X_D2 = df_rats_D2$Time
X_D2 = bind_cols(Intercept = 1,Time = X_D2)
Y_D2 = df_rats_D2$weight

X_D3 = df_rats_D3$Time
X_D3 = bind_cols(Intercept = 1,Time = X_D3)
Y_D3 = df_rats_D3$weight

Yfull = matrix(c(Y_D1,Y_D2,Y_D3),ncol=1)

Xfull = matrix(1,nrow=nrow(Yfull))
Xfull = cbind(Xfull,
              c(rep(0,nrow(X_D1)),
                rep(1,nrow(X_D2)),
                rep(0,nrow(X_D2))))
Xfull = cbind(Xfull,
              c(rep(0,nrow(X_D1)),
                rep(0,nrow(X_D2)),
                rep(1,nrow(X_D2))))

Xfull = cbind(Xfull,
              c(df_rats_D1$Time,
                df_rats_D2$Time,
                df_rats_D3$Time))

## ------ Now fitting: 
set.seed(1234)
P2_draw = P1_draw_posterior(X = Xfull,Y=Yfull,B = 1000)
colnames(P2_draw) = c('Alpha','Alpha1','Alpha2','Beta1','Sig2')
plot(P2_draw)


P2_draw %>% tidyr::pivot_longer(
  cols=c('Alpha','Alpha1','Alpha2','Beta1','Sig2'),
  names_to = 'Parameter',values_to = 'Posterior') %>% 
  ggplot(aes(Posterior))+
  geom_density()+
  facet_wrap(Parameter~.,scales = 'free',ncol = 2) +
  theme_minimal() +
  theme(legend.position = 'bottom')

#expected posterior
expecPost <- colMeans(P2_draw)

df_rats %>% mutate(Diet=as.factor(Diet)) %>% 
  ggplot(aes(y=weight,x=Time,shape=Diet,color=Diet))+
  geom_point()+
  geom_abline(slope = expecPost[4],
              intercept = expecPost[1],
              color='firebrick',linetype=1)+
  geom_abline(slope = expecPost[4],
              intercept = expecPost[2] + expecPost[1],
              color='forestgreen',linetype=1)+
  geom_abline(slope = expecPost[4],
              intercept = expecPost[3] + expecPost[1],
              color='darkblue',linetype=1)+
  theme_minimal()+theme(legend.position = 'bottom')
