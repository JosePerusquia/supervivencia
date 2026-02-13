##################################################################
# Log-logistic distribution                 
# Author: Jose Antonio Perusquia Cortes
# Afil: Facultad de Ciencias - UNAM
# Module: Survival analysis 
##################################################################

##################################################################
# Required libraries
library(ggplot2)        # Version 3.5.2  
library(ggthemes)       # Version 5.1.0
library(tidyr)          # Version 1.3.1
library(dplyr)          # Version 1.1.4
##################################################################

##################################################################
# Base theme for the plots
theme_plots = theme_minimal()+
  theme(axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))

theme_set(theme_plots)
##################################################################

##################################################################
# Density, distribution function, quantile function and random
# generation for the normal distribution with parameters
# kappa > 0 and theta in R

# Density with default kappa=theta=0
dloglogistic = function(t,kappa = 1, theta = 0){
  num = kappa*(t^(kappa-1))*exp(theta)
  den = (1 + exp(theta)*(t^(kappa)))^2
  return(num/den)
}

# Distribution function with default kappa=theta=1 and lower.tail=T
# otherwise it calculates P(T>t)
ploglogistic = function(t,kappa = 1, theta = 0, lower.tail = T){
  num = exp(theta)*(t^kappa)
  den = 1 + exp(theta)*(t^kappa)
  if(lower.tail){
    return(num/den)
  }else{
    return(1/den)
  }
}

# Quantile function with default values kappa = 1 = theta and
# lower.tail = T otherwise it obtains the value t such that
# P(T>t) = q
qloglogistic = function(q,kappa = 1, theta = 0, lower.tail = T){
  if(lower.tail){
    num = q/(1-q)
    den = exp(theta)
    return((num/den)^(1/kappa))
  }else{
    num = (1-q)/q
    den = exp(theta)
    return((num/den)^(1/kappa))
  }
}

# Random generation for the loglogistic distribution with default 
# kappa = 1 = theta
rloglogistic = function(n, kappa = 1, theta = 0){
  u = runif(n)
  t = ((u/(1-u))*exp(-theta))^(1/kappa)
  return(t)
}
##################################################################

##################################################################
# Loglogistic with kappa = 1 (no finite mean) and theta = 0 
t = seq(0,1,by=.01)

# Density, distribution, survival and hazard
ft = dloglogistic(t)
Ft = ploglogistic(t)
St = 1-Ft
ht = ft/St

# Plots
df_loglogistic = data.frame(t,ft,Ft,St,ht)

# Density
ggplot(data=df_loglogistic,aes(x=t,y=ft))+
  geom_line(col='cyan4')+
  labs(x=expression(t),y=expression(f(t)))

# Distribution
ggplot(data=df_loglogistic,aes(x=t,y=Ft))+
  geom_line(col='darkred')+
  labs(x=expression(t),y=expression(F(t)))

# Survival
ggplot(data=df_loglogistic,aes(x=t,y=St))+
  geom_line(col='blue4')+
  labs(x=expression(t),y=expression(S(t)))

# Hazard
ggplot(data=df_loglogistic,aes(x=t,y=ht))+
  geom_line(col='darkgreen')+
  labs(x=expression(t),y=expression(h(t)))

# Random sample
set.seed(314159)
t = data.frame(x=rloglogistic(250))
ggplot(data = t, aes(x=x,y=after_stat(density)))+
  geom_histogram(col='black',fill='cyan4')+
  labs(x='',y='')
##################################################################

##################################################################
# Log-logistic with kappa = 8 (finite mean) and theta = 0
kappa = 8
theta = 0
t = seq(0,10,by=.01)

# Density, distribution, survival and hazard
ft = dloglogistic(t,kappa = kappa, theta = theta)
Ft = ploglogistic(t,kappa = kappa, theta = theta)
St = 1-Ft
ht = ft/St;ht

# Plots
df_loglogistic = data.frame(t,ft,Ft,St,ht)

# Density
ggplot(data=df_loglogistic,aes(x=t,y=ft))+
  geom_line(col='cyan4')+
  labs(x=expression(t),y=expression(f(t)))+
  coord_cartesian(x=c(0,2.5))

# Distribution
ggplot(data=df_loglogistic,aes(x=t,y=Ft))+
  geom_line(col='darkred')+
  labs(x=expression(t),y=expression(F(t)))+
  coord_cartesian(x=c(0,2.5))

# Survival
ggplot(data=df_loglogistic,aes(x=t,y=St))+
  geom_line(col='blue4')+
  labs(x=expression(t),y=expression(S(t)))+
  coord_cartesian(x=c(0,2.5))

# Hazard
ggplot(data=df_loglogistic,aes(x=t,y=ht))+
  geom_line(col='darkgreen')+
  labs(x=expression(t),y=expression(h(t)))

##################################################################