##################################################################
# Weibull distribution                 
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
# Weibull distribution with constant hazard rate
lambda = 1
gamma = 1
t = seq(0,7.5,by=.01)

# Density, distribution, survival and hazard
ft = dweibull(t,gamma,1/lambda)
Ft = pweibull(t,gamma,1/lambda)
St = 1-Ft
ht = ft/St

# Plots
df_weibull = data.frame(t,ft,Ft,St,ht)

# Density
ggplot(data=df_weibull,aes(x=t,y=ft))+
  geom_line(col='cyan4')+
  labs(x=expression(t),y=expression(f(t)))

# Distribution
ggplot(data=df_weibull,aes(x=t,y=Ft))+
  geom_line(col='darkred')+
  labs(x=expression(t),y=expression(F(t)))

# Survival
ggplot(data=df_weibull,aes(x=t,y=St))+
  geom_line(col='blue4')+
  labs(x=expression(t),y=expression(S(t)))

# Hazard
ggplot(data=df_weibull,aes(x=t,y=ht))+
  geom_line(col='darkgreen')+
  labs(x=expression(t),y=expression(h(t)))
##################################################################

##################################################################
# Weibull distribution with increasing hazard rate
lambda = 1
gamma = 1.5
t = seq(0,5,by=.01)

# Density, distribution, survival and hazard
ft = dweibull(t,gamma,1/lambda)
Ft = pweibull(t,gamma,1/lambda)
St = 1-Ft
ht = ft/St

# Plots
df_weibull = data.frame(t,ft,Ft,St,ht)

# Density
ggplot(data=df_weibull,aes(x=t,y=ft))+
  geom_line(col='cyan4')+
  labs(x=expression(t),y=expression(f(t)))

# Distribution
ggplot(data=df_weibull,aes(x=t,y=Ft))+
  geom_line(col='darkred')+
  labs(x=expression(t),y=expression(F(t)))

# Survival
ggplot(data=df_weibull,aes(x=t,y=St))+
  geom_line(col='blue4')+
  labs(x=expression(t),y=expression(S(t)))

# Hazard
ggplot(data=df_weibull,aes(x=t,y=ht))+
  geom_line(col='darkgreen')+
  labs(x=expression(t),y=expression(h(t)))
##################################################################

##################################################################
# Weibull distribution with decreasing hazard rate
lambda = 1
gamma = .5
t = seq(0,5,by=.01)

# Density, distribution, survival and hazard
ft = dweibull(t,gamma,1/lambda)
Ft = pweibull(t,gamma,1/lambda)
St = 1-Ft
ht = ft/St

# Plots
df_weibull = data.frame(t,ft,Ft,St,ht)

# Density
ggplot(data=df_weibull,aes(x=t,y=ft))+
  geom_line(col='cyan4')+
  labs(x=expression(t),y=expression(f(t)))

# Distribution
ggplot(data=df_weibull,aes(x=t,y=Ft))+
  geom_line(col='darkred')+
  labs(x=expression(t),y=expression(F(t)))

# Survival
ggplot(data=df_weibull,aes(x=t,y=St))+
  geom_line(col='blue4')+
  labs(x=expression(t),y=expression(S(t)))

# Hazard
ggplot(data=df_weibull,aes(x=t,y=ht))+
  geom_line(col='darkgreen')+
  labs(x=expression(t),y=expression(h(t)))
##################################################################