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
kappa = 1
scale_weibull = lambda^(-1/kappa)
median_weibull = scale_weibull * (log(2))^(1/kappa)
mean_weibull   = scale_weibull * gamma(1 + 1/kappa)
t = seq(0.00001,7.5,by=.01)

# Density, distribution, survival and hazard
ft = dweibull(t,shape=kappa,scale=scale_weibull)
Ft = pweibull(t,shape=kappa,scale=scale_weibull)
St = pweibull(t,shape=kappa,scale=scale_weibull,lower.tail=FALSE)
ht = ft/St

# Plots
df_weibull = data.frame(t,ft,Ft,St,ht)

# Density with the median and the mean
ggplot(data=df_weibull,aes(x=t,y=ft))+
  geom_line(col='cyan4')+
  labs(x=expression(t),y=expression(f(t)),linetype=NULL,colour=NULL)+
  geom_vline(aes(xintercept = median_weibull,linetype = "Median",colour = "Median"))+
  geom_vline(aes(xintercept = mean_weibull,linetype = "Mean",colour = "Mean"))+
  scale_colour_manual(values = c("Median" = "black","Mean" = "darkred"))+
  scale_linetype_manual(values = c("Median" = "dashed","Mean" = "longdash"))

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
kappa = 1.5
scale_weibull = lambda^(-1/kappa)
median_weibull = scale_weibull * (log(2))^(1/kappa)
mean_weibull   = scale_weibull * gamma(1 + 1/kappa)
t = seq(0.00001,5,by=.01)

# Density, distribution, survival and hazard
ft = dweibull(t,shape=kappa,scale=scale_weibull)
Ft = pweibull(t,shape=kappa,scale=scale_weibull)
St = pweibull(t,shape=kappa,scale=scale_weibull,lower.tail=FALSE)
ht = ft/St

# Plots
df_weibull = data.frame(t,ft,Ft,St,ht)

# Density with the median and the mean
ggplot(data=df_weibull,aes(x=t,y=ft))+
  geom_line(col='cyan4')+
  labs(x=expression(t),y=expression(f(t)),linetype=NULL,colour=NULL)+
  geom_vline(aes(xintercept = median_weibull,linetype = "Median",colour = "Median"))+
  geom_vline(aes(xintercept = mean_weibull,linetype = "Mean",colour = "Mean"))+
  scale_colour_manual(values = c("Median" = "black","Mean" = "darkred"))+
  scale_linetype_manual(values = c("Median" = "dashed","Mean" = "longdash"))

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
kappa = .5
scale_weibull = lambda^(-1/kappa)
median_weibull = scale_weibull * (log(2))^(1/kappa)
mean_weibull   = scale_weibull * gamma(1 + 1/kappa)
t = seq(0.00001,5,by=.01)

# Density, distribution, survival and hazard
ft = dweibull(t,shape=kappa,scale=scale_weibull)
Ft = pweibull(t,shape=kappa,scale=scale_weibull)
St = pweibull(t,shape=kappa,scale=scale_weibull,lower.tail=FALSE)
ht = ft/St

# Plots
df_weibull = data.frame(t,ft,Ft,St,ht)

# Density with the median and the mean
ggplot(data=df_weibull,aes(x=t,y=ft))+
  geom_line(col='cyan4')+
  labs(x=expression(t),y=expression(f(t)),linetype=NULL,colour=NULL)+
  geom_vline(aes(xintercept = median_weibull,linetype = "Median",colour = "Median"))+
  geom_vline(aes(xintercept = mean_weibull,linetype = "Mean",colour = "Mean"))+
  scale_colour_manual(values = c("Median" = "black","Mean" = "darkred"))+
  scale_linetype_manual(values = c("Median" = "dashed","Mean" = "longdash"))+
  coord_cartesian(ylim=c(0,10))

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
  labs(x=expression(t),y=expression(h(t)))+
  coord_cartesian(ylim=c(0,10))
##################################################################