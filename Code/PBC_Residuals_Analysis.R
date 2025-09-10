####################################################################
# Residuals analysis Cox proportional hazards model                 
# Author: Jose Antonio Perusquia Cortes
# Afil: Facultad de Ciencias - UNAM
# Module: Survival analysis 
####################################################################

####################################################################
# Libraries
library(here)           # Version 1.0.1
library(ggplot2)        # Version 3.5.2  
library(ggthemes)       # Version 5.1.0
library(dplyr)          # Version 1.1.4
library(readxl)         # Version 1.4.4
library(survival)       # Version 3.8-3
library(survminer)      # Version 0.5.0
library(MASS)           # Version 7.3-60
####################################################################

####################################################################
# Base theme for the plots
theme_pbc = theme_minimal()+
  theme(axis.title.x = element_text(size = 11,face='bold'),
        axis.title.y = element_text(size = 11,face='bold'),
        axis.text.x = element_text(size=10,face='bold'),
        axis.text.y = element_text(size=10,face='bold'))

theme_set(theme_pbc)

####################################################################
# Mayo Clinic Primary Biliary Cholangitis Data
pbc=pbc
glimpse(pbc)
####################################################################

####################################################################
# Nonparametric survival estimates

# Kaplan-Meier survival estimate
km = survfit(Surv(time,status==2)~1,data=pbc)

# Plot of the survival function using survminer
ggsurvplot(km,ggtheme=theme_pbc,legend='none',
           conf.int = T,palette='darkred',
           conf.int.alpha = .3,
           conf.int.fill = "darkred",
           xlab='Days')
####################################################################

####################################################################
# Cox proportional hazards model

# Transform variables into factors
pbc_cox = pbc%>%
  dplyr::select(id,time,status,age,sex,edema,bili,albumin)%>%
  mutate(sex=factor(sex),
         edema=factor(edema))

# Cox model
cox_pbc = coxph(Surv(time,status==2)~sex+edema+
                  age+bili+albumin,
                 data=pbc_cox)
summary(cox_pbc)
####################################################################

####################################################################
# Cox - Snell 
res_mart = residuals(cox_pbc,type='martingale')
res_cs = (pbc_cox$status == 2) - res_mart_mod
pbc_cox$res_cs = res_cs

# qqplot against exponential
expo = rexp(418,1)
df = data.frame(x=expo,y=pbc_cox)

ggplot(data=pbc_cox,aes(sample=res_cs))+
  stat_qq(distribution=qexp,dparams=list(1))+
  labs(x="Theoretical",y="Sample")+
  stat_qq_line(distribution=qexp,dparams=list(1),
               col="red")

# Fit a Cox regression using Cox - Snell residuals
fit_cs = survfit(Surv(res_cs, status == 2) ~ 1, data = pbc_cox)

# Plot cumulative hazard
df = data.frame(x=fit_cs$time,y=-log(fit_cs$surv))

ggplot(data=df,aes(x=x,y=y))+
  geom_step()+
  labs(x="Cox-Snell residulas",y=expression(hat(H[t])))+
  geom_abline(intercept=0,slope=1,col="red")
####################################################################

####################################################################
# Martingale residuals vs continuous variables
pbc_cox$res_mart = res_mart

ggplot(data=pbc_cox,aes(x=age,y=res_mart))+
  geom_point()+
  geom_smooth(method = "loess", se = F, span = 0.75)+
  labs(x='Age',y='Martingale residuals')

ggplot(data=pbc_cox,aes(x=albumin,y=res_mart))+
  geom_point()+
  geom_smooth(method = "loess", se = F, span = 0.75)+
  labs(x='Serum albumin',y='Martingale residuals')

ggplot(data=pbc_cox,aes(x=bili,y=res_mart))+
  geom_point()+
  geom_smooth(method = "loess", se = F, span = 0.75)+
  labs(x='Serum bilirunbin',y='Martingale residuals')
####################################################################

####################################################################
# Deviance residuals
res_dev = residuals(cox_pbc,type='deviance')
pbc_cox$res_dev = res_dev

ggplot(data=pbc_cox,aes(x=id,y=res_dev))+
  geom_point()+
  geom_hline(yintercept=0,col='red')+
  labs(x="",y="Deviance residuals")
####################################################################


