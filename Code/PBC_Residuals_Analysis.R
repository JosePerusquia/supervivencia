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
res_cs = (pbc_cox$status == 2) - res_mart
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

####################################################################
# Delta-betas
res_db = residuals(cox_pbc, type = "dfbeta")

id = rep(c(1:418),6)

vars = c(rep('Women',418),rep('Edema 0.5', 418),
         rep('Edema 1',418), rep('Age',418),
         rep('Bilirubin',418),rep('Albumin',418))

vals = c(res_db[,1],res_db[,2],
         res_db[,3],res_db[,4],
         res_db[,5],res_db[,6])

db_df=data.frame(x=id,y=vals,name=vars)

ggplot(data=db_df,aes(x=x,y=y,col=name))+
  geom_line()+
  labs(x='',y=expression(Delta[i]*beta[j]),
       color = "Covariate")
####################################################################

####################################################################
# Likelihood displacement
LD = numeric(418)
res_u = residuals(cox_pbc,type='score')
var_beta = cox_pbc$var

for(i in 1:418){
  LD[i]=t(res_u[i,])%*%var_beta%*%res_u[i,]
}

id = c(1:418)
LD_df = data.frame(x=id,y=LD)

ggplot(data=LD_df,aes(x=id,y=y))+
  geom_point()+
  labs(x="",y=expression(LD[i]))

pbc_cox[which(LD>.2),]
LD[which(LD>.2)]
####################################################################

####################################################################
# Schoenfeld residuals 
res_sch = residuals(cox_pbc,type='schoenfeld')

s_res_sch = matrix(nrow = 161,ncol=6)

for(i in 1:161){
  s_res_sch[i,] = t(161*cox_pbc$var%*%res_sch[i,])+cox_pbc$coefficients
}

s_res_sch = as.data.frame(s_res_sch)
s_res_sch = s_res_sch%>%
  rename('Women'=V1,'Edema.5'=V2,'Edema1'=V3,
         'Age'=V4,'Bilirubin'=V5,'Albumin'=V6)%>%
  mutate(time=pbc$time[which(pbc$status==2)])

ggplot(data=s_res_sch,aes(x=time,y=Women))+
  geom_point()+
  geom_smooth(method = "loess", se = F, span = 0.75)+
  labs(x='Survival time',y='Schoenfeld residuals')

ggplot(data=s_res_sch,aes(x=time,y=Edema.5))+
  geom_point()+
  geom_smooth(method = "loess", se = F, span = 0.75)+
  labs(x='Survival time',y='Schoenfeld residuals')

ggplot(data=s_res_sch,aes(x=time,y=Edema1))+
  geom_point()+
  geom_smooth(method = "loess", se = F, span = 0.75)+
  labs(x='Survival time',y='Schoenfeld residuals')

ggplot(data=s_res_sch,aes(x=time,y=Age))+
  geom_point()+
  geom_smooth(method = "loess", se = F, span = 0.75)+
  labs(x='Survival time',y='Schoenfeld residuals')

ggplot(data=s_res_sch,aes(x=time,y=Bilirubin))+
  geom_point()+
  geom_smooth(method = "loess", se = F, span = 0.75)+
  labs(x='Survival time',y='Schoenfeld residuals')

ggplot(data=s_res_sch,aes(x=time,y=Albumin))+
  geom_point()+
  geom_smooth(method = "loess", se = F, span = 0.75)+
  labs(x='Survival time',y='Schoenfeld residuals')

# ZPH test
zph=cox.zph(cox_pbc,terms=F)
####################################################################