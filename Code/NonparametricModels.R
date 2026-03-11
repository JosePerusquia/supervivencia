####################################################################
# Kaplan-Meier and Cox model for PBC data                
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

# Nelson-Aalen survival estimate
na = survfit(Surv(time,status==2)~1,data=pbc,stype=2)

# Plot using survminer
ggsurvplot(na,ggtheme=theme_pbc,
           legend='none',conf.int = T,palette = 'navyblue',
           conf.int.alpha = .3,
           conf.int.fill = "navyblue",
           xlab='Days')
####################################################################

####################################################################
# Comparison by sex
pbc_sex=pbc%>%
  group_by(sex)%>%
  count()%>%
  ungroup()%>%
  mutate(n,Sex=case_when(sex=='m'~'Male',
                         sex=='f'~'Female'),.keep='none')

ggplot(data=pbc_sex,aes(x=Sex,y=n,fill=Sex))+
  geom_col(show.legend = F,col='black',alpha=.75)+
  geom_text(aes(label = n), vjust = 2, colour = "black") +
  labs(x='',y='')

# Boxplot stratified by sex
ggplot(pbc,aes(y=time,x=sex,col=sex)) +
  geom_boxplot(outlier.alpha=.75,show.legend = F)+
  labs(x='',y="Survival time")+
  scale_x_discrete("",labels=c("1" = "M", "2" = "F"))

# Descriptive statistics by genre
pbc_sex=pbc%>%
  group_by(sex)%>%
  summarise(Mean = mean(time),
            Min = min(time),
            Q1 = quantile(time,.25),
            Median = median(time),
            Q3 = quantile(time,.75),
            Max = max(time),
            StdDev = sd(time))
pbc_sex

# Kaplan-Meier survival estimate
km_sex = survfit(Surv(time,status==2)~sex,data=pbc)

# Plot of the survival function using survminer
ggsurvplot(km_sex,ggtheme=theme_pbc,legend='none',
           conf.int = F,
           conf.int.alpha = .3,
           xlab='Days')

# Log-rank and Peto-Peto tests
survdiff(Surv(time,status==2)~sex,data=pbc,rho=0)
survdiff(Surv(time,status==2)~sex,data=pbc,rho=1)
####################################################################

####################################################################
# Comparison by stage
pbc_stage=pbc%>%
  group_by(stage)%>%
  count()%>%
  ungroup()%>%
  mutate(n,Stage=case_when(stage==1~'Stage I',
                         stage==2~'Stage II',
                         stage==3~'Stage III',
                         stage==4~'Stage IV'),.keep='none')

ggplot(data=pbc_stage,aes(x=Stage,y=n,fill=Stage))+
  geom_col(show.legend = F,col='black',alpha=.75)+
  geom_text(aes(label = n), vjust = 2, colour = "black") +
  labs(x='',y='')

# Boxplot stratified by stage of disease
pbc_stage = pbc%>%
  filter(!is.na(stage))

ggplot(pbc_stage,aes(y=time,x=as.factor(stage),
                     col=as.factor(stage))) +
  geom_boxplot(outlier.alpha=.75,show.legend = F)+
  labs(x='',y="Survival time")

# Descriptive statistics by genre
pbc_stage=pbc%>%
  filter(!is.na(stage))%>%
  group_by(stage)%>%
  summarise(Mean = mean(time),
            Min = min(time),
            Q1 = quantile(time,.25),
            Median = median(time),
            Q3 = quantile(time,.75),
            Max = max(time),
            StdDev = sd(time))
pbc_stage

# Kaplan-Meier survival estimate
km_stage = survfit(Surv(time,status==2)~stage,data=pbc)

# Plot of the survival function using survminer
ggsurvplot(km_stage,ggtheme=theme_pbc,legend='none',
           conf.int = F,
           conf.int.alpha = .3,
           xlab='Days')

# Log-rank and Peto-Peto tests
survdiff(Surv(time,status==2)~stage,data=pbc,rho=0)
survdiff(Surv(time,status==2)~stage,data=pbc,rho=1)
####################################################################

####################################################################
# Comparison by stage (joining I and II)
pbc_mod = pbc%>%filter(!is.na(stage))
pbc_mod$stage[which(pbc_mod$stage==1)]=2

# Boxplot stratified by stage of disease
ggplot(pbc_mod,aes(y=time,x=as.factor(stage),
                   col=as.factor(stage))) +
  geom_boxplot(outlier.alpha=.75,show.legend = F)+
  labs(x='',y="Survival time")

# Descriptive statistics by genre
pbc_stage=pbc_mod%>%
  group_by(stage)%>%
  summarise(Mean = mean(time),
            Min = min(time),
            Q1 = quantile(time,.25),
            Median = median(time),
            Q3 = quantile(time,.75),
            Max = max(time),
            StdDev = sd(time))
pbc_stage

# Kaplan-Meier survival estimate
km_stage = survfit(Surv(time,status==2)~stage,data=pbc_mod)

# Plot of the survival function using survminer
ggsurvplot(km_stage,ggtheme=theme_pbc,legend='none',
           conf.int = F,
           conf.int.alpha = .3,
           xlab='Days')

# Log-rank and Peto-Peto tests
survdiff(Surv(time,status==2)~stage,data=pbc_mod,rho=0)
survdiff(Surv(time,status==2)~stage,data=pbc_mod,rho=1)
####################################################################

####################################################################
# Cox proportional hazards model using age, sex, edema,
# bili and albumin as covariates

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

# Cox model with log for serum bilirunbin
cox_pbc = coxph(Surv(time,status==2)~sex+edema+
                  age+log(bili)+albumin,
                data=pbc_cox)
summary(cox_pbc)
####################################################################

####################################################################
# Cox - Snell 
res_mart = residuals(cox_pbc,type='martingale')
res_cs = (pbc_cox$status == 2) - res_mart
pbc_cox$res_cs = res_cs

# KM using the Cox-Snell residuals as the times
fit_cs = survfit(Surv(res_cs, status == 2) ~ 1, data = pbc_cox)

# Plot cumulative hazard
df = data.frame(x=fit_cs$time,y=-log(fit_cs$surv))

ggplot(data=df,aes(x=x,y=y))+
  geom_line(col='darkgrey')+
  geom_point(size=.5)+
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

# number of observations (n) and covariates (k)
n = nrow(pbc_cox)
k = 6
id = rep(1:n, k)

vars = c(rep('Women',n),rep('Edema 0.5', n),
         rep('Edema 1',n), rep('Age',n),
         rep('Bilirubin',n),rep('Albumin',n))

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
LD = numeric(n)
res_u = residuals(cox_pbc,type='score')
var_beta = cox_pbc$var

for(i in 1:n){
  LD[i]=t(res_u[i,])%*%var_beta%*%res_u[i,]
}

id = c(1:n)
LD_df = data.frame(x=id,y=LD)

ggplot(data=LD_df,aes(x=id,y=y))+
  geom_point()+
  labs(x="",y=expression(LD[i]))

pbc_cox[which(LD>.2),]
LD[which(LD>.2)]
####################################################################

####################################################################
# Schoenfeld residuals 
res_sch = residuals(cox_pbc,type='scaledsch')
res_sch = as.data.frame(res_sch)

res_sch = res_sch%>%
  rename('Women'=sexf,'Edema.5'=edema0.5,'Edema1'=edema1,
         'Age'=age,'Bilirubin'='log(bili)','Albumin'=albumin)%>%
  mutate(time=pbc$time[which(pbc$status==2)])

ggplot(data=res_sch,aes(x=time,y=Women))+
  geom_point()+
  geom_smooth(method = "loess", se = F, span = 0.75)+
  labs(x='Survival time',y='Schoenfeld residuals')

ggplot(data=res_sch,aes(x=time,y=Edema.5))+
  geom_point()+
  geom_smooth(method = "loess", se = F, span = 0.75)+
  labs(x='Survival time',y='Schoenfeld residuals')

ggplot(data=res_sch,aes(x=time,y=Edema1))+
  geom_point()+
  geom_smooth(method = "loess", se = F, span = 0.75)+
  labs(x='Survival time',y='Schoenfeld residuals')

ggplot(data=res_sch,aes(x=time,y=Age))+
  geom_point()+
  geom_smooth(method = "loess", se = F, span = 0.75)+
  labs(x='Survival time',y='Schoenfeld residuals')

ggplot(data=res_sch,aes(x=time,y=Bilirubin))+
  geom_point()+
  geom_smooth(method = "loess", se = F, span = 0.75)+
  labs(x='Survival time',y='Schoenfeld residuals')

ggplot(data=res_sch,aes(x=time,y=Albumin))+
  geom_point()+
  geom_smooth(method = "loess", se = F, span = 0.75)+
  labs(x='Survival time',y='Schoenfeld residuals')

# ZPH test 
zph=cox.zph(cox_pbc,terms=F)
print(zph)

# Another way to plot scaled Schoenfeld residuals and check if the
# PH holds
ggcoxzph(zph[1],title='',ggtheme = theme_pbc,
         xlab='Time',ylab='Coefficient for women')
ggcoxzph(zph[2],title='',ggtheme = theme_pbc,
         xlab='Time',ylab='Coefficient for edema 0.5')
ggcoxzph(zph[3],title='',ggtheme = theme_pbc,
         xlab='Time',ylab='Coefficient for edema 1')
ggcoxzph(zph[4],title='',ggtheme = theme_pbc,
         xlab='Time',ylab='Coefficient for bilirunbin')
ggcoxzph(zph[5],title='',ggtheme = theme_pbc,
         xlab='Time',ylab='Coefficient for albumin')
####################################################################
