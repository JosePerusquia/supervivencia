####################################################################
# Kaplan-Meier and Cox regression model                    
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
theme_lung = theme_minimal()+
  theme(axis.title.x = element_text(size = 11,face='bold'),
        axis.title.y = element_text(size = 11,face='bold'),
        axis.text.x = element_text(size=10,face='bold'),
        axis.text.y = element_text(size=10,face='bold'))

theme_set(theme_lung)

####################################################################
# NCCTG Lung Cancer Data
lung=lung%>%
  mutate(status=status-1)
glimpse(lung)
####################################################################

####################################################################
# Nonparametric survival estimates

# Kaplan-Meier survival estimate
km = survfit(Surv(time,event=status)~1,data=lung)

# Plot of the survival function using survminer
ggsurvplot(km,ggtheme=theme_lung,legend='none',
           conf.int = T,palette='darkred',
           conf.int.alpha = .3,
           conf.int.fill = "darkred",
           xlab='Days')

# Nelson-Aalen survival estimate
na = survfit(Surv(time,event=status)~1,data=lung,stype=2)

# Plot using survminer
ggsurvplot(na,ggtheme=theme_lung,
           legend='none',conf.int = T,palette = 'navyblue',
           conf.int.alpha = .3,
           conf.int.fill = "navyblue",
           xlab='Days')
####################################################################

####################################################################
# Comparison by sex
lung_sex=lung%>%
  group_by(sex)%>%
  count()%>%
  ungroup()%>%
  mutate(n,Sex=case_when(sex==1~'Male',
                         sex==2~'Female'),.keep='none')

ggplot(data=lung_sex,aes(x=Sex,y=n,fill=Sex))+
  geom_col(show.legend = F,col='black',alpha=.75)+
  geom_text(aes(label = n), vjust = 2, colour = "black") +
  labs(x='',y='')

# Boxplot stratified by sex
ggplot(lung,aes(y=time,x=factor(sex),col=factor(sex))) +
  geom_boxplot(outlier.alpha=.75,show.legend = F)+
  labs(x='',y="Survival time")+
  scale_x_discrete("",labels=c("1" = "M", "2" = "F"))

# Descriptive statistics by genre
lung_sex=lung%>%
  group_by(sex)%>%
  summarise(Mean = mean(time),
            Min = min(time),
            Q1 = quantile(time,.25),
            Median = median(time),
            Q3 = quantile(time,.75),
            Max = max(time),
            StdDev = sd(time))
lung_sex

# Kaplan-Meier survival estimate
km_sex = survfit(Surv(time,event=status)~sex,data=lung)

# Plot of the survival function using survminer
ggsurvplot(km_sex,ggtheme=theme_lung,legend='none',
           conf.int = F,
           conf.int.alpha = .3,
           xlab='Days')

# Log-rank and Peto-Peto tests
survdiff(Surv(time,event=status)~sex,data=lung,rho=0)
survdiff(Surv(time,event=status)~sex,data=lung,rho=1)
####################################################################

####################################################################
# Comparison by ECOG performance 
lung_ecog=lung%>%
  group_by(ph.ecog)%>%
  count()%>%
  ungroup()%>%
  mutate(n,ECOG=case_when(ph.ecog==0~'Asymp.',
                          ph.ecog==1~'Symp.',
                          ph.ecog==2~'Bed <50%',
                          ph.ecog==3~'Bed>50%',
                          ph.ecog==4~'Bedbound'),
         .keep='none')

ggplot(data=lung_ecog,aes(x=ECOG,y=n,fill=ECOG))+
  geom_col(show.legend = F,col='black',alpha=.75)+
  geom_text(aes(label = n), vjust = -.25, colour = "black") +
  labs(x='',y='')

# Omit bed>50% and bedbound
lung_ecog = lung%>% 
  filter(ph.ecog%in%c(0,1,2))

# Boxplot stratified by ECOG performance
ggplot(lung_ecog,aes(y=time,x=factor(ph.ecog),col=factor(ph.ecog))) +
  geom_boxplot(outlier.alpha=.75,show.legend = F)+
  labs(x='',y="Survival time")+
  scale_x_discrete("",labels=c("0" = "Asymp.", "1" = "Symp.",
                               "2" = "Bed < 50%"))

# Descriptive statistics by ECOG performance
lung_ecog2=lung_ecog%>%
  group_by(ph.ecog)%>%
  summarise(Mean = mean(time),
            Min = min(time),
            Q1 = quantile(time,.25),
            Median = median(time),
            Q3 = quantile(time,.75),
            Max = max(time),
            StdDev = sd(time))
lung_ecog2

# Kaplan-Meier survival estimate
km_ecog = survfit(Surv(time,event=status)~ph.ecog,data=lung_ecog)

# Plot of the survival function using survminer
ggsurvplot(km_ecog,ggtheme=theme_lung,legend='none',
           conf.int = F,
           conf.int.alpha = .3,
           xlab='Days')

# Log-rank and Peto-Peto tests
survdiff(Surv(time,event=status)~ph.ecog,data=lung_ecog,rho=0)
survdiff(Surv(time,event=status)~ph.ecog,data=lung_ecog,rho=1)
####################################################################

####################################################################
# Cox proportional hazards model

# Transform variables into factors
lung_cox = lung%>%
  filter(ph.ecog%in%c("0","1","2"))%>%
  mutate(sex=factor(sex),
         ph.ecog=factor(ph.ecog))

# Complete model
cox_lung = coxph(Surv(time,event=status)~age+sex+ph.ecog+
                   +ph.karno+pat.karno+meal.cal+wt.loss,
                 data=lung_cox)
summary(cox_lung)

# Stepwise variable selection
lung_cox2 = lung_cox%>%
  filter(!is.na(meal.cal))%>%
  filter(!is.na(wt.loss))%>%
  filter(!is.na(pat.karno))

cox_lung = coxph(Surv(time,event=status)~age+sex+ph.ecog+
                   +ph.karno+pat.karno+meal.cal+wt.loss,
                 data=lung_cox2)
stepAIC(cox_lung, direction = "both")

# Model with Sex, ECOG, Ph. Karno and wt.loss
cox_lung = coxph(Surv(time,event=status)~sex+ph.ecog+
                   +ph.karno+wt.loss,
                 data=lung_cox2)
summary(cox_lung)

# Model with just sex and ECOG as covariates with the complete
# data and with the data used for variable selection
cox_lung = coxph(Surv(time,event=status)~sex+ph.ecog,
                 data=lung_cox2)
summary(cox_lung)

cox_lung = coxph(Surv(time,event=status)~sex+ph.ecog,
                 data=lung_cox)
summary(cox_lung)

# Likelihood ratio test
loglikrt = 2*(cox_lung$loglik[2]-cox_lung$loglik[1])

# Wald test
betas = cox_lung$coefficients
waldt = betas%*%solve(cox_lung$var)%*%betas
####################################################################
