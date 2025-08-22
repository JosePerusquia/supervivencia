####################################################################
# Exponential and Weibull survival models                    
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

# Censored and death
lung_status=lung%>%
  group_by(status)%>%
  count()%>%
  ungroup()%>%
  mutate(n,Status=case_when(status==0~'Censored',
                          status==1~'Death'),.keep='none')

ggplot(data=lung_status,aes(x=Status,y=n,fill=Status))+
  geom_col(show.legend = F,col='black',alpha=.75)+
  geom_text(aes(label = n), vjust = 2, colour = "black")+
  labs(x='',y='')


# Histogram of survival times
ggplot(data=lung,aes(x=time,y=after_stat(density)))+
  geom_histogram(fill='lightblue3',col='black',
                 breaks=hist(lung$time,plot=F)$breaks)+
  labs(x='Survival time',y='')

# Empirical quantiles of the survival times
mean(lung$time)
quantile(lung$time)
####################################################################

####################################################################
# Exponential model
n = dim(lung)[1]
num_death = sum(lung$status)
sum_times = sum(lung$time)

# Maximum likelihood estimator for lambda
lambda_emv = num_death/sum_times;lambda_emv

# Density
t = seq(0:max(lung$time))
ft_exp = lambda_emv*exp(-lambda_emv*t)
density_exp = data.frame(t,ft_exp)

ggplot(data=lung,aes(x=time,y=after_stat(density)))+
  geom_histogram(fill='lightblue3',col='black',
                 breaks=hist(lung$time,plot=F)$breaks)+
  geom_line(data=density_exp,aes(x=t,y=ft_exp))+  
  labs(x='Survival time',y='')

# Risk function
ht_exp = lambda_emv
risk_exp = data.frame(t,ht_exp)
ggplot(data=risk_exp,aes(x=t,y=ht_exp))+
  geom_line()+
  labs(x=expression(t),y=expression(h(t)))

# Survival function
St_exp = exp(-lambda_emv*t)
surv_exp = data.frame(t,St_exp)
ggplot(data=surv_exp,aes(x=t,y=St_exp))+
  geom_line()+
  labs(x=expression(t),y=expression(S(t)))

# Mean survival time
1/lambda_emv

# Median survival time
log(2)/lambda_emv

# Quantiles and confidence intervals (this method could potentially
# produce negative estimates)
p = seq(.01,.99,by=.01)

tp = log(1/(1-p))/lambda_emv
se_tp = tp/sqrt(num_death)

tp_li = tp - se_tp*qnorm(.975)
tp_ui = tp + se_tp*qnorm(.975)

# For the median
tp_li[50]
tp_ui[50]

# Quantiles and confidence intervals V2 
tp_li_v2 = tp*exp(-qnorm(.975)/sqrt(num_death))
tp_ui_v2 = tp*exp(qnorm(.975)/sqrt(num_death))

# Median
tp_li_v2[50]
tp_ui_v2[50]

# Survival function and confidence intervals
St = exp(-lambda_emv*tp)
St_li = exp(-lambda_emv*tp_li_v2)
St_ui = exp(-lambda_emv*tp_ui_v2)

survs = data.frame(tp,St,St_li,St_ui)
ggplot(data=survs,aes(x=tp,y=St))+
  geom_line()+
  geom_line(aes(x=tp,y=St_li),col="lightblue4")+
  geom_line(aes(x=tp,y=St_ui),col="lightblue4")+
  labs(x="",y="")

# AIC
2-2*(num_death*log(lambda_emv)-(lambda_emv*sum_times))

# BIC
log(n)-2*(num_death*log(lambda_emv)-(lambda_emv*sum_times))
####################################################################

####################################################################
# Weibull model
logLikeWb = function(par,status,times){
  lambda = par[1]
  gamma = par[2]
  m = sum(status)
  term1 = m*log(gamma*lambda)
  term2 = (gamma-1)*sum(status*log(times))
  term3 = lambda*sum(times^gamma)
  return(term1+term2-term3)
}

# Optimisation using a quasi-Newton method
result = optim(par=c(1,1),fn=logLikeWb,method="L-BFGS-B",
              control=list(fnscale = -1,ndeps=c(.000001,.000001),
                           factr=.0000000001),
              times=lung$time,status=lung$status,
              lower=c(0.0000001,0.0000001),
              upper=c(100,100),hessian=T)

# Estimates
lambda=result$par[1];lambda
gamma=result$par[2];gamma

# Density 
ft_wb = lambda*gamma*(t^(gamma-1))*exp(-lambda*(t^gamma))

density_wb = data.frame(t,ft_wb)
ggplot(data=lung,aes(x=time,y=after_stat(density)))+
  geom_histogram(fill='lightblue3',col='black',
                 breaks=hist(lung$time,plot=F)$breaks)+
  geom_line(data=density_wb,aes(x=t,y=ft_wb))+  
  labs(x='Days',y='')

# Risk function
ht_wb = lambda*gamma*t^(gamma-1)
risk_wb = data.frame(t,ht_wb)
ggplot(data=risk_wb,aes(x=t,y=ht_wb))+
  geom_line()+
  labs(x=expression(t),y=expression(h(t)))

# Survival
St_wb = exp(-lambda*(t^gamma))
surv_wb = data.frame(t,St_wb)
ggplot(data=surv_wb,aes(x=t,y=St_wb))+
  geom_line()+
  labs(x=expression(t),y=expression(S(t)))

# Mean survival time
gamma((1/gamma)+1)/(lambda^(1/gamma))

# Quantiles and confidence intervals
covMat=solve(-result$hessian)
var_l=covMat[1,1]
var_g=covMat[2,2]
cov_lg = covMat[1,2]

tp = (log(1/(1-p))/lambda)^(1/gamma) 

cp=log(log(1/(1-p)))
se_tp = 1/(lambda*(gamma^2))*(gamma^2*var_l+
                    lambda^2*(cp-log(lambda))^2*var_g+
                    2*lambda*gamma*(cp-log(lambda))*cov_lg)^(1/2)

tp_li = exp(log(tp) - se_tp*qnorm(.975))
tp_ui = exp(log(tp) + se_tp*qnorm(.975))

# Median confidence intervals
tp[50]
tp_li[50]
tp_ui[50]

# Survival function and confidence intervals
St = exp(-lambda*(tp^gamma))
St_li = exp(-lambda*(tp_li^gamma))
St_ui = exp(-lambda*(tp_ui^gamma))

survs = data.frame(tp,St,St_li,St_ui)
ggplot(data=survs,aes(x=tp,y=St))+
  geom_line()+
  geom_line(aes(x=tp,y=St_li),col="lightblue4")+
  geom_line(aes(x=tp,y=St_ui),col="lightblue4")+
  labs(x="",y="")

# AIC
4-(2*(logLikeWb(c(lambda,gamma),lung$status,lung$time)))

# BIC
2*log(n)-(2*(logLikeWb(c(lambda,gamma),lung$status,lung$time)))
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
glimpse(lung)

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

# Model with just sex and ECOG as covariates
cox_lung = coxph(Surv(time,event=status)~sex+ph.ecog,
                 data=lung_cox2)
summary(cox_lung)


loglikrt = 2*(cox_lung$loglik[2]-cox_lung$loglik[1])

betas = cox_lung$coefficients
waldt = betas%*%solve(cox_lung$var)%*%betas
####################################################################
