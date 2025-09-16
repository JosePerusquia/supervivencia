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
library(MASS)           # Version 7.3-60
library(flexsurv)       # Version 2.3.2
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
# Exponential survival regression
lung2 = lung%>%
  filter(ph.ecog%in%c("0","1","2"))%>%
  dplyr::select(time,status,sex,ph.ecog)%>%
  mutate(sex_m = if_else(sex == 2, 1, 0))%>%
  mutate(ecog2 = if_else(ph.ecog == 1, 1, 0))%>%
  mutate(ecog3 = if_else(ph.ecog == 2, 1, 0))%>%
  dplyr::select(time,status,sex_m,ecog2,ecog3)

logLikeExp = function(par,status,times,x,n){
  beta_sex = par[1]
  beta_ecog2=par[2]
  beta_ecog3=par[3]
  lambda = par[4]
  
  beta=c(beta_sex,beta_ecog2,beta_ecog3)
  m = sum(status)
  
  term1 = m*log(lambda)
  term2 = 0
  term3 = 0
  for(i in 1:n){
    term2 = term2 + status[i]*sum((beta*x[i,]))
    term3 = term3 + lambda*times[i]*exp(sum((beta*x[i,])))
  }
  
  return(term1+term2-term3)
}

# Optimisation using a quasi-Newton method
result = optim(par=c(0,0,0,1),fn=logLikeExp,method="L-BFGS-B",
               control=list(fnscale = -1,
                            ndeps=c(.000001,.000001,
                                    .000001,.000001),
                            factr=.0000000001),
               times=lung2$time,status=lung2$status,
               x = lung2[,c(3,4,5)], n =226,
               lower=c(-100,-100,-100,0.0000001),
               upper=c(100,100,100,100),hessian=T)

# Estimates
beta1=result$par[1];beta1
beta2=result$par[2];beta2
beta3=result$par[3];beta3
lambda = result$par[4];lambda

# Flexsurvreg estimates
lung = lung%>%
  filter(ph.ecog%in%c("0","1","2"))%>%
  mutate(sex=factor(sex),
         ph.ecog=factor(ph.ecog))
  
mod=flexsurvreg(Surv(time,event=status)~sex+ph.ecog,
        dist='exp',data=lung)
mod
####################################################################

####################################################################
# Cox - Snell residuals
res_cs = residuals(mod,type='coxsnell')

# KM estimate using Cox - Snell residuals
lung$CS = res_cs
fit_cs = survfit(Surv(CS, status) ~ 1, data = lung)

# Plot cumulative hazard
df = data.frame(x=fit_cs$time,y=-log(fit_cs$surv))

ggplot(data=df,aes(x=x,y=y))+
  geom_line(col='darkgrey')+
  geom_point(size=.5)+
  labs(x="Cox-Snell residuals",y=expression(hat(H[t])))+
  geom_abline(intercept=0,slope=1,col="red")
####################################################################

####################################################################
# Deviance residuals
res_mar = lung$status-res_cs
res_dev = sign(res_mar)*sqrt((-2*(res_mar+
                lung$status*log(lung$status-res_mar))))
lung$Dev = res_dev

ggplot(data=lung,aes(x=time,y=Dev))+
  geom_point()+
  geom_hline(yintercept=0,col='red')+
  labs(x="Survival time",y="Deviance residuals")
####################################################################
