####################################################################
# Time dependent variable using time transform               
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
# Transform variables into factors
pbc_cox = pbc%>%
  dplyr::select(id,time,status,age,sex,edema,bili,albumin)%>%
  mutate(sex=factor(sex),
         edema=factor(edema))

# Cox model without time depenent variables
cox_pbc = coxph(Surv(time,status==2)~sex+edema+
                age+bili+albumin,
                data=pbc_cox)
summary(cox_pbc)
AIC(cox_pbc)

# Cox model with bili as a function of time, cannot easily
# retrieve the residuals
cox_pbc_tt = coxph(Surv(time,status==2)~sex+edema+
                  age+tt(bili)+albumin,
                  data=pbc_cox,
                  tt = function(x, t, ...) x * log(t))
summary(cox_pbc_tt)
AIC(cox_pbc_tt)

# Plot of the coefficient of bili
times = seq(1, max(pbc$time), by=.5)
beta_bili_tt = 0.01931264 * log(times)
df_plot = data.frame(times,beta_bili_tt)

ggplot(df_plot, aes(x = times, y = beta_bili_tt)) +
  geom_line(color = "steelblue", size = 1) +
  labs(x = "Time", y = expression(beta(t)))+
  theme_minimal()

# Manually introduce cut points and the variable as a function
# of time to easily retrieve the residuals
pbc_cox2 = survSplit(Surv(time, status==2) ~ .,data=pbc_cox, 
                     cut=unique(pbc_cox$time))

pbc_cox2$bili_tt = pbc_cox2$bili * log(pbc_cox2$time)

cox_pbc_tt2 = coxph(Surv(tstart, time, event) ~ age + sex + 
                    edema + bili_tt + albumin,
                    data = pbc_cox2)
summary(cox_pbc_tt2)
AIC(cox_pbc_tt2)

# Martingale residuals
res_mart = residuals(cox_pbc_tt2, type="martingale", 
                     collapse=pbc_cox2$id)

# Deviance residuals
res_dev = residuals(cox_pbc_tt2, type="deviance", 
                    collapse=pbc_cox2$id)

# delta-beta residuals
res_db = residuals(cox_pbc_tt2, type="dfbeta",
                   collapse=pbc_cox2$id)
#######################################################################################################
