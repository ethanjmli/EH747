#####EH747 Code#####
##Ethan Li##
##Fall 2023##

######Packages#####
library(here)
library(dplyr)
library(readxl)
library(tidyverse)


####Case 1####
######Create Data#####
age <- c("<50","50-65","65+","<50","50-65","65+","<50","50-65","65+","<50","50-65","65+")
calendar_time <- c(rep("<2000",6),rep("+2000",6))
obs_cancer <- c(2,5,10,1,14,20,3,7,13,8,15,19)
P_years <- c(30000,23000,12000,32500,20000,14000,23000,26000,25000,25000,23000,18000)
blood_lead <- c("<30 ug/ml","<30 ug/ml","<30 ug/ml","30+ ug/ml","30+ ug/ml","30+ ug/ml","<30 ug/ml","<30 ug/ml","<30 ug/ml","30+ ug/ml","30+ ug/ml","30+ ug/ml")
cancer_data <- as.data.frame(cbind(age,calendar_time,obs_cancer,P_years,blood_lead))
cancer_data$P_years <- as.numeric(cancer_data$P_years)
cancer_data$obs_cancer <- as.numeric(cancer_data$obs_cancer)
cancer_data <- cancer_data %>% mutate(blood_lead_factor = factor(blood_lead,levels = c("<30 ug/ml","30+ ug/ml"),labels = c("Low","High")),
                                      cancer_rate = obs_cancer/P_years,
                                      age_factor = factor(age,levels = c("<50","50-65","65+")),
                                      calendar_time_factor = factor(calendar_time,levels = c("<2000","+2000")))
levels(cancer_data$blood_lead_factor)


######Poisson Modeling#####
poisson_model <- glm(obs_cancer ~ blood_lead_factor + age_factor + calendar_time_factor, family = "poisson",offset=log(P_years),data=cancer_data)
summary(poisson_model)
cbind(exp(coef(poisson_model)),exp(confint(poisson_model)))


quasipoisson_model <- glm(obs_cancer ~ blood_lead_factor + age_factor + calendar_time_factor, family = "quasipoisson",offset=log(P_years),data=cancer_data)
summary(quasipoisson_model)
cbind(exp(coef(quasipoisson_model)),exp(confint(quasipoisson_model)))



