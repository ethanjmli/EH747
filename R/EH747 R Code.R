#####EH747 Code#####
##Ethan Li##
##Fall 2023##

######Packages#####
library(here)
library(dplyr)
library(readxl)
library(tidyverse)
library(lmtest)
library(epitools)
library(ggplot2)

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


######Directly Standardized Rate Ratio######
cancer_data_std <- cancer_data %>% group_by(age, calendar_time) %>% mutate (sum_PT = sum(P_years))
dir_std_rate_ratio <- cancer_data_std %>% group_by(blood_lead) %>% summarise(std_rate = sum(cancer_rate*sum_PT))
161.87134/76.99434

######Poisson Modeling#####
poisson_model <- glm(obs_cancer ~ blood_lead_factor + age_factor + calendar_time_factor, family = "poisson",offset=log(P_years),data=cancer_data)
summary(poisson_model)

# Call:
#   glm(formula = obs_cancer ~ blood_lead_factor + age_factor + calendar_time_factor, 
#       family = "poisson", data = cancer_data, offset = log(P_years))
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)               -9.440399   0.313010 -30.160  < 2e-16 ***
#   blood_lead_factorHigh      0.768173   0.195324   3.933 8.40e-05 ***
#   age_factor50-65            1.297686   0.310202   4.183 2.87e-05 ***
#   age_factor65+              2.002321   0.297662   6.727 1.73e-11 ***
#   calendar_time_factor+2000 -0.009742   0.187521  -0.052    0.959    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for poisson family taken to be 1)
# 
# Null deviance: 87.963  on 11  degrees of freedom
# Residual deviance: 11.893  on  7  degrees of freedom
# AIC: 67.977
# 
# Number of Fisher Scoring iterations: 5

cbind(exp(coef(poisson_model)),exp(confint(poisson_model)))

# 2.5 %       97.5 %
#   (Intercept)               0.0000794487 4.124004e-05 1.415111e-04
# blood_lead_factorHigh     2.1558233317 1.479702e+00 3.189601e+00
# age_factor50-65           3.6608143432 2.044078e+00 6.965519e+00
# age_factor65+             7.4062292361 4.261441e+00 1.381325e+01
# calendar_time_factor+2000 0.9903050767 6.867210e-01 1.435100e+00


quasipoisson_model <- glm(obs_cancer ~ blood_lead_factor + age_factor + calendar_time_factor, family = "quasipoisson",offset=log(P_years),data=cancer_data)
summary(quasipoisson_model)

# Call:
#   glm(formula = obs_cancer ~ blood_lead_factor + age_factor + calendar_time_factor, 
#       family = "quasipoisson", data = cancer_data, offset = log(P_years))
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)               -9.440399   0.389525 -24.236 5.18e-08 ***
#   blood_lead_factorHigh      0.768173   0.243071   3.160   0.0159 *  
#   age_factor50-65            1.297686   0.386030   3.362   0.0121 *  
#   age_factor65+              2.002321   0.370425   5.405   0.0010 ** 
#   calendar_time_factor+2000 -0.009742   0.233360  -0.042   0.9679    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for quasipoisson family taken to be 1.548651)
# 
# Null deviance: 87.963  on 11  degrees of freedom
# Residual deviance: 11.893  on  7  degrees of freedom
# AIC: NA
# 
# Number of Fisher Scoring iterations: 5

cbind(exp(coef(quasipoisson_model)),exp(confint(quasipoisson_model)))

# 2.5 %       97.5 %
#   (Intercept)               0.0000794487 3.463638e-05 1.613066e-04
# blood_lead_factorHigh     2.1558233317 1.351743e+00 3.521407e+00
# age_factor50-65           3.6608143432 1.783341e+00 8.257240e+00
# age_factor65+             7.4062292361 3.752305e+00 1.631528e+01
# calendar_time_factor+2000 0.9903050767 6.279787e-01 1.573593e+00


####Case 2 ####
cadmium.data <- read_xlsx(here::here("Data","cadmium2.xlsx"), na = ".")
summary(cadmium.data)

#High B2PGC > 486 ug/g creatinine
#High creatinine > 1.4 mg/100 dL
#Use ln to transform

######Clean Data, Descriptive stats######
cadmium.data <- cadmium.data %>% mutate(b2pgc = as.numeric(b2pgc), scre = as.numeric(scre))
cadmium.data <- cadmium.data %>% mutate(high_B2PGC = ifelse(b2pgc > 486, 1, 0), 
                                        high_scre = ifelse(scre >= 1.4, 1, 0))
cadmium.data <- cadmium.data %>% mutate(high_B2PGC = factor(high_B2PGC, levels = c(0,1),labels = c("low","high")),
                                        high_scre = factor(high_scre, levels = c(0,1),labels = c("low","high")))
cadmium.data <- cadmium.data %>% mutate(exposure_factor = factor(exposure, levels= c(0,1), labels = c("unexposed","exposed")))
cadmium.data <- cadmium.data %>% mutate(ln_B2PGC = log(b2pgc),ln_scre = log(scre))



b2pgc.table <- table(cadmium.data$exposure_factor,cadmium.data$high_B2PGC,useNA = "ifany")
addmargins(b2pgc.table)
#           low high <NA> Sum
# unexposed  28    1    3  32
# exposed    27   18    0  45
# Sum        55   19    3  77

cbind(b2pgc.table,
      prop.table(table(cadmium.data$exposure_factor, cadmium.data$high_B2PGC),margin=1))
#           low high <NA>       low       high
# unexposed  28    1    3 0.9655172 0.03448276
# exposed    27   18    0 0.6000000 0.40000000


scre.table <- table(cadmium.data$exposure_factor,cadmium.data$high_scre,useNA = "ifany")
addmargins(scre.table)
#           low high <NA> Sum
# unexposed  30    1    1  32
# exposed    33   12    0  45
# Sum        63   13    1  77
cbind(scre.table,
      prop.table(table(cadmium.data$exposure_factor, cadmium.data$high_scre),margin=1))
#           low high <NA>       low       high
# unexposed  30    1    1 0.9677419 0.03225806
# exposed    33   12    0 0.7333333 0.26666667


######T test#####
#T test B2PGC exposed vs unexposed
t.test(ln_B2PGC ~ exposure_factor, data = cadmium.data)

# Welch Two Sample t-test
# 
# data:  ln_B2PGC by exposure_factor
# t = -4.1795, df = 53.055, p-value = 0.0001097
# alternative hypothesis: true difference in means between group unexposed and group exposed is not equal to 0
# 95 percent confidence interval:
#   -1.7352496 -0.6098502
# sample estimates:
#   mean in group unexposed   mean in group exposed 
# 5.243112                6.415662


#T test SCRE exposed vs unexposed
t.test(ln_scre ~ exposure_factor, data = cadmium.data)

# Welch Two Sample t-test
# 
# data:  ln_scre by exposure_factor
# t = -2.8135, df = 72.967, p-value = 0.006293
# alternative hypothesis: true difference in means between group unexposed and group exposed is not equal to 0
# 95 percent confidence interval:
#   -0.23199656 -0.03960338
# sample estimates:
#   mean in group unexposed   mean in group exposed 
# 0.01345381              0.14925378 


######Chi Square Test#####
#Chi Square Test BPGC
chisq.test(cadmium.data$high_B2PGC,cadmium.data$exposure_factor,correct=T)

# Pearson's Chi-squared test with Yates' continuity correction
# 
# data:  cadmium.data$high_B2PGC and cadmium.data$exposure_factor
# X-squared = 10.505, df = 1, p-value = 0.00119

#Chi Square Test SCRE
chisq.test(cadmium.data$high_scre,cadmium.data$exposure_factor,correct=T)

# Pearson's Chi-squared test with Yates' continuity correction
# 
# data:  cadmium.data$high_scre and cadmium.data$exposure_factor
# X-squared = 5.5559, df = 1, p-value = 0.01842


######Linear Regression#####
lm_BPGC <- glm(ln_B2PGC ~ dose + age, family = "gaussian",data=cadmium.data)
lm_BPGC
summary(lm_BPGC)
coef(lm_BPGC)
confint(lm_BPGC)
cbind(coef(lm_BPGC),confint(lm_BPGC))

#                                  2.5 %      97.5 %
#(Intercept)  4.8550563953  3.8699553828 5.840157408
# dose        0.0008986238  0.0006400692 0.001157178
# age         0.0102514505 -0.0093055217 0.029808423

lm_scre_withsbp <- glm(ln_scre ~ dose + sbp, family = "gaussian",data=cadmium.data)
lm_scre_withoutsbp <- glm(ln_scre ~ dose, family = "gaussian", data = cadmium.data)
lm_scre_EMMsbp <- glm(ln_scre ~ dose + sbp + dose*sbp, family = "gaussian",data = cadmium.data)

summary(lm_scre_withsbp)
coef(lm_scre_withsbp)
cbind(coef(lm_scre_withsbp),confint(lm_scre_withsbp))
#                                   2.5 %       97.5 %
#(Intercept)  -1.719738e-01 -5.192498e-01 0.1753021581
# dose         9.499613e-05  5.163189e-05 0.0001383604
# sbp          1.607016e-03 -1.095667e-03 0.0043096995
coef(lm_scre_withsbp)[2]

summary(lm_scre_withoutsbp)
coef(lm_scre_withoutsbp)
cbind(coef(lm_scre_withoutsbp),confint(lm_scre_withoutsbp))
#                                 2.5 %       97.5 %
#(Intercept)  0.032089224 -2.113088e-02 0.0853093261
# dose        0.000100951  5.865703e-05 0.0001432449
coef(lm_scre_withoutsbp)[2]

((coef(lm_scre_withsbp)[2] - coef(lm_scre_withoutsbp)[2])/coef(lm_scre_withsbp)[2])*100
lmtest::lrtest(lm_scre_withoutsbp,lm_scre_withsbp)

summary(lm_scre_EMMsbp)
cbind(coef(lm_scre_EMMsbp),confint(lm_scre_EMMsbp))
lmtest::lrtest(lm_scre_withoutsbp,lm_scre_EMMsbp)


######Odds Ratio#####
cadmium.data.noNA <- cadmium.data %>% drop_na()
oddsratio(cadmium.data.noNA$exposure_factor,cadmium.data.noNA$high_B2PGC)



######Logistic Regression High B2PGC vs Dose#####
logm_BPGC_Dose <- glm(high_B2PGC ~ dose, family = "binomial",data = cadmium.data)
summary(logm_BPGC_Dose)
cbind(coef(logm_BPGC_Dose),confint(logm_BPGC_Dose))
exp(coef(logm_BPGC_Dose)[2])

logm_BPGC_Dose_withage <- glm(high_B2PGC ~ dose + age, family = "binomial",data = cadmium.data)
summary(logm_BPGC_Dose_withage)
cbind(coef(logm_BPGC_Dose_withage),confint(logm_BPGC_Dose_withage))
exp(coef(logm_BPGC_Dose_withage)[2])

lmtest::lrtest(logm_BPGC_Dose,logm_BPGC_Dose_withage)
((coef(logm_BPGC_Dose)[2]-coef(logm_BPGC_Dose_withage)[2])/coef(logm_BPGC_Dose_withage)[2])*100

cadmium.data$Probability <- exp(logm_BPGC_Dose$coefficients[1]+cadmium.data$dose*logm_BPGC_Dose$coefficients[2])/(1+exp(logm_BPGC_Dose$coefficients[1]+cadmium.data$dose*logm_BPGC_Dose$coefficients[2]))

cadmium.data %>% 
  ggplot(aes(x=dose,y=Probability))+
    geom_point()+
    geom_smooth(method = "glm",method.args = list(family = "binomial"),se=FALSE)+
    labs(x = "Cadmium Dose (mg/m3-d)",y="Probability of Abnormal B2PGC")+
  theme_bw()+
  theme(
    plot.background = element_blank(),
    #panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  ) 
  
######Get mean dose for exposed vs unexposed#####
cadmium.data.meandose <- cadmium.data %>% group_by(exposure_factor) %>% summarize(meandose = mean(dose,na.rm=T))


####Case 3####
asthma.data <- read.csv(here::here("Data","asthma.csv"))
summary(asthma.data)

asthma.data <- asthma.data %>%
  mutate(Asthma_factor = factor(Asthma,levels = c(0,1),labels = c("No Asthma","Asthma")),
         Gender_factor = factor(Gender, levels= c(0,1)),
         Parental_Asthma_History_factor = factor(Parental_Asthma_History, levels = c(0,1),labels = c("No PAH","Yes PAH")),
         NO2_Quartile = cut(NO2, 
                            breaks = quantile(NO2,probs = c(0,0.25,0.5,0.75,1)), 
                            labels = c("First","Second","Third","Fourth"),
                            right=FALSE,
                            include.lowest=TRUE),
         Income_Quartile = cut(Income, 
                               breaks = quantile(Income,probs = c(0,0.25,0.5,0.75,1)), 
                               labels = c("First","Second","Third","Fourth"),
                               right=FALSE,
                               include.lowest=TRUE)
                           
         ) %>%
  mutate(NO2_first_quartile = ifelse(NO2_Quartile == "First",1,0),
         NO2_second_quartile = ifelse(NO2_Quartile == "Second",1,0),
         NO2_third_quartile = ifelse(NO2_Quartile == "Third",1,0),
         NO2_fourth_quartile = ifelse(NO2_Quartile == "Fourth",1,0),
         Income_first_quartile = ifelse(Income_Quartile == "First",1,0),
         Income_second_quartile = ifelse(Income_Quartile == "Second",1,0),
         Income_third_quartile = ifelse(Income_Quartile == "Third",1,0),
         Income_fourth_quartile = ifelse(Income_Quartile == "Fourth",1,0))
summary(asthma.data)

logistic_model_NO2cont <- glm(Asthma_factor ~ NO2 + Gender_factor + Parental_Asthma_History_factor + Income_second_quartile + Income_third_quartile + Income_fourth_quartile, 
                              family = binomial(link = "logit"),data = asthma.data)
logistic_model_NO2quant <- glm(Asthma_factor ~ NO2_second_quartile + NO2_third_quartile + NO2_fourth_quartile +  Gender_factor + Parental_Asthma_History_factor + Income_second_quartile + Income_third_quartile + Income_fourth_quartile, 
                               family = binomial(link = "logit"),data = asthma.data)

cbind(coef(logistic_model_NO2cont),confint(logistic_model_NO2cont))
#Income = continuous                                          2.5 %        97.5 %
# (Intercept)                           -8.947303e+00 -1.093491e+01 -7.115745e+00
# NO2                                    1.468365e-01  1.022815e-01  1.941701e-01
# Gender_factor1                         9.789114e-01  4.604003e-01  1.518549e+00
# Parental_Asthma_History_factorYes PAH  3.470165e+00  2.946150e+00  4.033142e+00
# Income                                 2.965801e-06 -8.842975e-06  1.449948e-05

#Income = categorical nominal                             2.5 %      97.5 %
# (Intercept)                           -8.7874515 -10.7823772 -6.95150978
# NO2                                    0.1525150   0.1066613  0.20133463
# Gender_factor1                         1.0009771   0.4788751  1.54469254
# Parental_Asthma_History_factorYes PAH  3.4925582   2.9644988  4.06026058
# Income_second_quartile                -0.1309688  -0.8542194  0.58914543
# Income_third_quartile                 -0.7099003  -1.4588138  0.02224836
# Income_fourth_quartile                -0.2196052  -0.9188928  0.47628701

cbind(coef(logistic_model_NO2quant),confint(logistic_model_NO2quant))
#Income = continuous                                          2.5 %        97.5 %
# (Intercept)                           -4.837983e+00 -5.858676e+00 -3.911337e+00
# NO2_second_quartile                    2.570119e-01 -5.993841e-01  1.126454e+00
# NO2_third_quartile                     1.450914e+00  6.642820e-01  2.284405e+00
# NO2_fourth_quartile                    2.076913e+00  1.337004e+00  2.875240e+00
# Gender_factor1                         9.894509e-01  4.697918e-01  1.530819e+00
# Parental_Asthma_History_factorYes PAH  3.593442e+00  3.061127e+00  4.167235e+00
# Income                                 2.297795e-06 -9.693807e-06  1.398346e-05

#Income = categoricla nominal                           2.5 %      97.5 %
# (Intercept)                           -4.5423028 -5.5410845 -3.64275036
# NO2_second_quartile                    0.2908322 -0.5674080  1.16332269
# NO2_third_quartile                     1.4406663  0.6474242  2.28127851
# NO2_fourth_quartile                    2.1286082  1.3805245  2.93691034
# Gender_factor1                         1.0051872  0.4825778  1.54975522
# Parental_Asthma_History_factorYes PAH  3.6099956  3.0744740  4.18755091
# Income_second_quartile                -0.1186325 -0.8414738  0.60083881
# Income_third_quartile                 -0.6345041 -1.3775093  0.09319591
# Income_fourth_quartile                -0.2054511 -0.9069424  0.49232693

cbind(exp(coef(logistic_model_NO2cont)),exp(confint(logistic_model_NO2cont)))
IQR(asthma.data$NO2)
coef(logistic_model_NO2cont)[2]*IQR(asthma.data$NO2)
exp(coef(logistic_model_NO2cont)[2]*IQR(asthma.data$NO2))


366+365*4
#1826 

asthma.data.losttofu <- asthma.data %>%
  mutate(lost_to_fu = ifelse(daysdiff < 1826, 1,0)) %>%
  group_by(NO2_Quartile) %>%
  summarise(lost_to_fu_percent = mean(lost_to_fu)*100)

head(asthma.data.losttofu)
#   NO2_Quartile lost_to_fu_percent
#   <fct>                     <dbl>
# 1 First                      37.5
# 2 Second                     34  
# 3 Third                      38  
# 4 Fourth                     35