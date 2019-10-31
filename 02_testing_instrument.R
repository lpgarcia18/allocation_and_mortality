#########################################################################
#Settings
#########################################################################
options(scipen=999)
set.seed(233)

#########################################################################
#Libraries
#########################################################################
library(readr)
library(tidyverse)
library(reshape2)
library(matrixStats)
library(wbstats)
library(grf)
library(sf)
library(spData)
library(mice)
library(deaR)
library(gt)
library(AER)
#########################################################################
#Reading the completed database
#########################################################################
completed_base <- read_csv("bases/completed_base.csv")
names(completed_base)[1] <- "LOCATION"
#completed_base$PUBLIC_EXP_LAGGED2 <- round(completed_base$PUBLIC_EXP_LAGGED2,0) 
#completed_base$TAX_PPP_LAGGED3 <- round(completed_base$TAX_PPP_LAGGED3,0) 


#########################################################################
#Impact of Public Expenditure on Mortality
#########################################################################
# Impact on Neontatal Mortality -------------------------------------------

demography <- c("FERTILITY_RATE_LAGGED","URBAN_RATE_LAGGED", "POP_DENS", "PLUS_65_YEARS_LAGGED")
geography <- c("LONG", "LAT")
women_empowerment <- c("UNEMPLOYMENT_FEM_LAGGED", "SCHOOL_FEM_LAGGED",            
"WOMEN_PARLIAMENT_LAGGED")
schooling <- c("SCHOOL_LIFE_EXP_LAGGED", "OUT_OF_SCHOOL_LAGGED")
governance <- c("CONTROL_CORRUPTION_LAGGED", "GOV_EFFECTIVENESS_LAGGED", 
                "POLITICAL_STABILITY_LAGGED", "REGULATORY_QUALITY_LAGGED", "RULE_OF_LAW_LAGGED")
sanitation <- c("BASIC_SANITATION_LAGGED", "BASIC_WATER_LAGGED")
economy_income <- c("GDP_PPP_LAGGED1", "GINI_LAGGED","POVERTY_GAP_LAGGED",           
"INFLATION_LAGGED", "UNEMPLOYMENT_LAGGED", "OOP_PPP_LAGGED3")
health <- c("DOCTORS_LAGGED", "DELIVERY_ASSISTANCE_LAGGED", "AIDS_PREVALENCE_LAGGED",       
"MALARIA_INCIDENCE_LAGGED", "DPT_LAGGED", "HOSPITAL_BEDS_LAGGED")
nutrition <- c("UNDERNOURISHMENT_LAGGED")
energy <- c("ELECTRICITY_LAGGED")

external_factors <- c(governance, economy_income, demography, geography, sanitation, 
                      women_empowerment,schooling,  
                      nutrition, health, energy)


Y_neo <- completed_base$MEAN_RATE_NEO
Y_neo_u5 <- completed_base$MEAN_RATE_NEO_U5
X <- dplyr::select(completed_base, external_factors)
W <- completed_base$PUBLIC_EXP_LAGGED2
Z <- completed_base$TAX_PPP_LAGGED2


#Analysing the endogeneity
##Neonatal
X1 <- as.matrix(X[,c(1:32)])
lm(Y_neo ~ W+X1) %>% summary()#Linear regression
lm(W~Z+X1)%>% summary()#First stage 2sls
ivreg(Y_neo ~ W+X1 | X1 + Z) %>%
   summary(vcov = sandwich, diagnostics=TRUE)#2sls whith Wu-Hausman test

##28d>5u
X1 <- as.matrix(X[,c(1:32)])
lm(Y_neo_u5 ~ W+X1) %>% summary()#Linear regression
lm(W~Z+X1)%>% summary()#First stage 2sls
ivreg(Y_neo_u5 ~ W+X1 | X1 + Z) %>%
   summary(vcov = sandwich, diagnostics=TRUE)#2sls whith Wu-Hausman test
