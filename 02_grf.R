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

demography <- c("FERTILITY_RATE_LAGGED","URBAN_RATE_LAGGED", "POP_DENS", "pop")
geography <- c("LONG", "LAT")
women_empowerment <- c("UNEMPLOYMENT_FEM_LAGGED", "SCHOOL_FEM_LAGGED",            
"WOMEN_PARLIAMENT_LAGGED")
schooling <- c("SCHOOL_LIFE_EXP_LAGGED", "OUT_OF_SCHOOL_LAGGED")
governance <- c("CONTROL_CORRUPTION_LAGGED", "GOV_EFFECTIVENESS_LAGGED", 
                "POLITICAL_STABILITY_LAGGED", "REGULATORY_QUALITY_LAGGED", "RULE_OF_LAW_LAGGED")
sanitation <- c("BASIC_SANITATION_LAGGED", "BASIC_WATER_LAGGED")
economy_income <- c("GDP_PPP_LAGGED1", "GINI_LAGGED","POVERTY_GAP_LAGGED",           
"INFLATION_LAGGED", "UNEMPLOYMENT_LAGGED", "OOP_PPP_LAGGED4")
health <- c("DOCTORS_LAGGED", "DELIVERY_ASSISTANCE_LAGGED", "AIDS_PREVALENCE_LAGGED",       
"MALARIA_INCIDENCE_LAGGED", "DPT_LAGGED", "HOSPITAL_BEDS_LAGGED")
nutrition <- c("UNDERNOURISHMENT_LAGGED")
energy <- c("ELECTRICITY_LAGGED")

#"PLUS_65_YEARS_LAGGED",
external_factors <- c(governance, economy_income, demography, geography, sanitation, 
                      women_empowerment,schooling,  
                      nutrition, health, energy)


Y_neo <- completed_base$DELTA_RATE_NEO
Y_neo_u5 <- completed_base$DELTA_RATE_NEO_U5
X <- dplyr::select(completed_base, external_factors)
W <- completed_base$PUBLIC_EXP_LAGGED2
Z <- completed_base$PLUS_65_YEARS_LAGGED
Y_neo.forest <- regression_forest(X, Y_neo, seed = 233)
Y_neo.hat <- predict(Y_neo.forest)$predictions
Y_neo_u5.forest <- regression_forest(X, Y_neo_u5, seed = 233)
Y_neo_u5.hat <- predict(Y_neo_u5.forest)$predictions
W.forest <- regression_forest (X , W, seed = 233)
W.hat <- predict(W.forest)$predictions
Z.forest <- regression_forest(X, Z, seed = 233)
Z.hat <- predict(Z.forest)$predictions
#Effect neo
iv_neo_raw <- instrumental_forest(X, Y_neo, W, Z,
                    Y_neo.hat, W.hat, Z.hat,
                    mtry = sqrt(ncol(X)),
                    num.trees = 500000,
                    honesty = T,
                    compute.oob.predictions = TRUE,
                    seed = 233)
varimp_neo <- variable_importance(iv_neo_raw)
selected_idx_neo <- which(varimp_neo > mean(varimp_neo))
iv_neo <- instrumental_forest(X[,selected_idx_neo], Y_neo, W, Z,
                    Y_neo.hat, W.hat, Z.hat,
                    mtry = sqrt(ncol(X)),
                    num.trees = 500000,
                    honesty = T,
                    compute.oob.predictions = TRUE,
                    seed = 233)
pred_neo <- predict(iv_neo, newdata = NULL, estimate.variance = TRUE, set.seed(233))
pred_neo$PRED_NEO <- pred_neo$predictions 
pred_neo$DESV_PADR_NEO <- sqrt(pred_neo$variance.estimates)
pred_neo$SE_NEO <- pred_neo$DESV_PADR/sqrt(500000)
pred_neo$IC_95_NEO <- pred_neo$predictions+1.96*(pred_neo$SE)
pred_neo$IC_05_NEO <- pred_neo$predictions-1.96*(pred_neo$SE)




 
#Effect neo_u5
iv_neo_u5_raw <- instrumental_forest(X, Y_neo_u5, W, Z,
                    Y_neo_u5.hat, W.hat, Z.hat,
                    mtry = sqrt(ncol(X)),
                    num.trees = 500000,
                    honesty = T,
                    compute.oob.predictions = TRUE,
                    seed = 233)
varimp_neo_u5 <- variable_importance(iv_neo_u5_raw)
selected_idx_neo_u5 <- which(varimp_neo_u5 > mean(varimp_neo_u5))
iv_neo_u5 <- instrumental_forest(X[,selected_idx_neo_u5], Y_neo_u5, W, Z,
                    Y_neo_u5.hat, W.hat, Z.hat,
                    mtry = sqrt(ncol(X)),
                    num.trees = 500000,
                    honesty = T,
                    compute.oob.predictions = TRUE,
                    seed = 233)
pred_neo_u5 <- predict(iv_neo_u5, newdata = NULL, estimate.variance = TRUE, set.seed(233))
pred_neo_u5$PRED_NEO_U5 <- pred_neo_u5$predictions 
pred_neo_u5$DESV_PADR_NEO_U5 <- sqrt(pred_neo_u5$variance.estimates)
pred_neo_u5$SE_NEO_U5 <- pred_neo_u5$DESV_PADR/sqrt(500000)
pred_neo_u5$IC_95_NEO_U5 <- pred_neo_u5$predictions+1.96*(pred_neo_u5$SE)
pred_neo_u5$IC_05_NEO_U5 <- pred_neo_u5$predictions-1.96*(pred_neo_u5$SE)



pred <- cbind(completed_base, pred_neo[,c(5,9,8)], pred_neo_u5[,c(5,9,8)])
pred$IMPACT_NEO <- pred$PRED_NEO * pred$PUBLIC_EXP_LAGGED2
pred$IMPACT_NEO_U5 <- pred$PRED_NEO_U5 * pred$PUBLIC_EXP_LAGGED2


#########################################################################
#Saving the databases
#########################################################################
write.csv(pred, "bases/impact_predictions.csv", row.names = F)
