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
completed_base <- subset(completed_base, completed_base$INCOME_CLASS %in% c("L", "LM"))


#########################################################################
#Impact of Public Expenditure on Mortality
#########################################################################
# Impact on Neontatal Mortality -------------------------------------------

demography <- c("FERTILITY_RATE_LAGGED","URBAN_RATE_LAGGED", "POP_DENS", "PLUS_65_YEARS_LAGGED", "pop")
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


#Causal_forest
Y_neo.forest <- regression_forest(X, Y_neo, seed = 233)
Y_neo.hat <- predict(Y_neo.forest)$predictions
Y_neo_u5.forest <- regression_forest(X, Y_neo_u5, seed = 233)
Y_neo_u5.hat <- predict(Y_neo_u5.forest)$predictions
W.forest <- regression_forest (X , W, seed = 233)
W.hat <- predict(W.forest)$predictions
clust <- as.numeric(as.factor(completed_base$INCOME_CLASS))

#Effect neo
cs_neo_raw <- causal_forest(X, Y_neo, W, 
                    Y_neo.hat, W.hat,
                    #clusters = clust,
                    tune.parameters = T,
                    honesty = T,
                    num.trees = 200000,   
                    compute.oob.predictions = TRUE,
                    seed = 233)
varimp_neo <- variable_importance(cs_neo_raw)
selected_idx_neo <- which(varimp_neo > mean(varimp_neo))
cs_neo <- causal_forest(X[,selected_idx_neo], Y_neo, W, 
                    Y_neo.hat, W.hat, 
                    #clusters = clust,
                    tune.parameters = T,
                    honesty = T,
                    num.trees = 200000, 
                    compute.oob.predictions = TRUE,
                    seed = 233)
pred_neo <- predict(cs_neo, newdata = NULL, estimate.variance = TRUE, set.seed(233))
pred_neo$PRED_NEO <- pred_neo$predictions 
pred_neo$IC_95_NEO <- pred_neo$predictions+(1.96*(sqrt(pred_neo$variance.estimates))) #"The square-root of column ’variance.estimates’ is the standard error of CATE" - https://cran.r-project.org/web/packages/grf/grf.pdf
pred_neo$IC_05_NEO <- pred_neo$predictions-(1.96*(sqrt(pred_neo$variance.estimates))) #"The square-root of column ’variance.estimates’ is the standard error of CATE" - https://cran.r-project.org/web/packages/grf/grf.pdf


#Effect neo_u5
cs_neo_u5_raw <- causal_forest(X, Y_neo_u5, W, 
                    Y_neo_u5.hat, W.hat, 
                    #clusters = clust,
                    tune.parameters = T,
                    honesty = T,
                    num.trees = 200000, 
                    compute.oob.predictions = TRUE,
                    seed = 233)
varimp_neo_u5 <- variable_importance(cs_neo_u5_raw)
selected_idx_neo_u5 <- which(varimp_neo_u5 > mean(varimp_neo_u5))
cs_neo_u5 <- causal_forest(X[,selected_idx_neo_u5], Y_neo_u5, W, 
                    Y_neo_u5.hat, W.hat,
                    #clusters = clust,
                    tune.parameters = T,
                    honesty = T,
                    num.trees = 200000, 
                    compute.oob.predictions = TRUE,
                    seed = 233)
pred_neo_u5 <- predict(cs_neo_u5, newdata = NULL, estimate.variance = TRUE, set.seed(233))
pred_neo_u5$PRED_NEO_U5 <- pred_neo_u5$predictions 
pred_neo_u5$IC_95_NEO_U5 <- pred_neo_u5$predictions+(1.96*(sqrt(pred_neo_u5$variance.estimates)))#"The square-root of column ’variance.estimates’ is the standard error of CATE" - https://cran.r-project.org/web/packages/grf/grf.pdf
pred_neo_u5$IC_05_NEO_U5 <- pred_neo_u5$predictions-(1.96*(sqrt(pred_neo_u5$variance.estimates)))#"The square-root of column ’variance.estimates’ is the standard error of CATE" - https://cran.r-project.org/web/packages/grf/grf.pdf

pred <- cbind(completed_base, pred_neo[,c(5,7,6)], pred_neo_u5[,c(5,7,6)])

pred$pol_ic05neo <- NA
for(i in seq_along(pred$IC_05_NEO)){
      pred$pol_ic05neo[i] <- ifelse(pred$IC_05_NEO[i] > 0, "pos", "neg")
}
pred$pol_ic95neo <- NA
for(i in seq_along(pred$IC_95_NEO)){
      pred$pol_ic95neo[i] <- ifelse(pred$IC_95_NEO[i] > 0, "pos", "neg")
}
pred$sig_neo <- NA
for(i in seq_along(pred$IC_95_NEO)){
      pred$sig_neo[i] <- ifelse(pred$pol_ic05neo[i] == pred$pol_ic95neo[i], "yes", "no")
}

pred$pol_ic05neo_u5 <- NA
for(i in seq_along(pred$IC_05_NEO_U5)){
      pred$pol_ic05neo_u5[i] <- ifelse(pred$IC_05_NEO_U5[i] > 0, "pos", "neg")
}
pred$pol_ic95neo_u5 <- NA
for(i in seq_along(pred$IC_95_NEO_U5)){
      pred$pol_ic95neo_u5[i] <- ifelse(pred$IC_95_NEO_U5[i] > 0, "pos", "neg")
}
pred$sig_neo_u5 <- NA
for(i in seq_along(pred$IC_95_NEO_U5)){
      pred$sig_neo_u5[i] <- ifelse(pred$pol_ic05neo_u5[i] == pred$pol_ic95neo_u5[i], "yes", "no")
}

pred$IMPACT_NEO <- NA 
for(i in seq_along(pred$PRED_NEO)){
   pred$IMPACT_NEO[i] <- ifelse(pred$sig_neo[i] == "yes", -1 * pred$PRED_NEO[i] * pred$PUBLIC_EXP_LAGGED2[i], NA)
}

pred$IMPACT_NEO_U5 <- NA
for(i in seq_along(pred$PRED_NEO_U5)){
   pred$IMPACT_NEO_U5[i] <- ifelse(pred$sig_neo_u5[i] == "yes", -1 * pred$PRED_NEO_U5[i] * pred$PUBLIC_EXP_LAGGED2[i], NA)
}

ape_neo <- average_partial_effect(cs_neo, subset = NULL)
ape_neo_ci95 <- as.numeric(ape_neo[1]) + 1.96 * as.numeric(ape_neo[2])
ape_neo_ci05 <- as.numeric(ape_neo[1]) - 1.96 * as.numeric(ape_neo[2])

paste0(round(ape_neo[1],6), " (CI95%", round(ape_neo_ci95,6)," ", round(ape_neo_ci05,6), ")")

ape_neo_u5 <- average_partial_effect(cs_neo_u5, subset = NULL)
ape_neo_u5_ci95 <- as.numeric(ape_neo_u5[1]) + 1.96 * as.numeric(ape_neo_u5[2])
ape_neo_u5_ci05 <- as.numeric(ape_neo_u5[1]) - 1.96 * as.numeric(ape_neo_u5[2])

paste0(round(ape_neo_u5[1],6), " (CI95%", round(ape_neo_u5_ci95,6)," ", round(ape_neo_u5_ci05,6), ")")



#########################################################################
#Saving the databases
#########################################################################
write.csv(pred, "bases/impact_predictions.csv", row.names = F)




