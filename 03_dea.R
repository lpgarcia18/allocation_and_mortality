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
#Reading the impact predicitions
#########################################################################
pred <- read_csv("bases/impact_predictions.csv")

#########################################################################
#DEA 
#########################################################################
pred$HEALTH <- pred$PUBLIC_EXP_LAGGED2 * pred$PROP_PUBLIC_HEALTH_EXP_LAGGED2
pred$NON_HEALTH <- pred$PUBLIC_EXP_LAGGED2 -  pred$HEALTH
#Mantaining only treatments with 0.05 of signf.
pred_cut <- subset(pred, pred$IC_95_NEO_U5 > 0)
pred_cut <- subset(pred_cut, pred$IC_95_NEO > 0)


base_dea <- pred_cut

variaveis <- c("LOCATION", "INCOME_CLASS", "IMPACT_NEO", 
               "IMPACT_NEO_U5", "HEALTH", "NON_HEALTH") 

base_dea <- dplyr::select(base_dea, variaveis) 
base_dea$IMPACT_NEO <- base_dea$IMPACT_NEO 
base_dea$IMPACT_NEO_U5 <- base_dea$IMPACT_NEO_U5 
base_dea_poor <- subset(base_dea, base_dea$INCOME_CLASS %in% c("L", "LM"))
base_dea_rich <- subset(base_dea, base_dea$INCOME_CLASS %in% c("UM", "H"))
base_dea <- rbind(base_dea_poor, base_dea_rich)

base_dea1 <- read_data(datadea = base_dea,
                        dmus = 1,
                        ni = 2,
                        no = 2,
                        inputs = c(5:6) ,
                        outputs = c(3:4))
model_voo <- model_sbmeff(datadea = base_dea1,
                              dmu_ref = 1:nrow(base_dea), #All countries
                              dmu_eval =  1:nrow(base_dea_poor), #Poor countries (L, LM)
                              orientation = "oo",
                              rts = "vrs",
                              compute_target = TRUE,
                              returnlp = FALSE)


efi_voo <- efficiencies(model_voo)

efi_targ_helth <- targets(model_voo)[[1]][,1]
efi_targ_non_helth <- targets(model_voo)[[1]][,2]
efi_targ_impact_neo <- targets(model_voo)[[2]][,1] 
efi_targ_impact_neo_u5 <- targets(model_voo)[[2]][,2] 

efi <- data.frame(EFI = efi_voo,
                  TARG_HEALTH = efi_targ_helth, 
                  TARG_NON_HEALTH = efi_targ_non_helth,
                  TARG_IMPACT_NEO = efi_targ_impact_neo,
                  TARG_IMPACT_NEO_U5 = efi_targ_impact_neo_u5)
efi$TARG_PROP <- efi$TARG_HEALT/(efi$TARG_HEALT+efi$TARG_NON_HEALTH)
efi$LOCATION <- row.names(efi)

base_efi <- merge(pred, efi, by = "LOCATION")


#########################################################################
#Saving the databases
#########################################################################
write.csv(base_efi, "bases/efficiency_analysis.csv", row.names = F)

