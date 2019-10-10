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

base_dea <- pred

variaveis <- c("LOCATION","PUBLIC_EXP_LAGGED2", "PROP_PUBLIC_HEALTH_EXP_LAGGED2", "PRED_NEO", 
               "PRED_NEO_U5", "INCOME_CLASS", "DELTA_RATE_NEO", "DELTA_RATE_NEO_U5", "HEALTH", "NON_HEALTH") 

base_dea <- dplyr::select(base_dea, variaveis) 
base_dea <- subset(base_dea, base_dea$PRED_NEO_U5 >= 0)
base_dea <- base_dea[,c(1,4,5,9,10)]
base_dea <- subset(base_dea, base_dea$LOCATION != "Madagascar")

base_dea1 <- read_data(datadea = base_dea,
                        dmus = 1,
                        ni = 2,
                        no = 2,
                        inputs = c(4:5) ,
                        outputs = c(2:3))
model_voo <- model_sbmeff(datadea = base_dea1,
                              dmu_ref = 1:104,
                              dmu_eval = 1:104,
                              orientation = "oo",
                              rts = "vrs",
                              compute_target = TRUE,
                              returnlp = FALSE)


efi_voo <- efficiencies(model_voo)

efi_targ_inp_helth <- targets(model_voo)[[1]][,1]
efi_targ_inp_non_helth <- targets(model_voo)[[1]][,2]
efi_targ_pred_neo <- targets(model_voo)[[2]][,1]
efi_targ_pred_neo_u5 <- targets(model_voo)[[2]][,2]

efi <- data.frame(EFI = efi_voo,
                  TARG_HEALTH = efi_targ_inp_helth, 
                  TARG_NON_HEALTH = efi_targ_inp_non_helth,
                  TARG_PRED_NEO = efi_targ_pred_neo,
                  TARG_PRED_NEO_U5 = efi_targ_pred_neo_u5)
efi$TARG_PROP <- efi$TARG_HEALT/(efi$TARG_HEALT+efi$TARG_NON_HEALTH)
efi$LOCATION <- row.names(efi)

base_efi <- merge(pred, efi, by = "LOCATION")

#########################################################################
#Saving the databases
#########################################################################
write.csv(base_efi, "bases/efficiency_analysis.csv", row.names = F)

