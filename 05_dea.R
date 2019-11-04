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

#Mantaining only treatments negatives with 0.05 of signf.
#neo
pred_cut_neo <- subset(pred, pred$sig_neo == "yes")



variables <- c("LOCATION", "INCOME_CLASS", "IMPACT_NEO", 
               "HEALTH", "NON_HEALTH") 

base_neo <- pred_cut_neo

base_neo <- dplyr::select(base_neo, variables) 


base_dea_neo <- read_data(datadea = base_neo,
                        dmus = 1,
                        ni = 2,
                        no = 2,
                        inputs = c(4,5), 
                        outputs = c(3))
model_voo_neo <- model_sbmeff(datadea = base_dea_neo,
                              dmu_ref = 1:nrow(base_neo), 
                              dmu_eval =  1:nrow(base_neo), 
                              orientation = "oo",
                              rts = "vrs",
                              compute_target = TRUE,
                              returnlp = FALSE)


efi_voo_neo <- efficiencies(model_voo_neo)

efi_neo <- data.frame(EFI_NEO = efi_voo_neo)

efi_neo$LOCATION <- row.names(efi_neo)

#Mantaining only treatments negatives with 0.05 of signf.
#28d>5y
pred_cut_neo_u5 <- subset(pred, pred$sig_neo_u5 == "yes")



variables <- c("LOCATION", "INCOME_CLASS", "IMPACT_NEO_U5",
               "HEALTH", "NON_HEALTH") 

base_neo_u5 <- pred_cut_neo_u5

base_neo_u5 <- dplyr::select(base_neo_u5, variables) 


base_dea_neo_u5 <- read_data(datadea = base_neo_u5,
                        dmus = 1,
                        ni = 2,
                        no = 2,
                        inputs = c(4,5), 
                        outputs = c(3))
model_voo_neo_u5 <- model_sbmeff(datadea = base_dea_neo_u5,
                              dmu_ref = 1:nrow(base_neo_u5), 
                              dmu_eval =  1:nrow(base_neo_u5), 
                              orientation = "oo",
                              rts = "vrs",
                              compute_target = TRUE,
                              returnlp = FALSE)


efi_voo_neo_u5 <- efficiencies(model_voo_neo_u5)

efi_neo_u5 <- data.frame(EFI_NEO_U5 = efi_voo_neo_u5)

efi_neo_u5$LOCATION <- row.names(efi_neo_u5)

#Merging
base_efi <- merge(pred, efi_neo, by = "LOCATION", all = T)
base_efi <- merge(base_efi, efi_neo_u5, by = "LOCATION", all = T)


#########################################################################
#Saving the databases
#########################################################################
write.csv(base_efi, "bases/efficiency_analysis.csv", row.names = F)

