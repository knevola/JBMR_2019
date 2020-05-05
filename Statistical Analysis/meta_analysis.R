rm(list = ls())
setwd("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/data")

library(xlsx)
library(meta)

data <- read.xlsx("Data_for_Meta_Analysis.xlsx", sheetIndex = 1)
data_miR19 <- data[data$miRNA == "miR-19a-3p",]
data_miR186 <- data[data$miRNA == "miR-186-5p",]

miR19_fixed<-metagen(TE = Estimate, seTE = SE,data= data_miR19,comb.fixed = T, studlab = paste(Cohort, Measure))
miR186_fixed <-metagen(TE = Estimate, seTE = SE,data= data_miR186,comb.fixed = T, studlab = paste(Cohort, Measure))

forest(miR19_fixed)
forest(miR186_fixed)

data_noEBMD <- data[-c(5:6,11:14,17:18),]
data_noEBMD_miR19 <- data_noEBMD[data_noEBMD$miRNA == "miR-19a-3p",]
data_noEBMD_miR186 <- data_noEBMD[data_noEBMD$miRNA == "miR-186-5p",]

miR19_fixed_no_eBMD<-metagen(TE = Estimate, seTE = SE,data= data_noEBMD_miR19,comb.fixed = T, studlab = paste(Cohort, Measure))
miR186_fixed_no_eBMD <-metagen(TE = Estimate, seTE = SE,data= data_noEBMD_miR186,comb.fixed = T, studlab = paste(Cohort, Measure))

forest(miR19_fixed_no_eBMD)
forest(miR186_fixed_no_eBMD)
