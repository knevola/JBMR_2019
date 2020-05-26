rm(list = ls())
setwd("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/data")

library(xlsx)
library(meta)

data <- read.xlsx("Data_for_Meta_Analysis.xlsx", sheetIndex = 1)
data_cor <- read.csv("Data_for_Meta_Analysis_cor.csv")
data_cor_z <- data_cor[5:8,] 
data_cor <- data_cor[1:4,]
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

data_cor_19 <- data_cor[data_cor$miRNA == "miR-19a-3p",]
data_cor_186 <- data_cor[data_cor$miRNA == "miR-186-5p",]

miR19_cor_fixed<-metagen(TE = corr, seTE = Standard.Error,data= data_cor_19,comb.fixed = T, studlab = paste(Cohort))
miR186_cor_fixed <-metagen(TE = corr, seTE = Standard.Error,data= data_cor_186,comb.fixed = T, studlab = paste(Cohort))

forest(miR19_cor_fixed)
forest(miR186_cor_fixed)

# T-score
data_Tscore <- data[c(9,10,15,16),]
data_Tscore_miR19 <- data_Tscore[data_Tscore$miRNA == "miR-19a-3p",]
data_Tscore_miR186 <- data_Tscore[data_Tscore$miRNA == "miR-186-5p",]

miR19_fixed_Tscore<-metagen(TE = Estimate, seTE = SE,data= data_Tscore_miR19,comb.fixed = T, studlab = paste(Cohort, Measure))
miR186_fixed_Tscore <-metagen(TE = Estimate, seTE = SE,data= data_Tscore_miR186,comb.fixed = T, studlab = paste(Cohort, Measure))

forest(miR19_fixed_Tscore)
forest(miR186_fixed_Tscore)

# Z-score correlation
data_cor_z_19 <- data_cor_z[data_cor_z$miRNA == "miR-19a-3p",]
data_cor_z_186 <- data_cor_z[data_cor_z$miRNA == "miR-186-5p",]
miR19_fixed_Zscore<-metagen(TE = corr, seTE = Standard.Error,data= data_cor_z_19,comb.fixed = T, studlab = Cohort)
miR186_fixed_Zscore <-metagen(TE = corr, seTE = Standard.Error,data= data_cor_z_186,comb.fixed = T, studlab = Cohort)

forest(miR19_fixed_Zscore)
forest(miR186_fixed_Zscore)
