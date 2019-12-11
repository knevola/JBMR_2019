# Table 2: Estimated BMD after Adjustment for Covariates
# Note: The output of this script is not Table 2 in its final form
# However all information needed to calculate the information in table 2 is available here
rm(list = ls())
library(tidyverse)
library(gt)
library(dplyr)
library(gtsummary)
library(ggplot2)
library(emmeans)
library(reshape2)
library(matrixStats)
library(haven)
# Read in data
setwd('/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/data')
miRNAdat <- read.csv('l_mrna_2011_m_0797s_17_c1.csv')
pheno <- read.csv('PhenoData_5_28.csv')
pheno1 <- pheno
pheno_female <- pheno1 %>% filter(., SEX == "Female")
pheno_male <- pheno1 %>%  filter(., SEX == "Male")

fnmodel <- lm(f8cbnbmd~BB + AGE8 + SEX + HGT8 + WGT8, data = pheno1)
ftomodel<-lm(f8cbtobmd~BB + AGE8 + SEX + HGT8 + WGT8, data = pheno1)
ftrmodel <- lm(f8cbtrbmd~BB + AGE8 + SEX + HGT8 + WGT8, data = pheno1)
s2model <- lm(s8cbl2bd~BB + AGE8 + SEX + HGT8 + WGT8, data = pheno1)
s3model <- lm(s8cbl3bd~BB + AGE8 + SEX + HGT8 + WGT8, data = pheno1)
s4model <- lm(s8cbl4bd~BB + AGE8 + SEX + HGT8 + WGT8, data = pheno1)
s24model <- lm(s8cbl24bd~BB + AGE8 + SEX + HGT8 + WGT8, data = pheno1)

summary(fnmodel)
summary(ftomodel)
summary(ftrmodel)
summary(s2model)
summary(s3model)
summary(s4model)
summary(s24model)

emmeans(fnmodel, specs = "BB")
fn <- emmeans(fnmodel, specs = "BB")
emmeans(ftomodel, specs = "BB")
fto<-emmeans(ftomodel, specs = "BB")
emmeans(ftrmodel, specs = "BB")
ftr<-emmeans(ftrmodel, specs = "BB")
s24<-emmeans(s24model, specs = "BB")
emmeans(s24model, specs = "BB")

fn <- data.frame(fn, BMD = "Femoral Neck")
fto <- data.frame(fto, BMD = "Total Femur")
ftr <- data.frame(ftr, BMD = "Femoral Trochanter")
s24 <- data.frame(s24, BMD = "Spine L2-L4")

bmd <- rbind(fn,fto,ftr,s24)

write.csv(bmd, "BMDEMMeans.csv")