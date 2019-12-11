# Table 1: Characteristic of Study Sample
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
pheno1 %>% tbl_summary(.,by = "BB",statistic = list(..continuous.. = "{mean} ({sd})",..categorical.. = "{n} / {N} ({p}%)")) 
pheno1 %>% tbl_summary(., statistic = list(..continuous.. = "{mean} ({sd})",..categorical.. = "{n} / {N} ({p}%)")) 
                  