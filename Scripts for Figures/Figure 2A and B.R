# Figure 2 A and B (Figure 2C,D, and E were all made in Prism)
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
mir_tech = "mirna_tech_17.sas7bdat"
mirnatech <- read_sas(mir_tech)

miRNA_delta_cq <- miRNAdat[-1]
miRNA_delta_cq <- -(miRNA_delta_cq-27)
miRNAdat <- cbind(miRNAdat[1], miRNA_delta_cq)

quantcon<-quantile(mirnatech$concentration, probs = seq(0,1, 0.1))
mirnatech$rankcon <- cut(mirnatech$concentration, quantcon)
quantqual <-quantile(mirnatech$RNA_quality, probs = seq(0,1, 0.1), na.rm = T)
mirnatech$rankqual <- cut(mirnatech$RNA_quality, quantqual)
quant260 <- quantile(mirnatech$`_260_280`, probs = seq(0,1, 0.1), na.rm = T)
mirnatech$rank260 <-cut(mirnatech$`_260_280`, quant260)
miRNA_pheno4 <- merge(pheno, mirnatech, by.x = "shareid", by.y = "shareid")
miRNA_pheno <- merge(miRNA_pheno4, miRNAdat, by.x = "shareid", by.y = "shareid")
miRNA_pheno$Isolation_Batch<-as.factor(miRNA_pheno$Isolation_Batch)

miRNA_pheno$miR19a <- miRNA_pheno$miR_19a_3p
miRNA_pheno$miR186 <- miRNA_pheno$miR_186_5p_a2

ggplot(data = miRNA_pheno, aes(x = f8cbtobmd, y = miR19a, color = BB)) + geom_smooth(method = lm, se = TRUE, fullrange = TRUE) + ggtitle(label =  "miR-19a-3p by Total Femur BMD") + 
  xlab("Total Femur BMD") + ylab ("| \u0394 Cq |") + theme_minimal() + scale_color_brewer(palette = "Reds") + theme(plot.title = element_text(hjust = 0.5)) +theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10), title = element_text(size = 10), legend.text = element_text(size = 10), legend.position = c(0.1,0.8), legend.background = element_rect(fill="white",colour = "white")) 
ggsave(filename = "miR19byTFbmd.tiff", width = 90, height = 60 ,unit = "mm")

ggplot(data = miRNA_pheno, aes(x = f8cbtobmd, y = miR186, color = BB)) + geom_smooth(method = lm, se = TRUE, fullrange = TRUE) + ggtitle(label =  "miR-186-5p by Total Femur BMD") + 
  xlab("Total Femur BMD") + ylab ("| \u0394 Cq |") + theme_minimal() + scale_color_brewer(palette = "Reds") + 
  theme(plot.title = element_text(hjust = 0.5)) + theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10), title = element_text(size = 10), legend.text = element_text(size = 10), legend.position = c(0.1,0.8), legend.background = element_rect(fill="white",colour = "white")) 
ggsave(filename = "miR186byTFbmd.tiff", width = 90, height = 60 ,unit = "mm")