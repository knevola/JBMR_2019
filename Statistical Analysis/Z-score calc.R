# Calculate Z-score in FHS
rm(list = ls())
library(dplyr)
setwd('/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/data')

pheno_0 <- read.csv('PhenoData_5_28.csv')
miRNAdat <- read.csv('l_mrna_2011_m_0797s_17_c1.csv')
# Corr test plus function
cor.test.plus <- function(x) {
  list(x, 
       Standard.Error = unname(sqrt((1 - x$estimate^2)/x$parameter)))
}
pheno_fem_0 <- pheno_0[pheno_0$SEX == "Female",]
pheno_male_0 <- pheno_0[pheno_0$SEX == "Male",]

# Format miRNA data
miRNA_delta_cq <- miRNAdat[-1]
miRNA_delta_cq <- -(miRNA_delta_cq-27)
miRNAdat <- cbind(miRNAdat[1], miRNA_delta_cq)

# Z-score function
source("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/Aim1_Potential_Mechanism_scripts/Z-Score_function.R")
pheno_fem_1 <- z_score(age = "AGE8", bone = "f8cbnbmd", data = pheno_fem_0, age_block = 5)
pheno_male_1 <- z_score(age = "AGE8", bone = "f8cbnbmd", data = pheno_male_0, age_block = 5)
pheno_1 <- rbind(pheno_fem_1, pheno_male_1)

# Merge pheno and miRNA
miRNA_pheno <- merge(miRNAdat, pheno_1, by = "shareid")

cor.test.plus(cor.test(miRNA_pheno$miR_19a_3p, miRNA_pheno$zscore))
cor.test.plus(cor.test(miRNA_pheno$miR_186_5p_a2, miRNA_pheno$zscore))
