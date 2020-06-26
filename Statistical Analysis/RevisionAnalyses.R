# Additional Analyses for Revisions
rm(list = ls())
library(tidyverse)
library(haven)
library(emmeans)
library(pwr)

setwd('/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/data')
miRNAdat <- read.csv('l_mrna_2011_m_0797s_17_c1.csv')
pheno <- read.csv('PhenoData_5_28.csv')
mir_tech = "mirna_tech_17.sas7bdat"
mirnatech <- read_sas(mir_tech)
ldl <- read.csv("LDL_cov.csv")

cor.test.plus <- function(x) {
  list(x, 
       Standard.Error = unname(sqrt((1 - x$estimate^2)/x$parameter)))
}

# Determine T-score differences between BB user vs non-users ####
pheno <- merge(pheno, ldl, by = "shareid")
pheno$FN_Tscore <- ((-0.023 + 0.939 * pheno$f8cbnbmd  - 0.019)/1.087 - 0.858) / 0.120

mean(pheno$FN_Tscore)
sd(pheno$FN_Tscore)

mean(pheno$FN_Tscore[pheno$BB == "Yes"])
sd(pheno$FN_Tscore[pheno$BB == "Yes"])

mean(pheno$FN_Tscore[pheno$BB == "No"])
sd(pheno$FN_Tscore[pheno$BB == "No"])

fnmodel <- lm(FN_Tscore~BB + AGE8 + SEX + HGT8 + WGT8, data = pheno)
summary(fnmodel)
fn<-emmeans(fnmodel, specs = "BB")

fn <- data.frame(fn, BMD = "Femoral Neck")

# New Boxplots ####
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

boxplot(miR_19a_3p~BB,data = miRNA_pheno)
miRNA_pheno$miR19a <- miRNA_pheno$miR_19a_3p
miRNA_pheno$miR186 <- miRNA_pheno$miR_186_5p_a2
boxplot(miR19a~BB, data = miRNA_pheno, xlab = "BB User Status", ylab = "| \u0394 Cq |", main = "miR-19a-3p by BB Use")

miR19adata <- miRNA_pheno %>% select(.,miR19a, BB )
miR186data <- miRNA_pheno %>% select(.,miR186, BB )

p <- ggplot(miR19adata, aes(x = BB, y = miR19a, fill= BB)) + geom_boxplot() + scale_fill_brewer(palette = "Reds") + theme_classic()
p + ylab("| \u0394 Cq |") + theme(text = element_text(size = 18)) + theme(plot.title = element_text(hjust = 0.5)) +ggtitle("miR-19a-3p Expression")


p <- ggplot(miR186data, aes(x = BB, y = miR186, fill= BB)) + geom_boxplot() + scale_fill_brewer(palette = "Reds") + theme_classic()
p + ylab("| \u0394 Cq |") + theme(text = element_text(size = 18)) + theme(plot.title = element_text(hjust = 0.5)) +ggtitle("miR-186-5p Expression")

# Correlation between t-score/BMD and top miRNA ####
cor.test(miRNA_pheno$miR19a, miRNA_pheno$FN_Tscore, use = "complete", method = "spearman")
cor.test(miRNA_pheno$miR186, miRNA_pheno$FN_Tscore, use = "complete", method = "spearman")

cor.test(miRNA_pheno$miR19a, miRNA_pheno$f8cbtobmd, use = "complete", method = "spearman")
cor.test(miRNA_pheno$miR186, miRNA_pheno$f8cbtobmd, use = "complete", method = "spearman")

cor.test(miRNA_pheno$miR19a, miRNA_pheno$s8cbl24bd, use = "complete", method = "spearman")
cor.test(miRNA_pheno$miR186, miRNA_pheno$s8cbl24bd, use = "complete", method = "spearman")

model19aBB <- lm(miR19a ~ BB + AGE8 + SEX + HGT8 + WGT8 + rank260 + rankcon + rankqual, data = miRNA_pheno)
summary(model19aBB)

# Sensitivity analysis adjusting for treatment for hypertension ####
miRNA_pheno$HRX_BB[miRNA_pheno$HRX8 == "Yes" & miRNA_pheno$BB == "Yes"] <-"BB Hypertension Treatment"
miRNA_pheno$HRX_BB[miRNA_pheno$HRX8 == "Yes" & miRNA_pheno$BB == "No"] <-"Non-BB Hypertension Treatment"
miRNA_pheno$HRX_BB[miRNA_pheno$HRX8 == "No"] <-"Non-Hypertension Treatment"
miRNA_pheno$HRX_BB <- relevel(miRNA_pheno$HRX_BB, "Non-Hypertension")
table(miRNA_pheno$HRX_BB)
miR19a_BB_hrx <- lm(miR19a ~ HRX_BB + AGE8 + SEX + HGT8 + WGT8 + rank260 + rankcon + rankqual, data = miRNA_pheno)
summary(miR19a_BB_hrx)
miR186_BB_hrx <- lm(miR186 ~ HRX_BB + AGE8 + SEX + HGT8 + WGT8 + rank260 + rankcon + rankqual, data = miRNA_pheno)
summary(miR186_BB_hrx)

miR19a_Fto_hrx <- lm(miR19a ~ f8cbtobmd + AGE8 + SEX + HGT8 + WGT8 + rank260 + rankcon + rankqual + HRX_BB , data = miRNA_pheno)
summary(miR19a_Fto_hrx)
miR186_Fto_hrx <- lm(miR186 ~ f8cbtobmd + AGE8 + SEX + HGT8 + WGT8 + rank260 + rankcon + rankqual + HRX_BB , data = miRNA_pheno)
summary(miR186_Fto_hrx)

miR19a_LS_hrx <- lm(miR19a ~ s8cbl24bd + AGE8 + SEX + HGT8 + WGT8 + rank260 + rankcon + rankqual + HRX_BB , data = miRNA_pheno)
summary(miR19a_LS_hrx)
miR186_LS_hrx <- lm(miR186 ~ s8cbl24bd + AGE8 + SEX + HGT8 + WGT8 + rank260 + rankcon + rankqual + HRX_BB , data = miRNA_pheno)
summary(miR186_LS_hrx)

chisq.test(miRNA_pheno$BB, miRNA_pheno$HRX8)

# miRNAs in only hypertension treated patients ####
miRNA_pheno_tr_ht <- miRNA_pheno[miRNA_pheno$HRX8 == "Yes",]
miR19a_BB <- lm(miR19a ~ BB + AGE8 + SEX + HGT8 + WGT8 + rank260 + rankcon + rankqual, data = miRNA_pheno_tr_ht)
summary(miR19a_BB)
miR186_BB <- lm(miR186 ~ BB + AGE8 + SEX + HGT8 + WGT8 + rank260 + rankcon + rankqual , data = miRNA_pheno_tr_ht)
summary(miR186_BB)

miR19a_Fto <- lm(miR19a ~ f8cbtobmd + AGE8 + SEX + HGT8 + WGT8 + rank260 + rankcon + rankqual , data = miRNA_pheno_tr_ht)
summary(miR19a_Fto)
miR186_Fto <- lm(miR186 ~ f8cbtobmd + AGE8 + SEX + HGT8 + WGT8 + rank260 + rankcon + rankqual , data = miRNA_pheno_tr_ht)
summary(miR186_Fto)

miR19a_LS <- lm(miR19a ~ s8cbl24bd + AGE8 + SEX + HGT8 + WGT8 + rank260 + rankcon + rankqual , data = miRNA_pheno_tr_ht)
summary(miR19a_LS)
miR186_LS <- lm(miR186 ~ s8cbl24bd + AGE8 + SEX + HGT8 + WGT8 + rank260 + rankcon + rankqual , data = miRNA_pheno_tr_ht)
summary(miR186_LS)

# miRNA ~ B1 ####
miR19a_B1 <- lm(miR19a ~ B1 + AGE8 + SEX + HGT8 + WGT8 + rank260 + rankcon + rankqual , data = miRNA_pheno)
summary(miR19a_B1)
miR186_B1 <- lm(miR186 ~ B1 + AGE8 + SEX + HGT8 + WGT8 + rank260 + rankcon + rankqual , data = miRNA_pheno)
summary(miR186_B1)

# Trivariate BB use ####
miRNA_pheno$BB_3[miRNA_pheno$B1 == "Yes"] <- "B1 User"
miRNA_pheno$BB_3[(miRNA_pheno$BB == "Yes") & (miRNA_pheno$B1 == "No")] <- "Non-specific BB user"
miRNA_pheno$BB_3[miRNA_pheno$BB == "No"] <- "Non BB-user"
miRNA_pheno$BB_3 <- as.factor(miRNA_pheno$BB_3)
table(miRNA_pheno$BB_3)

miR19a_BB_3 <- lm(miR19a ~ BB_3 + AGE8 + SEX + HGT8 + WGT8 + rank260 + rankcon + rankqual , data = miRNA_pheno)
summary(miR19a_B1)
miR186_B1 <- lm(miR186 ~ B1 + AGE8 + SEX + HGT8 + WGT8 + rank260 + rankcon + rankqual , data = miRNA_pheno)
summary(miR186_B1)

# Sample size determination/Power Calculation ####
sd<- sd(pheno$f8cbtobmd)
toF_sd_miRNA <-sd(miRNA_pheno$f8cbtobmd)
table(pheno$BB)
sd_miRNA<-sd(na.omit(miRNA_pheno$miR19a))
table(miRNA_pheno$BB)
sd(na.omit(miRNA_pheno$miR186))
table(is.na(miRNA_pheno$miR19a))
table(is.na(miRNA_pheno$miR186))

# power for BB effect on TF BMD
# BB users in total cohort: 411
# BB non-users in total cohort: 1234 users

d = pwr.t2n.test(n1=411,n2=1234,d=,power=0.8)$d
eff = d*sd
print(eff)

d = pwr.t2n.test(n1=411,n2=1234,d=,power=0.9)$d
eff = d*sd
print(eff)
# 0.02880039 in the BMD difference between BB and non-BB we can detect with 90% power

# power for BB effect on miRNA
#SD for miR-19a-3p: 2.883572, 53 NAs (NAs omitted from SD calculation)
#SD for miR-186-5p: 2.526584, 24 NAs (NAs omitted from SD calculation)
# BB users in overlapping cohort: 357
# BB non-users in overlapping cohort: 1059
n1 = 357
n2 = 1059
d = pwr.t2n.test(n1=n1,n2=n2,d=,power=0.9)$d
eff = d*sd_miRNA
print(eff)
# BB use vs. no BB use
# 0.572427 is the miRNA difference between BB and non-BB we have 90% power to detect

# correlation between BMD and miRNA
r = pwr.r.test(n=(n1+n2),power=0.9,r=)$r
sd_bmd = sd
r*sd_miRNA/sd_bmd

# 1.590035 is the beta coefficient for bmd on miRNA that can be detected with 90% power

# Adjusting for Blood Pressure ####
miR19a_BB_BP <- lm(miR19a ~ BB + AGE8 + SEX + HGT8 + WGT8 + rank260 + rankcon + rankqual + SBP8 , data = miRNA_pheno)
summary(miR19a_BB_BP)
miR186_BB_BP <- lm(miR186 ~ BB + AGE8 + SEX + HGT8 + WGT8 + rank260 + rankcon + rankqual + SBP8 , data = miRNA_pheno)
summary(miR186_BB_BP)

miR19a_Fto_BP <- lm(miR19a ~ f8cbtobmd + AGE8 + SEX + HGT8 + WGT8 + rank260 + rankcon + rankqual + SBP8 , data = miRNA_pheno)
summary(miR19a_Fto_BP)
miR186_Fto_BP <- lm(miR186 ~ f8cbtobmd + AGE8 + SEX + HGT8 + WGT8 + rank260 + rankcon + rankqual + SBP8 , data = miRNA_pheno)
summary(miR186_Fto_BP)

miR19a_LS_BP <- lm(miR19a ~ s8cbl24bd + AGE8 + SEX + HGT8 + WGT8 + rank260 + rankcon + rankqual + SBP8 , data = miRNA_pheno)
summary(miR19a_LS_BP)
miR186_LS_BP <- lm(miR186 ~ s8cbl24bd + AGE8 + SEX + HGT8 + WGT8 + rank260 + rankcon + rankqual + SBP8 , data = miRNA_pheno)
summary(miR186_LS_BP)

miR19a_BB_DBP <- lm(miR19a ~ BB + AGE8 + SEX + HGT8 + WGT8 + rank260 + rankcon + rankqual + DBP8 , data = miRNA_pheno)
summary(miR19a_BB_DBP)
miR186_BB_DBP <- lm(miR186 ~ BB + AGE8 + SEX + HGT8 + WGT8 + rank260 + rankcon + rankqual + DBP8 , data = miRNA_pheno)
summary(miR186_BB_DBP)

miR19a_Fto_DBP <- lm(miR19a ~ f8cbtobmd + AGE8 + SEX + HGT8 + WGT8 + rank260 + rankcon + rankqual + DBP8 , data = miRNA_pheno)
summary(miR19a_Fto_DBP)
miR186_Fto_DBP <- lm(miR186 ~ f8cbtobmd + AGE8 + SEX + HGT8 + WGT8 + rank260 + rankcon + rankqual + DBP8 , data = miRNA_pheno)
summary(miR186_Fto_DBP)

miR19a_LS_DBP <- lm(miR19a ~ s8cbl24bd + AGE8 + SEX + HGT8 + WGT8 + rank260 + rankcon + rankqual + DBP8 , data = miRNA_pheno)
summary(miR19a_LS_DBP)
miR186_LS_DBP <- lm(miR186 ~ s8cbl24bd + AGE8 + SEX + HGT8 + WGT8 + rank260 + rankcon + rankqual + DBP8 , data = miRNA_pheno)
summary(miR186_LS_DBP)

# Adjusting for Both BP and Treatment for hypertension ####

miR19a_BB_BP3 <- lm(miR19a ~ BB + AGE8 + SEX + HGT8 + WGT8 + rank260 + rankcon + rankqual + SBP8 + DBP8 + HRX_BB  , data = miRNA_pheno)
summary(miR19a_BB_BP3)
miR186_BB_BP3 <- lm(miR186 ~ BB + AGE8 + SEX + HGT8 + WGT8 + rank260 + rankcon + rankqual + SBP8 + DBP8 + HRX_BB , data = miRNA_pheno)
summary(miR186_BB_BP3)

miR19a_Fto_BP3 <- lm(miR19a ~ f8cbtobmd + AGE8 + SEX + HGT8 + WGT8 + rank260 + rankcon + rankqual + SBP8 + DBP8 + HRX_BB , data = miRNA_pheno)
summary(miR19a_Fto_BP3)
miR186_Fto_BP3 <- lm(miR186 ~ f8cbtobmd + AGE8 + SEX + HGT8 + WGT8 + rank260 + rankcon + rankqual + SBP8 + DBP8 + HRX_BB , data = miRNA_pheno)
summary(miR186_Fto_BP3)

miR19a_LS_BP3 <- lm(miR19a ~ s8cbl24bd + AGE8 + SEX + HGT8 + WGT8 + rank260 + rankcon + rankqual + SBP8 + DBP8 + HRX_BB , data = miRNA_pheno)
summary(miR19a_LS_BP3)
miR186_LS_BP3 <- lm(miR186 ~ s8cbl24bd + AGE8 + SEX + HGT8 + WGT8 + rank260 + rankcon + rankqual + SBP8 + DBP8 + HRX_BB , data = miRNA_pheno)
summary(miR186_LS_BP3)

# Correlation between BP and Treatment for Hypertension/BB use ####
t.test(SBP8~HRX8, pheno)
t.test(DBP8~HRX8, pheno)
t.test(SBP8~BB, pheno)
t.test(DBP8~BB, pheno)
t.test(SBP8~B1, pheno)
t.test(DBP8~B1, pheno)

sbp_hrx_lm <- lm(SBP8~HRX8+AGE8+SEX+HGT8+WGT8, data = pheno)
summary(sbp_hrx_lm)
dbp_hrx_lm <- lm(DBP8~HRX8+AGE8+SEX+HGT8+WGT8, data = pheno)
summary(dbp_hrx_lm)

sbp_BB_lm <- lm(SBP8~BB+AGE8+SEX+HGT8+WGT8, data = pheno)
summary(sbp_BB_lm)
dbp_BB_lm <- lm(DBP8~BB+AGE8+SEX+HGT8+WGT8, data = pheno)
summary(dbp_BB_lm)

sbp_B1_lm <- lm(SBP8~B1+AGE8+SEX+HGT8+WGT8, data = pheno)
summary(sbp_B1_lm)
dbp_B1_lm <- lm(DBP8~B1+AGE8+SEX+HGT8+WGT8, data = pheno)
summary(dbp_BB_lm)

# Explained Variance ####
ToF_lm_miR19 = lm(f8cbtobmd~miR19a + AGE8 + SEX + HGT8 + WGT8 + rank260 + rankcon + rankqual, data = miRNA_pheno)
summary(ToF_lm_miR19)
summary(ToF_lm_miR19)$r.squared
mod = anova(ToF_lm_miR19)
ss = mod$`Sum Sq`
tot_rsq = sum(ss[1:length(ss)-1])/sum(ss)
print(tot_rsq)
x1_rsq = ss[1]/sum(ss)
print(x1_rsq)

# % explained variance for miR-19a#
exv_19 <- x1_rsq/tot_rsq

LS_lm_miR19 = lm(s8cbl24bd~miR19a + AGE8 + SEX + HGT8 + WGT8 + rank260 + rankcon + rankqual, data = miRNA_pheno)
summary(LS_lm_miR19)
summary(LS_lm_miR19)$r.squared
mod = anova(LS_lm_miR19)
ss = mod$`Sum Sq`
tot_rsq = sum(ss[1:length(ss)-1])/sum(ss)
print(tot_rsq)

x1_rsq = ss[1]/sum(ss)
print(x1_rsq)

# % explained variance for miR-19a#
exv_19 <- x1_rsq/tot_rsq
print(exv_19)

ToF_lm_miR186 = lm(f8cbtobmd~miR186 + AGE8 + SEX + HGT8 + WGT8 + rank260 + rankcon + rankqual, data = miRNA_pheno)
summary(ToF_lm_miR186)
summary(ToF_lm_miR186)$r.squared
mod = anova(ToF_lm_miR186)
ss = mod$`Sum Sq`
tot_rsq = sum(ss[1:length(ss)-1])/sum(ss)
print(tot_rsq)
x1_rsq = ss[1]/sum(ss)
print(x1_rsq)

# % explained variance for miR-186#
exv_186 <- x1_rsq/tot_rsq
print(exv_186)

LS_lm_miR186 = lm(s8cbl24bd~miR186 + AGE8 + SEX + HGT8 + WGT8 + rank260 + rankcon + rankqual, data = miRNA_pheno)
summary(LS_lm_miR186)
summary(LS_lm_miR186)$r.squared
mod = anova(LS_lm_miR186)
ss = mod$`Sum Sq`
tot_rsq = sum(ss[1:length(ss)-1])/sum(ss)
print(tot_rsq)

x1_rsq = ss[1]/sum(ss)
print(x1_rsq)

# % explained variance for miR-186#
exv_186 <- x1_rsq/tot_rsq
print(exv_186)

# Correlation between miRNA
cor.test(miRNA_pheno$miR19a, miRNA_pheno$miR186, method = "spearman")

# % difference in BMD between BB users and non-users ####
tscore_yes<- mean(pheno$FN_Tscore[pheno$BB == "Yes"])
tscore_no <- mean(pheno$FN_Tscore[pheno$BB == "No"])
pc_diff <- (tscore_yes - tscore_no) /tscore_no 
pc_diff
fnmodel_T <- lm(FN_Tscore~BB + AGE8 + SEX + HGT8 + WGT8, data = pheno)

fn<-emmeans(fnmodel, specs = "BB")
fn <- data.frame(fn, BMD = "Femoral Neck")
(fn$emmean[2] - fn$emmean[1])/fn$emmean[1]

fnmodel <- lm(f8cbnbmd~BB + AGE8 + SEX + HGT8 + WGT8, data = pheno)
ftomodel<-lm(f8cbtobmd~BB + AGE8 + SEX + HGT8 + WGT8, data = pheno)
ftrmodel <- lm(f8cbtrbmd~BB + AGE8 + SEX + HGT8 + WGT8, data = pheno)
s2model <- lm(s8cbl2bd~BB + AGE8 + SEX + HGT8 + WGT8, data = pheno)
s3model <- lm(s8cbl3bd~BB + AGE8 + SEX + HGT8 + WGT8, data = pheno)
s4model <- lm(s8cbl4bd~BB + AGE8 + SEX + HGT8 + WGT8, data = pheno)
s24model <- lm(s8cbl24bd~BB + AGE8 + SEX + HGT8 + WGT8, data = pheno)

summary(fnmodel)
summary(ftomodel)
summary(ftrmodel)
summary(s2model)
summary(s3model)
summary(s4model)
summary(s24model)

fn<-as.data.frame(emmeans(fnmodel, specs = "BB"))
fto<-as.data.frame(emmeans(ftomodel, specs = "BB"))
ftr<-as.data.frame(emmeans(ftrmodel, specs = "BB"))
s24<-as.data.frame(emmeans(s24model, specs = "BB"))

(fn$emmean[2] - fn$emmean[1])/fn$emmean[1]
(fto$emmean[2] - fto$emmean[1])/fto$emmean[1]
(ftr$emmean[2] - ftr$emmean[1])/ftr$emmean[1]
(s24$emmean[2] - s24$emmean[1])/s24$emmean[1]

# Analyses for meta analysis
pheno_fem <- miRNA_pheno[miRNA_pheno$SEX == "Female",]

summary(lm(miR186~f8cbtobmd+AGE8+HGT8+WGT8,data = pheno_fem))
summary(lm(miR19a~f8cbtobmd+AGE8+HGT8+WGT8,data = pheno_fem))

summary(lm(miR186~s8cbl24bd+AGE8+HGT8+WGT8,data = pheno_fem))
summary(lm(miR19a~s8cbl24bd+AGE8+HGT8+WGT8,data = pheno_fem))

summary(lm(miR186~FN_Tscore+AGE8+HGT8+WGT8,data = pheno_fem))
summary(lm(miR19a~FN_Tscore+AGE8+HGT8+WGT8,data = pheno_fem))

summary(lm(miR186~FN_Tscore+SEX+AGE8+HGT8+WGT8,data = miRNA_pheno))
summary(lm(miR19a~FN_Tscore+SEX+AGE8+HGT8+WGT8,data = miRNA_pheno))

cor.test.plus(cor.test(miRNA_pheno$FN_Tscore, miRNA_pheno$miR19a))
cor.test.plus(cor.test(miRNA_pheno$FN_Tscore, miRNA_pheno$miR186))

summary(lm(f8cbtobmd~miR186+AGE8+HGT8+WGT8,data = pheno_fem))
summary(lm(f8cbtobmd~miR19a+AGE8+HGT8+WGT8,data = pheno_fem))
sd(miRNA_pheno$miR_19a_3p, na.rm = T)
sd(miRNA_pheno$miR_186_5p_a2, na.rm = T)
