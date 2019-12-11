#Phenotype Data#
# Setup R ####
rm(list = ls())
library(tidyverse)
setwd("/home/clary@mmcf.mehealth.org/Framingham")
# Read in Data ####
offdate <- read.delim("phs000007.v30.pht003099.v5.p11.c1.vr_dates_2014_a_0912s.HMB-IRB-MDS.txt",
                      skip=10,header=T,stringsAsFactors = F)
cov <- read.delim("phs000007.v30.pht006027.v2.p11.c1.vr_wkthru_ex09_1_1001s.HMB-IRB-MDS.txt",
                  skip=10,header=T,stringsAsFactors = F)
est <- read.delim("phs000007.v30.pht000307.v8.p11.c1.meno1_8s.HMB-IRB-MDS.txt",
                  skip=10,header=T,stringsAsFactors = F)
mi <- read.delim("phs000007.v30.pht003335.v6.p11.c1.vr_soe4srv_2014_a_1027s.HMB-IRB-MDS.txt",
                     skip=10,header=T,stringsAsFactors = F)
drug <- read.delim("phs000007.v30.pht000828.v6.p11.c1.meds1_8s.HMB-IRB-MDS.txt", 
                   skip=10,header=T,stringsAsFactors = F)
bmd <- read.delim("phs000007.v30.pht003096.v2.p11.c1.t_bmdhs_2008_1_0748s.HMB-IRB-MDS.txt",
                  skip=10,header=T,stringsAsFactors = F)

# Filter Exam 8 people ####
off8 <- offdate %>% 
  select(shareid,age8,att8,date8,idtype) %>% 
  filter(idtype==1) %>% # In Offspring Cohort
  filter(att8==1) # Attended Exam 8
# Exam 8 covariates
cov8 = cov %>% # 
  select(shareid,SEX,AGE8,BMI8,HGT8,WGT8,CURRSMK8,CPD8,DMRX8,HRX8,LIPRX8,IDTYPE) %>% 
  filter(IDTYPE==1) # In Offspring Cohort
# Exam 8 Estrogen Levels
est8 <- est %>%
  select(shareid,OVREM8,EST8,STOP_AGE)
# MI Offspring Cohort
mi_off <- mi[mi$IDTYPE ==1,]
# BB Use
drug <- drug %>% 
  select(shareid,ther_gp1,MEDNAME,MEDPRN,chem_gp1,chem_nm1) %>% 
  filter(ther_gp1=="BETA BLOCKING AGENTS") %>% 
  filter(!(MEDNAME %in% c("ISTALOL","IC ATENOLOL"))) %>% 
  filter(MEDPRN==0) %>% # 4 RN, 50 unknown
  mutate(B1drug=ifelse(chem_gp1=="Beta blocking agents, selective",1,0)) %>% 
  mutate(BBdrug=1)  %>% 
  group_by(shareid) %>%
  summarize(BB=max(BBdrug),B1=max(B1drug))
# BMD data
bmd8 <- merge(bmd,off8,by.x="shareid",by.y="shareid")
bmd8 <- bmd8 %>%
  mutate(f_dpast8 = f8cbscdt-date8) %>% #Days femurs scan after exam 8 date
  mutate(s_dpast8 = s8cbscdt-date8) %>% # Days spine scan after exam 8 date
  mutate(mis_f8 = is.na(f8cbscdt) | f_dpast8<0) %>% # Is femur data before exam 8 or missing
  mutate(mis_s8 = is.na(s8cbscdt) | s_dpast8<0) %>% # Is spine data before exam 8 or missing
  filter(!mis_f8 | !mis_s8) # Filter for non missing or before exam 8 femur of spine data

# Merge Phenotype data together ####
cohort <- merge(off8,cov8,by.x="shareid",by.y="shareid")
cohort <- merge(cohort,est8,by.x="shareid",by.y="shareid",all.x=T)
cohort <- merge(cohort, drug,by.x="shareid",by.y="shareid",all.x=T)
cohort <- merge(cohort, bmd8, by.x ="shareid",by.y="shareid")
cohort$BB[is.na(cohort$BB)] = 0
cohort$B1[is.na(cohort$B1)] = 0

# Create menov variable ####
cohort <- cohort %>% mutate(menov = (OVREM8==2 | STOP_AGE<age8.x)) # 2 ovaries removed or menopause
cohort$menov[cohort$SEX==1] <- 2 # If male <- 2
cohort$EST8[cohort$SEX==1] <- 2 # If male <- 2

# Create priorcvd variable ####
cohortmi <- merge(cohort,mi_off,by.x="shareid",by.y="shareid",all.x=T)
cvd <- cohortmi %>%
  mutate(cvd_event = ifelse(is.na(EVENT),0,ifelse(EVENT %in% c(1:19,30:49),1,0))) %>%
  mutate(event_days8 = DATE-date8.x) %>%
  mutate(priorcvd_event = ifelse(event_days8<0 & cvd_event==1,1,0)) %>%
  group_by(shareid) %>%
  summarize(priorcvd = max(priorcvd_event))
cohort <- merge(cohort, cvd,by.x="shareid",by.y="shareid",all.x=T)
cohort$priorcvd[is.na(cohort$priorcvd)] <- 0 # Set NAs to 0

pheno<- cohort

# Subset analytic dataset/phenotype data ####
keep <- c("shareid", "SEX", "AGE8", "HGT8", "WGT8", "CURRSMK8", "DMRX8", "HRX8",
          "LIPRX8", "priorcvd", "menov","EST8", "BB", "B1", "f8cbnbmd",
          "s8cbl24bd", "f8cbtobmd", "f8cbtrbmd", "s8cbl2bd","s8cbl3bd","s8cbl4bd")
pheno1 <- pheno[,names(pheno) %in% keep]
pheno3 <- pheno1[is.na(pheno1$f8cbnbmd) == F,] # Filter out if bmd is missing
pheno3$DMRX8[is.na(pheno3$DMRX8) == T]<- 2 #Missing Data set to 2
pheno3$EST8[is.na(pheno3$EST8) == T]<- 3 # Missing data set to 3

# Set necessary variables as factors ####
pheno3$SEX<- as.factor(pheno3$SEX)
levels(pheno3$SEX)<- c("Male", "Female")

pheno3$CURRSMK8 <- as.factor(pheno3$CURRSMK8)
levels(pheno3$CURRSMK8)<- c("No", "Yes")

pheno3$DMRX8 <- as.factor(pheno3$DMRX8)
levels(pheno3$DMRX8)<- c("No", "Yes", "Missing")

pheno3$HRX8 <- as.factor(pheno3$HRX8)
levels(pheno3$HRX8)<- c("No", "Yes")

pheno3$LIPRX8 <- as.factor(pheno3$LIPRX8)
levels(pheno3$LIPRX8) <-c("No", "Yes")

pheno3$EST8 <- as.factor(pheno3$EST8)
levels(pheno3$EST8)<- c("No", "Yes", "Male", "Missing")

pheno3$menov <- as.factor(pheno3$menov)
levels(pheno3$menov)<- c("No", "Yes", "Male")

pheno3$priorcvd <- as.factor(pheno3$priorcvd)
levels(pheno3$priorcvd) <-c("No", "Yes")

pheno3$BB <- as.factor(pheno3$BB)
levels(pheno3$BB) <-c("No", "Yes")

pheno3$B1 <- as.factor(pheno3$B1)
levels(pheno3$B1) <-c("No", "Yes")

pheno3$B1B <- pheno3$BB
levels(pheno3$B1B) <- c("No", "Non-selective", "B1-Selective")
pheno3$B1B[pheno3$B1 == "Yes"]<- "B1-Selective"
table(pheno3$B1B)
table(pheno3$BB)
setwd("/home/clary@mmcf.mehealth.org/Framingham/OmicData/data")
write.csv(pheno3,"PhenoData_5_28.csv", row.names = F)
