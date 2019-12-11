# 7/9 MLR with Technical Variables Rank Transformed
# 8/20 Added glm and filtered < 0.05
# Setup ####
rm(list = ls())
library(tidyverse)
library(dplyr)
library(VennDiagram)
library(metaseqR)
library(ggplot2)
library(pheatmap)
library(haven)
library(lmerTest)
setwd('/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/data')
miRNAdat <- read.csv('l_mrna_2011_m_0797s_17_c1.csv')
pheno <- read.csv('PhenoData_5_28.csv')
mir_tech = "mirna_tech_17.sas7bdat"
mirnatech <- read_sas(mir_tech)

miRNA_delta_cq <- miRNAdat[-1]
miRNA_delta_cq <- -(miRNA_delta_cq-27)
miRNAdat <- cbind(miRNAdat[1], miRNA_delta_cq)

# Merge pheno with miRNA ####
quantcon<-quantile(mirnatech$concentration, probs = seq(0,1, 0.1))
mirnatech$rankcon <- cut(mirnatech$concentration, quantcon)
quantqual <-quantile(mirnatech$RNA_quality, probs = seq(0,1, 0.1), na.rm = T)
mirnatech$rankqual <- cut(mirnatech$RNA_quality, quantqual)
quant260 <- quantile(mirnatech$`_260_280`, probs = seq(0,1, 0.1), na.rm = T)
mirnatech$rank260 <-cut(mirnatech$`_260_280`, quant260)
miRNA_pheno4 <- merge(pheno, mirnatech, by.x = "shareid", by.y = "shareid")
miRNA_pheno <- merge(miRNA_pheno4, miRNAdat, by.x = "shareid", by.y = "shareid")
miRNA_pheno$Isolation_Batch<-as.factor(miRNA_pheno$Isolation_Batch)
Colsums <- data.frame(variables = names(miRNA_pheno), NAs = colSums(is.na(miRNA_pheno))/nrow(miRNA_pheno))
Colsums$NAs[1:30]<- NA 
drop <- c('cvdpair', 'casecontrol', 'idtype.x', "idtype.y")
miRNA_pheno <- miRNA_pheno[, !(names(miRNA_pheno) %in% drop)]

# Setup model specific data ####
# Model 1 Only (linear)
Colsums1 <- Colsums[Colsums$NAs < 0.1 | is.na(Colsums$NAs) == T,]
miRNA_pheno1 <- miRNA_pheno[names(miRNA_pheno) %in% Colsums1$variables] # Filter out miRNA with more than 10% NAs
names1 <- names(miRNA_pheno1)[-c(1:29)]
# Model 2 only (log/binary)
Colsums2 <- Colsums[(Colsums$NAs > 0.9) & (Colsums$NAs < 0.95) | is.na(Colsums$NAs) == T,]
miRNA_pheno2 <- miRNA_pheno[names(miRNA_pheno) %in% Colsums2$variables]
miRNAdatna <- miRNA_pheno2[-c(1:29)]
miRNAdatna[!is.na(miRNAdatna)] <- 0 #Detectable
miRNAdatna[is.na(miRNAdatna)] <- 1 #Not Detectable
miRNAdatna[sapply(miRNAdatna, is.numeric)] <- lapply(miRNAdatna[sapply(miRNAdatna, is.numeric)], as.logical)
miRNA_pheno2 <- cbind(miRNA_pheno2[c(1:29)], miRNAdatna)
names2 <- names(miRNA_pheno2)[-c(1:29)]

# Model 3:Both Model 1 and 2
Colsums3 <- Colsums[(Colsums$NAs > 0.1 & Colsums$NAs < 0.9)| is.na(Colsums$NAs) == T,]
miRNA_pheno3 <- miRNA_pheno[names(miRNA_pheno) %in% Colsums3$variables]
miRNAdatna <- miRNA_pheno3[-c(1:29)]
miRNAdatna[!is.na(miRNAdatna)] <- 0 #Detectable
miRNAdatna[is.na(miRNAdatna)] <- 1 #Not Detectable
miRNAdatna[sapply(miRNAdatna, is.numeric)] <- lapply(miRNAdatna[sapply(miRNAdatna, is.numeric)], as.logical)
miRNA_pheno3_2 <- cbind(miRNA_pheno3[c(1:29)], miRNAdatna)
names3 <- names(miRNA_pheno3)[-c(1:29)]

# MLR Function ####
mlr<- function(names, cov, data, var){
  name <- NULL
  mlr1 <- NULL
  for (i in 1:length(names)){
    #tryCatch({
    print(i)
    x <- lm(as.formula(paste(names[i], cov, sep = "")), data = data)
    y <- summary(x)
    z<- as.data.frame(y$coefficients)
    z$variable <- row.names(z)
    row.names(z)<- NULL
    if (i == 1){
      mlr1 <- z
      name <- rep(names[i], nrow(z))
    }else{
      mlr1 <- rbind(mlr1, z)
      name <- c(name, rep(names[i], nrow(z)))
    }
    #}, error = function(e){})
  }
  mlr_results <- cbind(miRNA = name, mlr1)
  mlr_fdr <- mlr_results %>% filter(., variable == var)
  mlr_fdr$FDR <- p.adjust(p = mlr_fdr$`Pr(>|t|)`, method = "BH", n = length(names)) 
  mlr_results <- merge(x = mlr_results, y = mlr_fdr, all.x = T)
  return(mlr_results)
}

mlrg<- function(names, cov, data, var){
  name <- NULL
  mlr1 <- NULL
  for (i in 1:length(names)){
    #tryCatch({
    print(i)
    x <- glm(as.formula(paste(names[i], cov, sep = "")), data = data)
    y <- summary(x)
    z<- as.data.frame(y$coefficients)
    z$variable <- row.names(z)
    row.names(z)<- NULL
    if (i == 1){
      mlr1 <- z
      name <- rep(names[i], nrow(z))
    }else{
      mlr1 <- rbind(mlr1, z)
      name <- c(name, rep(names[i], nrow(z)))
    }
    #}, error = function(e){})
  }
  mlr_results <- cbind(miRNA = name, mlr1)
  mlr_fdr <- mlr_results %>% filter(., variable == var)
  mlr_fdr$FDR <- p.adjust(p = mlr_fdr$`Pr(>|t|)`, method = "BH", n = length(names)) 
  mlr_results <- merge(x = mlr_results, y = mlr_fdr, all.x = T)
  return(mlr_results)
}
fisherpvalue <- function(data1, data2, var){
  names(data1) <- c("miRNA", "Estimate", "SE", "tvalue", "pvalue", "variable", "FDR")
  names(data2) <- c("miRNA", "Estimate", "SE", "tvalue", "pvalue", "variable", "FDR")
  data1 <- data1 %>% filter(.,variable == var)
  data2 <- data2 %>%  filter(.,variable == var)
  model1a2<- merge(data1, data2, by.x = 'miRNA', by.y = "miRNA")
  row.names(model1a2)<- model1a2$miRNA
  keep1 <- c("pvalue.x","pvalue.y") 
  model1a26<- model1a2[names(model1a2) %in% keep1]
  results<-fisher.method(model1a26, p.corr = "BH")
  model3<-merge(model1a2,results, by.x = 0, by.y = 0)
  return(model3)
}

Ftocov <- "~f8cbtobmd + AGE8 + SEX + HGT8 + WGT8 + rankcon + rankqual + rank260"
S24cov <- "~s8cbl24bd + AGE8 + SEX + HGT8 + WGT8 + rankcon + rankqual + rank260"
BBcov <- "~BB + AGE8 + SEX + HGT8 + WGT8 + rankcon + rankqual + rank260"

Ftoresults1 <- mlr(names = names1, cov = Ftocov, data = miRNA_pheno1, var = "f8cbtobmd")
S24results1 <- mlr(names = names1, cov = S24cov, data = miRNA_pheno1, var = "s8cbl24bd")
BBresults1 <- mlr(names = names1, cov = BBcov, data = miRNA_pheno1, var = "BBYes")

Ftoresults2 <- mlrg(names = names2, cov = Ftocov, data = miRNA_pheno2, var = "f8cbtobmd")
S24results2 <- mlrg(names = names2, cov = S24cov, data = miRNA_pheno2, var = "s8cbl24bd")
BBresults2 <- mlrg(names = names2, cov = BBcov, data = miRNA_pheno2, var = "BBYes")

Ftoresults3_1 <- mlr(names = names3, cov = Ftocov, data = miRNA_pheno3, var = "f8cbtobmd")
S24results3_1 <- mlr(names = names3, cov = S24cov, data = miRNA_pheno3, var = "s8cbl24bd")
BBresults3_1 <- mlr(names = names3, cov = BBcov, data = miRNA_pheno3, var = "BBYes")

Ftoresults3_2 <- mlrg(names = names3, cov = Ftocov, data = miRNA_pheno3_2, var = "f8cbtobmd")
S24results3_2 <- mlrg(names = names3, cov = S24cov, data = miRNA_pheno3_2, var = "s8cbl24bd")
BBresults3_2 <- mlrg(names = names3, cov = BBcov, data = miRNA_pheno3_2, var = "BBYes")

Ftoresults3 <- fisherpvalue(data1 = Ftoresults3_1, data2 = Ftoresults3_2, var = "f8cbtobmd")
S24results3 <- fisherpvalue(data1 = S24results3_1, data2 = S24results3_2, var = "s8cbl24bd")
BBresults3 <- fisherpvalue(data1 = BBresults3_1, data2 = BBresults3_2, var = "BBYes")

Ftoresults1_p <- Ftoresults1 %>%  filter(., variable == "f8cbtobmd") %>% filter(., `Pr(>|t|)` < 0.05)
#Ftoresults1_f <- Ftoresults1 %>%  filter(., variable == "f8cbtobmd") %>% filter(., FDR < 0.1)
S24results1_p <- S24results1 %>%  filter(., variable == "s8cbl24bd") %>% filter(., `Pr(>|t|)` < 0.05)
#S24results1_f <- S24results1 %>%  filter(., variable == "s8cbl24bd") %>% filter(., FDR < 0.1)
BBresults1_f <- BBresults1 %>% filter(., variable == "BBYes") %>% filter(.,FDR < 0.1)
BBresults1_p <- BBresults1 %>% filter(., variable == "BBYes") %>% filter(., `Pr(>|t|)` < 0.05)

Ftoresults2_p <- Ftoresults2 %>%  filter(., variable == "f8cbtobmd") %>% filter(., `Pr(>|t|)` < 0.05)
Ftoresults2_f <- Ftoresults2 %>%  filter(., variable == "f8cbtobmd") %>% filter(., FDR < 0.1)
S24results2_p <- S24results2 %>%  filter(., variable == "s8cbl24bd") %>% filter(., `Pr(>|t|)` < 0.05)
#S24results2_f <- S24results2 %>%  filter(., variable == "s8cbl24bd") %>% filter(., FDR < 0.1)
BBresults2_f <- BBresults2 %>% filter(., variable == "BBYes") %>% filter(.,FDR < 0.1)
BBresults2_p <- BBresults2 %>% filter(., variable == "BBYes") %>% filter(., `Pr(>|t|)` < 0.05)

Ftoresults3_p <- Ftoresults3  %>% filter(., p.value < 0.05)
#Ftoresults3_f <- Ftoresults3 %>%  filter(., p.adj < 0.1)
S24results3_p <- S24results3  %>% filter(., p.value < 0.05)
#S24results3_f <- S24results3 %>% filter(., p.adj < 0.1)
BBresults3_f <- BBresults3 %>% filter(.,p.adj < 0.1)
BBresults3_p <- BBresults3 %>% filter(.,p.value < 0.1)

Ftoall <- c(as.character(Ftoresults1_p$miRNA), as.character(Ftoresults2_p$miRNA), as.character(Ftoresults3_p$miRNA))
S24all <- c(as.character(S24results1_p$miRNA), as.character(S24results2_p$miRNA), as.character(S24results3_p$miRNA))
BBall <- c(as.character(BBresults1_f$miRNA), as.character(BBresults2_f$miRNA), as.character(BBresults3_f$miRNA))
BBall_p <- c(as.character(BBresults1_p$miRNA), as.character(BBresults2_p$miRNA), as.character(BBresults3_p$miRNA))

int<- intersect(intersect(Ftoall, S24all), BBall)
int_p<- intersect(intersect(Ftoall, S24all), BBall_p)

a <- data.frame(Ftoresults1_p[Ftoresults1_p$miRNA %in% int,], model = "linear")
#b <- data.frame(Ftoresults2_p[Ftoresults2_p$miRNA %in% int,], model = "logistic")
c <- data.frame(Ftoresults3_p[Ftoresults3_p$miRNA %in% int,2:8], model = "linear")
names(c)<-names(a)
d <- data.frame(Ftoresults3_p[Ftoresults3_p$miRNA %in% int,c(2,9:14)], model = "logistic")
names(d)<-names(a)
Femurtable<-rbind(a,c,d)
Femurtable <- Femurtable[order(Femurtable$miRNA),]
write.csv(Femurtable, "Femuroverlapresults_7_9_rank.csv")

a <- data.frame(S24results1_p[S24results1_p$miRNA %in% int,], model = "linear")
#b <- data.frame(S24results2_p[S24results2_p$miRNA %in% int,], model = "logistic")
c <- data.frame(S24results3_p[S24results3_p$miRNA %in% int,2:8], model = "linear")
names(c)<-names(a)
d <- data.frame(S24results3_p[S24results3_p$miRNA %in% int,c(2,9:14)], model = "logistic")
names(d)<-names(a)
Spinetable<-rbind(a,c,d)
Spinetable <- Spinetable[order(Spinetable$miRNA),]
write.csv(Spinetable, "Spineoverlapresults_7_9_rank.csv")


a <- data.frame(BBresults1_f[BBresults1_f$miRNA %in% int,], model = "linear")
#b <- data.frame(BBresults2_f[BBresults2_f$miRNA %in% int,], model = "logistic")
c <- data.frame(BBresults3_f[BBresults3_f$miRNA %in% int,2:8], model = "linear")
names(c)<-names(a)
d <- data.frame(BBresults3_f[BBresults3_f$miRNA %in% int,c(2,9:14)], model = "logistic")
names(d)<-names(a)
BBtable<-rbind(a,c,d)
BBtable <- BBtable[order(BBtable$miRNA),]
write.csv(BBtable, "BBoverlapresults_7_9_rank.csv")


setwd("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/Clustering_miRNA/Figures")
WGCNA <- read.csv("MiRNA_significance_module_membership_nofilter.csv")
setwd('/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/data')
WGCNA_select<-WGCNA %>%  filter(., mergedColors == "blue")
WGCNA_brown <- WGCNA %>% filter(., mergedColors == "brown")

listall <- list( `Total Femur \n BMD` = Ftoall, `Spine L2-L4 \n BMD` =S24all, `BB Use` = BBall)
listall_wgcna <- list(`Total Femur \n BMD` = Ftoall, `Spine L2-L4 \n BMD` =S24all, `BB Use` = BBall, WGCNA = as.character(WGCNA_select$X))
setwd('/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/figure')
venn.diagram(listall, filename = "SpineFemurBBVenn_7_9_rank.tiff", fill = c("slateblue1", "lightskyblue1", "violetred1"), main = "Total Femur BMD, Spine L2-L4 BMD, \n and BB use Overlap")
venn.diagram(listall_wgcna, filename = "SpineFemurBBVenn_WGCNA_7_9_rank.tiff", fill = c("slateblue1", "lightskyblue1", "violetred1", "goldenrod1"), main = "Total Femur BMD, Spine L2-L4 BMD, \n BB use and WGCNA Overlap")

# Heatmap ####
miRNA_interest <- int
intvariables <- c("SEX", "AGE8", "HGT8", "WGT8", "BB", "s8cbl24bd", "f8cbtobmd")
int1 <- c(int, intvariables)
pheno_int <-miRNA_pheno[names(miRNA_pheno) %in% int1]
pheno_int$BB <- as.numeric(pheno_int$BB)
pheno_int$SEX <- as.numeric(pheno_int$SEX)
pheno_int1 <- t(pheno_int)
pheatmap(pheno_int1, scale = "row")

intersect(int, WGCNA_select$X)

FemurInt <- intersect(Ftoall, BBall)
intersect(FemurInt, WGCNA_select$X)
SpineInt <- intersect(S24all, BBall)
intersect(SpineInt, WGCNA_select$X)
SpineFem <- intersect(Ftoall, S24all)
intersect(SpineFem, WGCNA_select$X)
intersect(Ftoall, WGCNA_select$X)
intersect(S24all, WGCNA_select$X)
intersect(BBall, WGCNA_select$X)

x <- list(`Total Femur \n BMD` = Ftoall, `Lumbar Spine \n BMD` = S24all, `BB Use` = BBall, `Blue WGCNA \n Cluster` = WGCNA_select$X)
venn.diagram(x, filename = "WGCNABMDBBVenn_8_19.tiff", fill = c("steelblue", "lightblue", "red", "goldenrod"))

y <- list(`Total Femur \n BMD` = Ftoall, `Lumbar Spine \n BMD` = S24all, `BB Use` = BBall, `Brown WGCNA \n Cluster` = WGCNA_brown$X)
venn.diagram(y, filename = "WGCNABMDBBVenn_8_19_brown.tiff", fill = c("steelblue", "lightblue", "red", "goldenrod"))

# Supplementary Table Results

#Filter BMD Results no significance filter
TFResults1 <- Ftoresults1 %>% filter(., variable == "f8cbtobmd")
TFResults2 <- Ftoresults2 %>% filter(., variable == "f8cbtobmd")
LSResults1 <- S24results1 %>% filter(., variable == "s8cbl24bd")
LSResults2 <- S24results2 %>% filter(., variable == "s8cbl24bd")
BBResults1 <- BBresults1 %>% filter(., variable == "BBYes")
BBResults2 <- BBresults2 %>% filter(., variable == "BBYes")

TFResults1$Phenotype <- "Total Femur"
TFResults1$Analysis <- "Linear"
TFResults2$Phenotype <- "Total Femur"
TFResults2$Analysis <- "Logistic"
LSResults1$Phenotype <- "Lumbar Spine"
LSResults1$Analysis <- "Linear"
LSResults2$Phenotype <- "Lumbar Spine"
LSResults2$Analysis <- "Logistic"
BBResults1$Phenotype <- "BB Use"
BBResults1$Analysis <- "Linear"
BBResults2$Phenotype <- "BB Use"
BBResults2$Analysis <- "Logistic"

LnLResults <- rbind(TFResults1, TFResults2, LSResults1, LSResults2, BBResults1, BBResults2)
Ftoresults3$Phenotype <- "Total Femur"
S24results3$Phenotype <- "Lumbar Spine"
BBresults3$Phenotype <- "BB Use"

Ftoresults3$Analysis <- "Combined"
S24results3$Analysis <- "Combined"
BBresults3$Analysis <- "Combined"

CombinedResults <- rbind(Ftoresults3,S24results3, BBresults3)

write.csv(LnLResults, "LinearAndLogisticResults.csv")
write.csv(CombinedResults, "CombinedResults.csv")
