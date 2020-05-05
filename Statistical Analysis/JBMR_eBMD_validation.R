# Working with Ines and Barbara's data
rm(list = ls())
options(stringsAsFactors = F)
setwd("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/data")
library(xlsx)
library(tidyverse)

miR_data <- read.xlsx("2020_01_16_miR19_186_in_SOS_hip_for_Kristine - Kopie (2).xlsx",sheetIndex = 1)
bmd_data <- read.xlsx("eBMD heel ultrsound for rebuttal_Ines 26042020_incl. upper limb data_age_sex.xlsx", sheetIndex = 1)
bmd_data$Sex <- as.factor(bmd_data$Sex)
column_names <- as.data.frame(miR_data[1:4,c(1,8:42)])
info<- as.data.frame(t(column_names))
names(info) <- info[1,]
info <- info[-1,]
info1<- separate(info,col = 1, into = c("Group", "ID"), sep = ":")
info1$ID <- as.numeric(info1$ID)
info1$`hsa-miR-19a-3p` <- as.numeric(info1$`hsa-miR-19a-3p`)
info1$`hsa-miR-186-5p` <- as.numeric(info1$`hsa-miR-186-5p`)


bmd_info <- merge(bmd_data, info1, by.x = "PatNr", by.y = "ID")

t.test(eBMD..Lunar.heel.ultrasound. ~ Group, bmd_info)

bmd_info_complete <- bmd_info[!is.na(bmd_info$eBMD..Lunar.heel.ultrasound.),]
table(bmd_info_complete$Group)
Tbmd_info_complete <- bmd_info[!is.na(bmd_info$T.Score..Lunar.heel.ultrasound.),]
table(Tbmd_info_complete$Group,Tbmd_info_complete$Sex)

t.test(T.Score..Lunar.heel.ultrasound. ~ Group, bmd_info)
t.test(Z.Score..Lunar.heel.ultrasound. ~ Group, bmd_info)
t.test(BUA..Lunar.heel.ultrasound. ~ Group, bmd_info)
t.test(SOS..Lunar.heel.ultrasound. ~ Group, bmd_info)

cor.test(bmd_info$eBMD..Lunar.heel.ultrasound., bmd_info$`hsa-miR-19a-3p`, method = "s")
cor.test(bmd_info$eBMD..Lunar.heel.ultrasound., bmd_info$`hsa-miR-186-5p`, method = "s")

cor.test(bmd_info$Z.Score..Lunar.heel.ultrasound., bmd_info$`hsa-miR-19a-3p`)
cor.test(bmd_info$Z.Score..Lunar.heel.ultrasound., bmd_info$`hsa-miR-186-5p`)

summary(lm(`hsa-miR-19a-3p`~eBMD..Lunar.heel.ultrasound. + Height.cm + Weight.kg + Age. , data = bmd_info))
summary(lm(`hsa-miR-186-5p`~eBMD..Lunar.heel.ultrasound. + Height.cm + Weight.kg + Age. , data = bmd_info))

summary(lm(`hsa-miR-19a-3p`~Z.Score..Lunar.heel.ultrasound. + Height.cm + Weight.kg + Age. + Sex, data = bmd_info))
summary(lm(`hsa-miR-186-5p`~Z.Score..Lunar.heel.ultrasound. + Height.cm + Weight.kg + Age. + Sex, data = bmd_info))

summary(lm(`hsa-miR-19a-3p`~T.Score..Lunar.heel.ultrasound. + Height.cm + Weight.kg + Age. + Sex, data = bmd_info))
summary(lm(`hsa-miR-186-5p`~T.Score..Lunar.heel.ultrasound. + Height.cm + Weight.kg + Age. + Sex, data = bmd_info))

hist(bmd_info$T.Score..Lunar.heel.ultrasound.)
hist(bmd_info$Age.)
mean(bmd_info$Age.)
mean_sdl(bmd_info$Age.)
min(bmd_info$Age.)
max(bmd_info$Age.)
bmd_info$TGroup[bmd_info$T.Score..Lunar.heel.ultrasound. > -1] <- "normal"
bmd_info$TGroup[bmd_info$T.Score..Lunar.heel.ultrasound. < -1] <- "osteopenia"
bmd_info$TGroup[bmd_info$T.Score..Lunar.heel.ultrasound. < -2.5] <- "osteoporosis"

summary(lm(`hsa-miR-19a-3p`~TGroup + Height.cm + Weight.kg + Age. + Sex, data = bmd_info))
summary(lm(`hsa-miR-186-5p`~TGroup+ Height.cm + Weight.kg + Age. + Sex, data = bmd_info))
summary(bmd_info$T.Score..Lunar.heel.ultrasound.)
