# miRNA TissueAtlas Bar Graph#
rm(list = ls())
setwd("/home/clary@mmcf.mehealth.org/Framingham/OmicData/data")
library(ggplot2)
tissue <- read.csv("miRTissueAtlasExpression.csv")

ggplot(data = tissue, aes(x = tissue, y=vsn, fill = miR))+geom_bar(stat = "identity", position = position_dodge())+ scale_fill_brewer(palette="Paired")
ggplot(data = tissue, aes(x = tissue, y=quantile.normalized, fill = miR))+geom_bar(stat = "identity", position = position_dodge())+ scale_fill_brewer(palette="Paired")
ggplot(data = tissue, aes(x = tissue, y=raw, fill = miR))+geom_bar(stat = "identity", position = position_dodge())+ scale_fill_brewer(palette="Paired")
