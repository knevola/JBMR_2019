# Figure 5
# Figure GO and KEGG Enrichment, merged GO and KEGG, separate by species 
rm(list = ls())
setwd("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/data")
library(VennDiagram)
library(dplyr)
library(ggplot2)

GO_19 <- read.csv("GO_miR19a.csv", stringsAsFactors = F)
KEGG_19<-read.csv("KEGG_miR19a.csv", stringsAsFactors = F)
GO_186<-read.csv("GO_miR186.csv")
KEGG_186 <-read.csv( "KEGG_miR186.csv")
GO_blue <- read.csv("BlueClusterGO.csv")
KEGG_blue <- read.csv("BlueClusterKEGG.csv")
GO_terms <- read.csv("GOterms.csv", stringsAsFactors = F)
GO_blue_all <-read.csv("BlueClusterunfilterGO.csv")
GO_19_all <- read.csv("GO_miR19a_unfiltered.csv")
GO_186_all <- read.csv("GO_miR186_unfiltered.csv")


GO_19_noont <- GO_19[-3]
names(KEGG_19)[2] <- "Term"
GOKEGG19 <- rbind(GO_19_noont, KEGG_19)
terms19a <- c("osteoblast differentiation","Insulin signaling pathway", "response to parathyroid hormone", "TGF-beta signaling pathway", "beta-adrenergic receptor kinase activity",
              "insulin-like growth factor receptor signaling pathway", "Cellular senescence", "Thyroid hormone signaling pathway")
GOKEGG19_f <- GOKEGG19 %>% filter(.,Term %in% terms19a)
GOKEGG19_f$log10p <- -log10(GOKEGG19_f$P.DE)

GOKEGG19_f <- GOKEGG19_f[order(GOKEGG19_f$log10p, decreasing = T),]
row.names(GOKEGG19_f) <- NULL
p19a<-ggplot(data = GOKEGG19_f, aes(x = reorder(Term, log10p), y = log10p)) + geom_bar(stat = "identity", position = position_dodge(), fill = "steelblue") + 
  geom_abline(slope = 0, intercept = -log10(0.05)) + coord_flip() + ylab(label = "-log10(p-value)") + xlab(label = "Term")+
  ggtitle(label = "miR-19a-3p") + ylim(0,5) + theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10), title = element_text(size = 10), legend.text = element_text(size = 10))

p19a

GO_186_noont <- GO_186[-3]
names(KEGG_186)[2] <- "Term"
GOKEGG186 <- rbind(GO_186_noont, KEGG_186)
terms186 <- c("TGF-beta signaling pathway", "Cellular senescence", "Thyroid hormone signaling pathway", "insulin metabolic process")
GOKEGG186_f <- GOKEGG186 %>% filter(.,Term %in% terms186)
GOKEGG186_f$log10p <- -log10(GOKEGG186_f$P.DE)

GOKEGG186_f <- GOKEGG186_f[order(GOKEGG186_f$log10p, decreasing = T),]
row.names(GOKEGG186_f) <- NULL
p186 <- ggplot(data = GOKEGG186_f, aes(x = reorder(Term, log10p), y = log10p)) + geom_bar(stat = "identity", position = position_dodge(), fill = "steelblue") + 
  geom_abline(slope = 0, intercept = -log10(0.05)) + coord_flip() + ylab(label = "-log10(p-value)") + xlab(label = "Term")+
  ggtitle(label = "miR-186-5p") + ylim (0,5) + theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10), title = element_text(size = 10), legend.text = element_text(size = 10))
p186


GO_blue_noont <- GO_blue[-3]
names(KEGG_blue)[2] <- "Term"
GOKEGGblue <- rbind(GO_blue_noont, KEGG_blue)
termsbluea <- c("osteoblast differentiation","Insulin signaling pathway", "bone development", "TGF-beta signaling pathway", "beta-adrenergic receptor kinase activity",
                "insulin-like growth factor receptor signaling pathway", "Cellular senescence", "Thyroid hormone signaling pathway", "osteoclast differentiation",
                "Estrogen signaling pathway", "beta-1 adrenergic receptor binding")
GOKEGGblue_f <- GOKEGGblue %>% filter(.,Term %in% termsbluea)
GOKEGGblue_f$log10p <- -log10(GOKEGGblue_f$P.DE)

GOKEGGblue_f <- GOKEGGblue_f[order(GOKEGGblue_f$log10p, decreasing = T),]
row.names(GOKEGGblue_f) <- NULL
p_blue <- ggplot(data = GOKEGGblue_f, aes(x = reorder(Term, log10p), y = log10p)) + geom_bar(stat = "identity", position = position_dodge(), fill = "steelblue", width = 0.9) + 
  geom_abline(slope = 0, intercept = -log10(0.05)) + coord_flip() + ylab(label = "-log10(p-value)") + xlab(label = "Term")+
  ggtitle(label = "Blue miRNA Module") + ylim(0,5) + theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10), title = element_text(size = 10), legend.text = element_text(size = 10))


library(ggpubr)
ggarrange(p19a, p186, p_blue, 
          labels = c("A", "B", "C"),
          ncol = 1, nrow = 3, align = "hv", heights = c(8,6,11))
ggsave(filename = "GOKeggCombined.tiff", width = 177, height = 150 ,unit = "mm")
