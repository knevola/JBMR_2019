# GSEA of Cor results
# mRNA correlation with miRNA of interest, edited 7/12/19
# Setup ####
rm(list = ls())
library(readr)
library(edgeR)
library(clusterProfiler)
library(enrichplot)
library(DOSE)
library(ggplot2)
library(org.Hs.eg.db)
setwd('/home/clary@mmcf.mehealth.org/Framingham/OmicData/data')
GPL5175 <- read_table2("GPL5175.txt", skip = 14)
cortarget19<-read.csv("miR19targetscor.csv")
cortarget186<-read.csv("miR186targetscor.csv")
hs <- org.Hs.eg.db
entrez_universe <- select(hs,keys = GPL5175$category, 
                          columns = c("ENTREZID", "SYMBOL"),
                          keytype = "SYMBOL")
universe_unique <- unique(entrez_universe$ENTREZID)

GO19<-goana(cortarget19$ENTREZID, universe = universe_unique)
KEGG19 <- kegga(cortarget19$ENTREZID,universe = universe_unique)
GO186<-goana(cortarget186$ENTREZID,universe = universe_unique)
KEGG186<-kegga(cortarget186$ENTREZID,universe = universe_unique)

GO19_p <- GO19[GO19$P.DE < 0.05,]
KEGG19_p <-KEGG19[KEGG19$P.DE< 0.05,]
GO186_p <- GO186[GO186$P.DE < 0.05,]
KEGG186_p <- KEGG186[KEGG186$P.DE < 0.05,]

write.csv(GO19, "GO_miR19a_unfiltered.csv")
write.csv(GO186, "GO_miR186_unfiltered.csv")
write.csv(GO19_p, "GO_miR19a.csv")
write.csv(KEGG19_p, "KEGG_miR19a.csv")
write.csv(GO186_p, "GO_miR186.csv")
write.csv(KEGG186_p, "KEGG_miR186.csv")

int_GO <- merge(GO19_p, GO186_p, by.x = "Term", by.y = "Term")
intersect(KEGG186_p$Pathway, KEGG19_p$Pathway)


edo_19 <- enrichDGN(na.omit(cortarget19$ENTREZID),pvalueCutoff = 0.1,pAdjustMethod = "BH",universe = universe_unique)
library(enrichplot)


edo_186 <- enrichDGN(na.omit(cortarget186$ENTREZID),pvalueCutoff = 0.1,pAdjustMethod = "BH",universe = universe_unique)

library(dplyr)
library(multiMiR)
library(edgeR)
library(clusterProfiler)
setwd("/home/clary@mmcf.mehealth.org/Framingham/OmicData/Clustering_miRNA")
WGCNA <- read.csv("MiRNA_significance_module_membership_nofilter.csv")
WGCNA.blue <- WGCNA %>% filter(., mergedColors == "blue")
blue_miRNA <- as.character(WGCNA.blue$X)
blue_miRNA <- gsub(x = blue_miRNA,pattern ="_",replacement = "-" )
blue_miRNA <- gsub(x = blue_miRNA, pattern = "miR", replacement = "hsa-miR")
blue_miRNA <- gsub(x = blue_miRNA, pattern = "5p-a2", replacement = "5p")
blue_miRNA <- gsub(x = blue_miRNA, pattern = "5p-a1", replacement = "5p")
genes<- get_multimir(mirna = blue_miRNA, summary = TRUE)
genesdf <- genes@data
setdiff(blue_miRNA,unique(genesdf$mature_mirna_id))
genesdfu <- genesdf[3:5]
genesdfuni <- unique(genesdfu)

GO<-goana(unique(genesdf$target_entrez))
KEGG <- kegga(unique(genesdf$target_entrez))

GO_p <- GO[GO$P.DE < 0.05,]
KEGG_p <- KEGG[KEGG$P.DE < 0.05,]

write.csv(GO_p, "BlueModuleGO_filtered.csv")
write.csv(KEGG_p, "BlueModuleKEGG_filtered.csv")
length(unique(genesdf$target_entrez))

table<-as.data.frame(table(genesdfuni$target_entrez))
table <- table[order(table$Freq),]
barplot(table$Freq )

write.csv(table, "FrequencyofGenes.csv")
quantile(table$Freq, probs = seq(0, 1, 0.1))
table2<-table[table$Freq > 2 & table$Freq < 15,]
tablesym <- as.data.frame(table(genesdfuni$target_symbol))
tablesym2 <- tablesym[tablesym$Freq > 2,]

GO2 <- goana(as.vector(table2$Var1))
KEGG2 <- kegga(as.vector(table2$Var1))

GO2_p <- GO2[GO2$P.DE < 0.05,]
KEGG2_p <- KEGG2[KEGG2$P.DE < 0.05,]

setwd("/home/clary@mmcf.mehealth.org/Framingham/OmicData/data")
GO186 <- read.csv("GO_miR186.csv")
GO19 <- read.csv("GO_miR19a.csv")
KEGG186 <- read.csv("KEGG_miR186.csv")
KEGG19 <- read.csv("KEGG_miR19a.csv")

intersect(intersect(GO186$Term, GO19$Term), GO2_p$Term)

library(enrichplot)
library(DOSE)
edo_WGCNA <- enrichDGN(table2$Var1)
barplot(edo, showCategory=30)
dotplot(edo, showCategory=10) + ggtitle("dotplot for ORA")

library(dplyr)
WGCNA_result<-edo_WGCNA@result
WGCNA_result1 <- WGCNA_result %>% filter(.,p.adjust < 0.1)
Result19 <- edo_19@result
Result191 <- Result19 %>%  filter(., p.adjust < 0.1)
Result186 <- edo_186@result
Result1861 <- Result186 %>%  filter(., p.adjust < 0.1)



int<-intersect(intersect(WGCNA_result1$Description, Result191$Description), Result1861$Description)

int19 <- Result191[Result191$Description %in% int,]
int19$Group <- "miR-19a-3p"
int186 <- Result1861[Result1861$Description %in% int,]
int186$Group <- "miR-186-5p"
intWGCNA <- WGCNA_result1[WGCNA_result1$Description %in% int,]
intWGCNA$Group <- "Blue Cluster"
intall<-rbind(int19, int186, intWGCNA)
intall$log_padj <- -log(intall$p.adjust)

ggplot(data = intall, aes(x = Description, y = log_padj, fill = Group)) + geom_bar(stat = "identity", position = position_dodge())+
  theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(y = "-log(FDR)") + scale_fill_brewer(palette = "Paired") 

IntGO <- intersect(intersect(GO19_p$Term, GO186_p$Term), GO2_p$Term)
IntKEGG <- intersect(intersect(KEGG19_p$Pathway, KEGG186_p$Pathway), KEGG2_p$Pathway)

