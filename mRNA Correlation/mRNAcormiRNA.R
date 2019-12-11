# mRNA correlation with miRNA of interest, edited 7/12/19
# Setup ####
rm(list = ls())
library(tidyverse)
library(VennDiagram)
library(multiMiR)
setwd('/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/data')
miRNAdat <- read.csv('l_mrna_2011_m_0797s_17_c1.csv')
pheno <- read.csv('PhenoData_5_28.csv')
mRNA_sample <- read.delim("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/data/phe000002.v7_release_manifest.txt", comment.char="#")
mRNAdat <- read_delim("FinalFile_Gene_OFF_2446_Adjusted_c1.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
pheno_samp <- merge(pheno, mRNA_sample[,c(1,2)], by.x = "shareid", by.y = "Subject_ID")
int <- read.csv("Spineoverlapresults_7_9_rank.csv")
GPL5175 <- read_table2("GPL5175.txt", skip = 14)

miRNA_delta_cq <- miRNAdat[-1]
miRNA_delta_cq <- -(miRNA_delta_cq-27)
miRNAdat <- cbind(miRNAdat[1], miRNA_delta_cq)

data_adj_t <-as.data.frame(t(mRNAdat))
data_adj_t <- data_adj_t[-1,]
colnames(data_adj_t)<- mRNAdat$transcript_cluster_id

dat<- merge(pheno_samp, data_adj_t, by.x = "Sample_ID", by.y = 0)
dat$ADRB1 <- dat$`3265140`
dat$ADRB2 <- dat$`2834743`
dat$HDAC4 <- dat$`2606026`
dat$RUNX2 <- dat$`2908762`
dat$CTHRC1 <- dat$`3110317`
dat$APOE <- dat$`3835879`
hist(dat$CTHRC1)
hist(dat$ADRB2)
hist(dat$APOE)
summary(dat$ADRB1)
summary(dat$ADRB2)
keep <- c("shareid",as.character(int$miRNA))
miRNA_subset <- miRNAdat[,names(miRNAdat) %in%keep]

mmdat <- merge(dat, miRNA_subset, by.x = "shareid", by.y = "shareid")
mRNAs <- names(dat)[-c(1:23)]
cor.test(~miR_19a_3p + ADRB1,data = mmdat, method = "spearman")
cor.test(~miR_19a_3p + ADRB1,data = mmdat)
cor.test(~miR_19a_3p + ADRB2,data = mmdat, method = "spearman")
cor.test(~miR_19a_3p + HDAC4,data = mmdat, method = "spearman")
cor.test(~miR_19a_3p + RUNX2,data = mmdat, method = "spearman")

plot(mmdat$miR_19a_3p, mmdat$ADRB1)
 x<-lm(ADRB1~miR_19a_3p+AGE8 + SEX+HGT8 +WGT8, data = mmdat)
 summary(x)
 hist(mmdat$ADRB1)
 ggplot(data = mmdat, aes(x = miR_19a_3p, y = ADRB1)) + geom_smooth(method = lm, se = TRUE, fullrange = TRUE) + ggtitle(label =  "miR-19a-3p by Total Femur BMD") + 
   xlab("miR-19a") + ylab ("ADRB1") + theme_minimal() + scale_color_brewer(palette = "Dark2") + theme(plot.title = element_text(hjust = 0.5)) +theme(axis.text = element_text(size = 6), axis.title = element_text(size = 8), title = element_text(size = 10), legend.text = element_text(size = 8), legend.position = c(0.1,0.8), legend.background = element_rect(fill="white",colour = "white")) 
 

cormir<-function(mRNAs, miRNA, data)
{
  for(i in 1:length(mRNAs)){
    print(i)
    x <- cor.test(as.formula(paste("~", miRNA, " + `", mRNAs[i],"`", sep = ""  )), data = data, method = "s")
    estimate<- x$estimate
    p_value <- x$p.value
    if (i == 1){
      df <- data.frame(Estimate = estimate, Pvalue = p_value)
    }else{
      df2 <- data.frame(Estimate = estimate, Pvalue = p_value)
      df <- rbind(df,df2)
    }
  }
  row.names(df)<- mRNAs
  df$fdr <- p.adjust(df$Pvalue)
  return(df)
}
# mRNAs correlated with miRNA of interest (paper uses lmekin need to research more)
cor19 <- cormir(mRNAs = mRNAs, miRNA = "miR_19a_3p", data = mmdat)
cor19_fdr <- cor19[cor19$Pvalue< 0.05,]
cor186 <- cormir(mRNAs = mRNAs, miRNA = "miR_186_5p_a2", data = mmdat)
cor186_fdr <- cor186[cor186$Pvalue < 0.05,]

cor19_f_pos <- cor19_fdr[cor19_fdr$Estimate < 0,]
cor186_f_pos <- cor186_fdr[cor186_fdr$Estimate < 0,]

library(readr)
GPL5175 <- read_table2("GPL5175.txt", skip = 14)
cor19_geneinfo <- merge(cor19_f_pos, GPL5175, by.x = 0, by.y = 1)
cor186_geneinfo <- merge(cor186_f_pos, GPL5175, by.x = 0, by.y =1)

cor_19_sym <- cor19_geneinfo$category
cor_186_sym <- cor186_geneinfo$category
intersect(cor_19_sym, cor_186_sym)
library(org.Hs.eg.db)
hs <- org.Hs.eg.db
entrez19 <- select(hs,keys = cor_19_sym, 
                         columns = c("ENTREZID", "SYMBOL"),
                         keytype = "SYMBOL")
entrez186 <- select(hs,keys = cor_186_sym, 
                    columns = c("ENTREZID", "SYMBOL"),
                    keytype = "SYMBOL")
entrez_universe <- select(hs,keys = GPL5175$category, 
                    columns = c("ENTREZID", "SYMBOL"),
                    keytype = "SYMBOL")
universe_unique <- unique(entrez_universe$ENTREZID)
# MiRNA targets
targets19<-get_multimir(mirna = 'hsa-miR-19a-3p', summary = TRUE)
targets186 <- get_multimir(mirna = 'hsa-miR-186-5p', summary = T)

df19<-targets19@data
df186 <- targets186@data

int19 <- intersect(df19$target_symbol, entrez19$SYMBOL)
int186 <- intersect(df186$target_symbol, entrez186$SYMBOL)

cortarget19 <- cor19_geneinfo[cor19_geneinfo$category %in% int19,]
cortarget19 <- merge(cortarget19, entrez19, by.x = "category", by.y = "SYMBOL")
cortarget186 <- cor186_geneinfo[cor186_geneinfo$category %in% int186,]
cortarget186 <- merge(cortarget186, entrez186, by.x = "category", by.y = "SYMBOL")
write.csv(cortarget19,"miR19targetscor.csv")
write.csv(cortarget186, "miR186targetscor.csv")

cortarget19 <-read.csv("miR19targetscor.csv")
cortarget186 <- read.csv("miR186targetscor.csv")
library(edgeR)
GO19<-goana(de = na.omit(cortarget19$ENTREZID), universe = universe_unique)
KEGG19 <- kegga(de =cortarget19$ENTREZID,universe = universe_unique)
GO186<-goana(cortarget186$ENTREZID,universe = universe_unique)
KEGG186<-kegga(cortarget186$ENTREZID,universe = universe_unique)

GO19_p <- GO19[GO19$P.DE < 0.05,]
KEGG19_p <-KEGG19[KEGG19$P.DE< 0.05,]
GO186_p <- GO186[GO186$P.DE < 0.05,]
KEGG186_p <- KEGG186[KEGG186$P.DE < 0.05,]

write.csv(GO19_p, "GO_miR19a.csv")
write.csv(KEGG19_p, "KEGG_miR19a.csv")
write.csv(GO186_p, "GO_miR186.csv")
write.csv(KEGG186_p, "KEGG_miR186.csv")

int_GO <- merge(GO19_p, GO186_p, by.x = "Term", by.y = "Term")
intersect(KEGG186_p$Pathway, KEGG19_p$Pathway)

library(biomaRt)
ensembl = useMart(biomart = "ensembl",dataset="hsapiens_gene_ensembl") #uses human ensembl annotations
#gets gene symbol, transcript_id and go_id for all genes annotated with GO:0007507
gene.data <- getBM(attributes=c('ensembl_gene_id','entrezgene_id','hgnc_id','hgnc_symbol'),
                   filters = 'go', values = 'GO:0047696', mart = ensembl)
insulin <- getBM(attributes=c('ensembl_gene_id','entrezgene_id','hgnc_id','hgnc_symbol'),
                 filters = 'go', values = 'GO:0032868', mart = ensembl)
insulin186 <- getBM(attributes=c('ensembl_gene_id','entrezgene_id','hgnc_id','hgnc_symbol'),
                 filters = 'go', values = 'GO:1901142', mart = ensembl)
PTH <- getBM(attributes=c('ensembl_gene_id','entrezgene_id','hgnc_id','hgnc_symbol'),
             filters = 'go', values = 'GO:0071107', mart = ensembl)
  
intersect(PTH$entrezgene_id, cortarget186$ENTREZID)

GK <- getGeneKEGGLinks(species.KEGG = "hsa")
#TGF-B
cortarget19$category[cortarget19$ENTREZID %in% GK$GeneID[GK$PathwayID=="path:hsa04350"]]
cortarget186$category[cortarget186$ENTREZID %in% GK$GeneID[GK$PathwayID=="path:hsa04350"]]
#PTH
cortarget19$category[cortarget19$ENTREZID %in% GK$GeneID[GK$PathwayID=="path:hsa04928"]]
cortarget186$category[cortarget186$ENTREZID %in% GK$GeneID[GK$PathwayID=="path:hsa04928"]]
#Insulin
cortarget19$category[cortarget19$ENTREZID %in% GK$GeneID[GK$PathwayID=="path:hsa04910"]]
cortarget186$category[cortarget186$ENTREZID %in% GK$GeneID[GK$PathwayID=="path:hsa04910"]]

#Estrogen
cortarget19$category[cortarget19$ENTREZID %in% GK$GeneID[GK$PathwayID=="path:hsa04915"]]
cortarget186$category[cortarget186$ENTREZID %in% GK$GeneID[GK$PathwayID=="path:hsa04915"]]
