# Enrichment of all clusters mRNA and miRNA
rm(list = ls())
setwd("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MergingClustering")
library(dplyr)
library(multiMiR)
library(edgeR)
library(clusterProfiler)
library(readr)
library(org.Hs.eg.db)
hs <- org.Hs.eg.db
GPL5175 <- read_table2("/home/clary@mmcf.mehealth.org/Framingham/OmicData/data/GPL5175.txt", skip = 14)
miRNA <- read.csv("/home/clary@mmcf.mehealth.org/Framingham/OmicData/Clustering_miRNA/MiRNA_significance_module_membership_nofilter.csv")
mRNA <- read.csv("/home/clary@mmcf.mehealth.org/Framingham/OmicData/Clustering_mRNA/mRNA_geneModuleMembership.csv")

miRtarget <- function(miRNAs){
  miRNAs <- gsub(x = miRNAs,pattern ="_",replacement = "-" )
  miRNAs <- gsub(x = miRNAs, pattern = "miR", replacement = "hsa-miR")
  miRNAs <- gsub(x = miRNAs, pattern = "5p-a2", replacement = "5p")
  miRNAs <- gsub(x = miRNAs, pattern = "5p-a1", replacement = "5p")
  genes<- get_multimir(mirna = miRNAs, summary = TRUE)
  genes <- genes@data
  return(genes)
}

miRfilter <- function(genes){
  table <- as.data.frame(table(genes$target_entrez))
  table2 <- table[table$Freq > 2,]
  genes2 <- genes %>% filter(., target_entrez %in% table2$Var1)
}

'%!in%' <- function(x,y)!('%in%'(x,y))

GOmRNA <- function(mRNAlist, universe){
  GO <- goana(de = mRNAlist, universe = universe, FDR = 0.05)
  GO_p <- GO[GO$P.DE < 0.001,]
}

KEGGmRNA <- function(mRNAlist, universe){
  KEGG <- kegga(de = mRNAlist, universe = universe, FDR = 0.05)
  KEGG_p <- KEGG[KEGG$P.DE < 0.001,]
}

# MIR ENRICHMENT ####
# Separated in miRNA clusters
table(miRNA$mergedColors)
miRNA_blue <- miRNA %>% filter(., mergedColors == "blue")
miRNA_brown <- miRNA %>% filter(., mergedColors == "brown")
miRNA_green <- miRNA %>% filter(., mergedColors == "green")
miRNA_turquoise <- miRNA %>% filter(., mergedColors == "turquoise")
miRNA_yellow <- miRNA %>% filter(., mergedColors == "yellow")

# Gene Targets of Each Cluster
genetargets_blue <- miRtarget(miRNAs = as.character(miRNA_blue$miRNA))
genetargets_brown <- miRtarget(miRNAs = as.character(miRNA_brown$miRNA))
genetargets_green <- miRtarget(miRNAs = as.character(miRNA_green$miRNA))
genetargets_turquoise <- miRtarget(miRNAs = as.character(miRNA_turquoise$miRNA))
genetargets_yellow <- miRtarget(miRNAs = as.character(miRNA_yellow$miRNA))

# Gene Filter
genetargets_blue2 <- miRfilter(genetargets_blue)
genetargets_brown2 <- miRfilter(genetargets_brown)
genetargets_green2 <- miRfilter(genetargets_green)
genetargets_turquoise2 <- miRfilter(genetargets_turquoise)
genetargets_yellow2 <- miRfilter(genetargets_yellow)


#GO miRs
GO_blue <- goana(genetargets_blue2$target_entrez, FDR = 0.05)
GO_brown <- goana(genetargets_brown2$target_entrez, FDR = 0.05)
GO_green <- goana(genetargets_green2$target_entrez, FDR = 0.05)
GO_turquoise <- goana(genetargets_turquoise2$target_entrez, FDR = 0.05)
GO_yellow <- goana(genetargets_yellow2$target_entrez, FDR = 0.05)

# Filter GO
GO_blue_p <- GO_blue[GO_blue$P.DE < 0.05,]
GO_brown_p <- GO_brown[GO_brown$P.DE < 0.05,]
GO_green_p <- GO_green[GO_green$P.DE < 0.05,]
GO_turquoise_p <- GO_turquoise[GO_turquoise$P.DE < 0.05,]
GO_yellow_p <- GO_yellow[GO_yellow$P.DE < 0.05,]

GO_blue_only<-GO_blue_p[GO_blue_p$Term %!in% c(GO_brown_p$Term, GO_green_p$Term, GO_turquoise_p$Term, GO_yellow_p$Term),]
GO_brown_only <- GO_brown_p[GO_brown_p$Term %!in% c(GO_blue_p$Term, GO_green_p$Term, GO_turquoise_p$Term, GO_yellow_p$Term),]
GO_green_only <- GO_green_p[GO_green_p$Term %!in% c(GO_blue_p$Term, GO_brown_p$Term, GO_turquoise_p$Term, GO_yellow_p$Term),]
GO_turquoise_only <- GO_turquoise_p[GO_turquoise_p$Term %!in% c(GO_blue_p$Term, GO_brown_p$Term, GO_green_p$Term, GO_yellow_p$Term),]
GO_yellow_only <- GO_yellow_p[GO_yellow_p$Term %!in% c(GO_blue_p$Term, GO_brown_p$Term, GO_green_p$Term, GO_turquoise_p$Term),]

# KEGG miRs 
KEGG_blue <- kegga(genetargets_blue2$target_entrez)
KEGG_brown <- kegga(genetargets_brown2$target_entrez)
KEGG_green <- kegga(genetargets_green2$target_entrez)
KEGG_turquoise <- kegga(genetargets_turquoise2$target_entrez)
KEGG_yellow <- kegga(genetargets_yellow2$target_entrez)

# Filter KEGG
KEGG_blue_p <- KEGG_blue[KEGG_blue$P.DE < 0.05,]
KEGG_brown_p <- KEGG_brown[KEGG_brown$P.DE < 0.05,]
KEGG_green_p <- KEGG_green[KEGG_green$P.DE < 0.05,]
KEGG_turquoise_p <- KEGG_turquoise[KEGG_turquoise$P.DE < 0.05,]
KEGG_yellow_p <- KEGG_yellow[KEGG_yellow$P.DE < 0.05,]

KEGG_blue_only <- KEGG_blue_p[KEGG_blue_p$Pathway %!in% c(KEGG_brown_p$Pathway, KEGG_green_p$Pathway, KEGG_turquoise_p$Pathway, KEGG_yellow_p$Pathway),]
KEGG_brown_only <- KEGG_brown_p[KEGG_brown_p$Pathway %!in% c(KEGG_blue_p$Pathway, KEGG_green_p$Pathway, KEGG_turquoise_p$Pathway, KEGG_yellow_p$Pathway),]
KEGG_green_only <- KEGG_green_p[KEGG_green_p$Pathway %!in% c(KEGG_brown_p$Pathway, KEGG_blue_p$Pathway, KEGG_turquoise_p$Pathway, KEGG_yellow_p$Pathway),]
KEGG_turquoise_only <- KEGG_turquoise_p[KEGG_turquoise_p$Pathway %!in% c(KEGG_brown_p$Pathway, KEGG_green_p$Pathway, KEGG_blue_p$Pathway, KEGG_yellow_p$Pathway),]
KEGG_yellow_only <- KEGG_yellow_p[KEGG_yellow_p$Pathway %!in% c(KEGG_brown_p$Pathway, KEGG_green_p$Pathway, KEGG_turquoise_p$Pathway, KEGG_blue_p$Pathway),]


# GENE ENRICHMENT
table(mRNA$NA.)
mRNA_GPL<-merge(mRNA, GPL5175, by.x = "X", by.y = "ID")
entrez_universe <- select(hs,keys = GPL5175$category, 
                          columns = c("ENTREZID", "SYMBOL"),
                          keytype = "SYMBOL")
mRNA_entrez <- merge(mRNA_GPL, entrez_universe, by.x = "category", by.y = "SYMBOL")

# Clusters
mRNA_black <- mRNA_entrez %>%  filter(.,NA. == "black")
mRNA_blue <- mRNA_entrez %>% filter(., NA. == "blue")
mRNA_brown<- mRNA_entrez %>% filter(., NA. == "brown")
mRNA_green <- mRNA_entrez %>% filter(., NA. == "green")
mRNA_greenyellow <- mRNA_entrez %>% filter(., NA. == "greenyellow")
mRNA_magenta <- mRNA_entrez %>% filter(., NA. == "magenta")
mRNA_pink <- mRNA_entrez %>% filter(., NA. == "pink")
mRNA_purple <- mRNA_entrez %>% filter(., NA. == "purple")
mRNA_red <- mRNA_entrez %>% filter(., NA. == "red")
mRNA_turquoise <- mRNA_entrez %>% filter(., NA. == "turquoise")
mRNA_yellow <- mRNA_entrez %>% filter(., NA. == "yellow")

# GO
GO_black_m_p <- GOmRNA(mRNAlist = mRNA_black$ENTREZID, universe = entrez_universe$ENTREZID)
GO_blue_m_p <- GOmRNA(mRNAlist = mRNA_blue$ENTREZID, universe = entrez_universe$ENTREZID)
GO_brown_m_p <- GOmRNA(mRNAlist = mRNA_brown$ENTREZID, universe = entrez_universe$ENTREZID)
GO_green_m_p <- GOmRNA(mRNAlist = mRNA_green$ENTREZID, universe = entrez_universe$ENTREZID)
GO_greenyellow_m_p <- GOmRNA(mRNAlist = mRNA_greenyellow$ENTREZID, universe = entrez_universe$ENTREZID)
GO_magenta_m_p <- GOmRNA(mRNAlist = mRNA_magenta$ENTREZID, universe = entrez_universe$ENTREZID)
GO_pink_m_p <- GOmRNA(mRNAlist = mRNA_pink$ENTREZID, universe = entrez_universe$ENTREZID)
GO_purple_m_p <- GOmRNA(mRNAlist = mRNA_purple$ENTREZID, universe = entrez_universe$ENTREZID)
GO_red_m_p <- GOmRNA(mRNAlist = mRNA_red$ENTREZID, universe = entrez_universe$ENTREZID)
GO_turquoise_m_p <- GOmRNA(mRNAlist = mRNA_turquoise$ENTREZID, universe = entrez_universe$ENTREZID)
GO_yellow_m_p <- GOmRNA(mRNAlist = mRNA_yellow$ENTREZID, universe = entrez_universe$ENTREZID)

KEGG_bluem <- kegga(mRNA_blue$ENTREZID, universe = entrez_universe$ENTREZID)
KEGG_bluem_p <- KEGG_bluem[KEGG_bluem$P.DE < 0.05,]
KEGG_brownm <- kegga(mRNA_brown$ENTREZID, universe = entrez_universe$ENTREZID)
KEGG_brownm_p <- KEGG_brownm[KEGG_brownm$P.DE < 0.05,]
KEGG_turquoisem <- kegga(mRNA_turquoise$ENTREZID, universe = entrez_universe$ENTREZID)
KEGG_turquoisem_p <- KEGG_turquoisem[KEGG_turquoisem$P.DE < 0.05,]

GO_bluem <- goana(mRNA_blue$ENTREZID, universe = entrez_universe$ENTREZID)
GO_bluem_p <- GO_bluem[GO_bluem$P.DE < 0.05,]
GO_brownm <- goana(mRNA_brown$ENTREZID, universe = entrez_universe$ENTREZID)
GO_brownm_p <- GO_brownm[GO_brownm$P.DE < 0.05,]
GO_turquoisem <- goana(mRNA_turquoise$ENTREZID, universe = entrez_universe$ENTREZID)
GO_turquoisem_p <- GO_turquoisem[GO_turquoisem$P.DE < 0.05,]
# Filter GO
GO_black_m_only <- GO_black_m_p[GO_black_m_p$Term %!in% c(GO_blue_m_p$Term, GO_brown_m_p$Term, GO_green_m_p$Term, GO_greenyellow_m_p$Term, GO_magenta_m_p$Term, GO_pink_m_p$Term, GO_purple_m_p$Term, GO_red_m_p$Term, GO_turquoise_m_p$Term, GO_yellow_m_p$Term), ]
GO_blue_m_only <- GO_blue_m_p[GO_blue_m_p$Term %!in% c(GO_black_m_p$Term, GO_brown_m_p$Term, GO_green_m_p$Term, GO_greenyellow_m_p$Term, GO_magenta_m_p$Term, GO_pink_m_p$Term, GO_purple_m_p$Term, GO_red_m_p$Term, GO_turquoise_m_p$Term, GO_yellow_m_p$Term), ]
GO_brown_m_only <- GO_brown_m_p[GO_brown_m_p$Term %!in% c(GO_black_m_p$Term, GO_blue_m_p$Term, GO_green_m_p$Term, GO_greenyellow_m_p$Term, GO_magenta_m_p$Term, GO_pink_m_p$Term, GO_purple_m_p$Term, GO_red_m_p$Term, GO_turquoise_m_p$Term, GO_yellow_m_p$Term), ]
GO_green_m_only <- GO_green_m_p[GO_green_m_p$Term %!in% c(GO_black_m_p$Term, GO_brown_m_p$Term, GO_blue_m_p$Term, GO_greenyellow_m_p$Term, GO_magenta_m_p$Term, GO_pink_m_p$Term, GO_purple_m_p$Term, GO_red_m_p$Term, GO_turquoise_m_p$Term, GO_yellow_m_p$Term), ]
GO_greenyellow_m_only <- GO_greenyellow_m_p[GO_greenyellow_m_p$Term %!in% c(GO_black_m_p$Term, GO_brown_m_p$Term, GO_blue_m_p$Term, GO_green_m_p$Term, GO_magenta_m_p$Term, GO_pink_m_p$Term, GO_purple_m_p$Term, GO_red_m_p$Term, GO_turquoise_m_p$Term, GO_yellow_m_p$Term), ]
GO_magenta_m_only <- GO_magenta_m_p[GO_magenta_m_p$Term %!in% c(GO_black_m_p$Term, GO_brown_m_p$Term, GO_blue_m_p$Term, GO_greenyellow_m_p$Term, GO_green_m_p$Term, GO_pink_m_p$Term, GO_purple_m_p$Term, GO_red_m_p$Term, GO_turquoise_m_p$Term, GO_yellow_m_p$Term), ]
GO_pink_m_only <- GO_pink_m_p[GO_pink_m_p$Term %!in% c(GO_black_m_p$Term, GO_brown_m_p$Term, GO_blue_m_p$Term, GO_greenyellow_m_p$Term, GO_magenta_m_p$Term, GO_green_m_p$Term, GO_purple_m_p$Term, GO_red_m_p$Term, GO_turquoise_m_p$Term, GO_yellow_m_p$Term), ]
GO_purple_m_only <- GO_purple_m_p[GO_purple_m_p$Term %!in% c(GO_black_m_p$Term, GO_brown_m_p$Term, GO_blue_m_p$Term, GO_greenyellow_m_p$Term, GO_magenta_m_p$Term, GO_pink_m_p$Term, GO_green_m_p$Term, GO_red_m_p$Term, GO_turquoise_m_p$Term, GO_yellow_m_p$Term), ]
GO_red_m_only <- GO_red_m_p[GO_red_m_p$Term %!in% c(GO_black_m_p$Term, GO_brown_m_p$Term, GO_blue_m_p$Term, GO_greenyellow_m_p$Term, GO_magenta_m_p$Term, GO_pink_m_p$Term, GO_green_m_p$Term, GO_purple_m_p$Term, GO_turquoise_m_p$Term, GO_yellow_m_p$Term), ]
GO_turquoise_m_only <- GO_turquoise_m_p[GO_turquoise_m_p$Term %!in% c(GO_black_m_p$Term, GO_brown_m_p$Term, GO_blue_m_p$Term, GO_greenyellow_m_p$Term, GO_magenta_m_p$Term, GO_pink_m_p$Term, GO_green_m_p$Term, GO_red_m_p$Term, GO_purple_m_p$Term, GO_yellow_m_p$Term), ]
GO_yellow_m_only <- GO_yellow_m_p[GO_yellow_m_p$Term %!in% c(GO_black_m_p$Term, GO_brown_m_p$Term, GO_blue_m_p$Term, GO_greenyellow_m_p$Term, GO_magenta_m_p$Term, GO_pink_m_p$Term, GO_green_m_p$Term, GO_red_m_p$Term, GO_turquoise_m_p$Term, GO_purple_m_p$Term), ]

intersect(GO_turquoisem_p$Term,intersect(GO_bluem_p$Term, GO_brownm_p$Term))
