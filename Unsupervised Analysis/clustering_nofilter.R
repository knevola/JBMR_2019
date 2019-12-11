# WGCNA no filtering
rm(list = ls())
library(tidyverse)
library(emmeans)
library(limma)
library(WGCNA)
library(tidyverse)
library(VennDiagram)
options(stringsAsFactors = T)
setwd('/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/data')
miRNA <- read.csv('l_mrna_2011_m_0797s_17_c1.csv')
pheno <- read.csv('PhenoData_5_28.csv')

miRNA_delta_cq <- miRNA[-1]
miRNA_delta_cq <- -(miRNA_delta_cq-27)
miRNA <- cbind(miRNA[1], miRNA_delta_cq)

# Subset analytic dataset/phenotype data ####
pheno3 <- pheno[is.na(pheno$f8cbnbmd) == F,]

# Merge pheno3 with miRNA ####
miRNA_pheno <- merge(pheno3, miRNA, by.x = 1, by.y = 1)
drop <- c('cvdpair', 'casecontrol', 'idtype')
miRNA_pheno <- miRNA_pheno[, !(names(miRNA_pheno) %in% drop)]
dataExpr <- miRNA_pheno[,-c(1:22)]

gag <- goodSamplesGenes(dataExpr, verbose = 0)
gag$allOK # False so  need to remove any offending genes
if(!gag$allOK){
  dataExpr = dataExpr[gag$goodSamples, gag$goodGenes]
}

# Sample Clustering -------------------------------------------------------

sampleTree = hclust(dist(dataExpr), method = "average")


# Graph Sample Clustering
sizeGrWindow(12,9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample Clustering to detect outliers", sub = "",
     xlab = "", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

# Look at clustering by traits
datTraits<-miRNA_pheno[c(1:22)]
datTraits1<- datTraits[gag$goodSamples,]
datTraits1$SEX <- as.numeric(datTraits1$SEX)
datTraits1$CURRSMK8 <- as.numeric(datTraits1$CURRSMK8)
datTraits1$DMRX8 <- as.numeric(datTraits1$DMRX8)
datTraits1$HRX8 <- as.numeric(datTraits1$HRX8)
datTraits1$LIPRX8 <- as.numeric(datTraits1$LIPRX8)
datTraits1$EST8 <- as.numeric(datTraits1$EST8)
datTraits1$BB <- as.numeric(datTraits1$BB)
datTraits1$B1 <- as.numeric(datTraits1$B1)
datTraits1$priorcvd <- as.numeric(datTraits1$priorcvd)
datTraits1$menov <- as.numeric(datTraits1$menov)
is.na(datTraits1$menov)<-2

traitColors = numbers2colors(datTraits1, signed = FALSE)

plotDendroAndColors(sampleTree, traitColors, groupLabels = names(datTraits1),
                    main = "Sample dendrogram and trait heatmap")
# No clear clustering of samples by traits

save(dataExpr, datTraits1, file = "FHS_miRNA_1.RData")



# Variance Filtering ------------------------------------------------------
thevar <- diag(var(dataExpr, na.rm = T))
vars<- quantile(thevar, na.rm = T)
setwd("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/Clustering_miRNA/Figures")
write.csv(vars, "Var_miRNA_nofilter.csv")


# module detection --------------------------------------------------------
#Automatic network construction and module detection
#Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
#Call the network topology analysis function
sft = pickSoftThreshold(dataExpr, powerVector = powers, verbose = 5)
#Plot the results
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9 
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2", type = "n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels = powers, cex = cex1, col = "blue")
#This line corresponds to used an R^2 cutoff of h
abline(h=0.9, col = "red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="blue")

# Choosing a power of 8 based on mean graphs

#Constructing the gene network and identifying modules: auto
net = blockwiseModules(dataExpr, power = 8,
                       TOMType = "unsigned", minModuleSize = 10,
                       reassignThreshold = 0, mergeCutHeight = 0.05,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "femaleMouseTOM",
                       verbose = 3)
#How many modules were identified? (6/5)
modulebkdn<- table(net$colors)
write.csv(modulebkdn, "Genes_per_module_nofilter.csv")

#View heirarchical clustering
#Open a graphics window
sizeGrWindow(12,9)
#convert labels to colors for ploting
mergedColors = labels2colors(net$colors)
#Plot the dendrogram
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors", dendroLabels = F, hang = 0.03,
                    addGuide = T, guideHang = 0.05)
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "miRNA-networkConstruction-auto_nofilter.RData")

# Relate Modules to Traits ----------------------------------------
#Recalculate MEs with color labels
nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)

MEs0 = moduleEigengenes(dataExpr, moduleColors)$eigengenes
row.names(MEs0)<- datTraits1$shareid
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits1, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

#Graphical Representation
sizeGrWindow(10,6)
#Display correlations and their Pvalue
textMatrix = paste(signif(moduleTraitCor,2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3,3))

#Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor, 
               xLabels = names(datTraits1),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = F,
               cex.text = 0.5, 
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

#Names of modules (colors)
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(dataExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep = "")
names(MMPvalue) = paste("p.MM", modNames, sep = "")


write.csv(MEs, "miRNA_MEs_nofilter.csv")
write.csv(geneModuleMembership, "miRNA_geneModuleMembership_nofilter.csv")

genes_colors <- cbind(geneModuleMembership,mergedColors)
miRNA <- row.names(genes_colors)
genes_colors<- cbind(genes_colors, miRNA)
write.csv(genes_colors, "MiRNA_significance_module_membership_nofilter.csv")
write.csv(moduleTraitCor, "miRNA_Module_trait_cor_nofilter.csv")
write.csv(moduleTraitPvalue, "miRnA_Module_trait_pvalue_nofilter.csv")

# brown and blue are sig with BMD and BB
geneofinterest <- genes_colors %>%  filter(., mergedColors %in% c("blue", "brown"))
#Connection between blue and brown clusters
dissTOM <- 1- TOMsimilarityFromExpr(datExpr = dataExpr, power = 8)
plotTOM <- dissTOM^9
diag(plotTOM) =NA
sizeGrWindow(9,9)
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")

FNBMD <- datTraits1$f8cbtobmd
MET = orderMEs(cbind(MEs, FNBMD))
sizeGrWindow(5,7.5);
par(cex = 0.9)
plotEigengeneNetworks(MEs, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90)
write.csv(geneofinterest, "GenesofInterest_nofilter.csv")

# Recalculate topological overlap if needed
TOM = TOMsimilarityFromExpr(dataExpr, power = 8);
# Read in the annotation file
# Select modules
modules = c("brown", "blue");
modules2 = c("brown", "blue", "green", "turquoise", "yellow", "grey");
modules3 = "blue"
# Select module probes
probes = names(dataExpr)
inModule = is.finite(match(moduleColors, modules));
inModule2 = is.finite(match(moduleColors, modules2));
inModule3 =is.finite(match(moduleColors, modules3));
modProbes = probes[inModule];
modProbes2 = probes[inModule2];
modProbes3 = probes[inModule3];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
modTOM2 = TOM[inModule2, inModule2];
modTOM3 = TOM[inModule3, inModule3];
dimnames(modTOM) = list(modProbes, modProbes)
dimnames(modTOM2) = list(modProbes2, modProbes2)
dimnames(modTOM3) = list(modProbes3, modProbes3)

# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,
                               nodeAttr = moduleColors[inModule])
vis = exportNetworkToVisANT(modTOM,
                            file = paste("VisANTInput.txt", sep=""),
                            weighted = TRUE,
                            threshold = 0 )
cyt = exportNetworkToCytoscape(modTOM2,
                               edgeFile = paste("CytoscapeInput-edges-all", ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-all", ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes2,
                               nodeAttr = moduleColors[inModule2])
cyt = exportNetworkToCytoscape(modTOM3,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules3, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules3, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes3,
                               nodeAttr = moduleColors[inModule3])

chooseTopHubInEachModule(dataExpr, colorh = genes_colors$mergedColors,power = 8 )

colorh1 <- genes_colors$mergedColors
ADJ1=abs(cor(dataExpr,use="p"))^8
Alldegrees1=intramodularConnectivity(ADJ1, colorh1)
head(Alldegrees1)
y = datTraits1$s8cbl24bd
GS1=as.numeric(cor(y,dataExpr, use="p"))
GeneSignificance=abs(GS1)

colorlevels=unique(colorh1)
sizeGrWindow(9,6)
par(mfrow=c(2,as.integer(0.5+length(colorlevels)/2)))
par(mar = c(4,5,3,1))
for (i in c(1:length(colorlevels))) 
{
  whichmodule=colorlevels[[i]]; 
  restrict1 = (colorh1==whichmodule);
  verboseScatterplot(Alldegrees1$kWithin[restrict1], 
                     GeneSignificance[restrict1], col=colorh1[restrict1],
                     main=whichmodule, 
                     xlab = "Connectivity", ylab = "Gene Significance", abline = TRUE)
}
imconnect <- merge(Alldegrees1, genes_colors, by.x = 0, by.y = 0)
imconnectblue <- imconnect %>% filter(., mergedColors == "blue")
imconnectblue1 <- imconnectblue[,c(1:5,8)]

write.csv(imconnectblue1, "IMConnect_MM_blue.csv", row.names = F)

