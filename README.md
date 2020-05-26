# JBMR_2019
Code for JBMR 2020 Paper

Analysis Steps:
1. Format Phenotype Data: Compiles and Filters Phenotype data
2. Statistical Analysis: miRNA Multiple Linear Regression Analysis
3. Unsupervised Analysis: miRNA clustering and enrichment of clusters
4. mRNA Correlation: Determine mRNA that are targets of and inversely correlated with miRNA of interest
5. Enrichment analysis: Enrichment of miR-19a and miR-186
6. Scripts for Figures: Scripts for data for specific tables and figures


## Format Phenotype Data
### PhenoData_ISPE_KN.R 
* creates phenotype file used in data analysis
* integrates medication data, bone mineral density, and clinical covariates data from Exam 8
* filters for BMD after  Clinical Exam 8

## Statistical Analysis
### MLRTechRank.R
* Statistical analysis of miRNA association with BMD and BB use
* Using linear and logistic regression, adjusting for age, sex, height, weight, and miRNA technical variables 
* Model 1: Multiple linear regression, miRNAs < 10% NA (expressed in 90% of samples)
* Model 2: Multiple logistic regression, miRNAs with 90-95% NA (expressed in 5-10% of samples)
* Model 3: linear and logistic regression, miRNAs with 10-90% NAs (expressed in 10 to 90% of samples)
* Combining linear and logistic p-value using Fisher's method for Model 3 miRNAs

### JBMR_eBMD_validation.R
* analyze eBMD data from external cohort
* calculation correlation between miRNAs and Z-score for meta-analysis

### Z-score_function.R
* function to calculate Z-score in FHS data
* group people by age
* calculate mean and standard deviation (sd) BMD measure for each age group
* Z-score = (BMD for each individual - mean for group)/sd for group

### Z-score calc.R
* calculate Z-score in males and females for FN BMD
* calculate correlation between Z-score and miRNAs

### meta_analysis.R
* perform meta-analysis for Z-score correlation with miRNAs for external cohort and FHS
* create forest plots

### RevisionAnalyses.R
* Additional statistical analyses for revisions
* T-score differences between BB users and non-users
* New Boxplots for figure 2
* Correlation between miRNA and BMD/T-score
* Sensitivity analysis for hypertension
* miRNAs in hypertension treated patients
* miRNAs association with B1-specific BBs
* miRNAs association with B1 vs non-specific BBs vs non-BB users
* Perform power calculation
* Sensitivity analysis for BP
* Sensitivity analysis for BP and hypertension
* Correlation between BP and BB use/Hypertension-treatment
* Explained variance
* % difference in BMD between BB users and non-users

## Unsupervised Analysis
### clustering_nofilter.R
* Perform WGCNA clustering on miRNA data
* Determined association of clusters with BB use and BMD
* Create network files for Cytoscape

### AllClusterEnrichment.R
* Perform enrichment analysis for each cluster

## mRNA Correlation
### mRNAcormiRNA
* determine validated target mRNAs of top miRNAs that are negatively correlated

## Enrichment Analysis
### GSEACor
* Perform enrichment of mRNA targets of top miRNAs in GO and KEGG

## Scripts for Figures
### Figure 2A and B.R
* Create figure 2A and 2B: association between miRNAs and BMD in BB users and non-users

### Figure 5
* Create figure 5: Enrichment of top bone-related terms in miR-19a, miR-186, and blue cluster

### Figure S2A
* Plot miRNA expression in tissues using data from Tissue Atlas

### Table 1
* Create Table 1 with the exception of T-score which was added in RevisionsAnalyses.R

### Table 2 
* Create Table 2 with the exception of T-score which was added in RevisionsAnalyses.R

### Table 3
* Output data for Table 3, similar to MLRTechRank.R


