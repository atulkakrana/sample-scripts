### FOr developemntal time point analysis Illumina - 22-Aug-2014

library(lumi)
library(limma)
library(GO.db)
library(affycoretools)
library(annotate)
library(R2HTML)
library(genefilter)

#Illumina------------------------------------------------------------------------------------
getwd()
setwd("/data2/homes/kakrana/3.ProjectR/8.Collab/3.iSyTE2/data_IL/Mut")

filelistA <- c('background_P4_Tdrd7_Lens.txt','background_P30_Tdrd7_Lens.txt',
               'unnormalized_E9.5_P300.txt','unnormalized_E10.5_Pax6.txt','unnormalized_E9.5_Pax6.txt','unnormalized_E10.5_P300.txt') ## Background corrected files 'GSE49223_non-normalized.txt' was used in dev for E 9.5
IL.mydataA <- lumiR.batch(filelistA)

filelistB <- c('unnormalized_iSyTE2.txt','unnormalized_P60_Mafg_2.txt')## Needs BG correction 
IL.mydataB <- lumiR.batch(filelistB)

# filelistC = c('unnormalized_E9.5_P300.txt','unnormalized_E10.5_Pax6.txt','unnormalized_E9.5_Pax6.txt','unnormalized_E10.5_P300.txt')
# IL.mydataC <- lumiR.batch(filelistC)

setwd("../")

## NORMALIZATION #############################################################################
##Illumina------------------------------------------------------------------------------------
sampleNames(IL.mydataA)
IL.eset.Rank.A <- lumiExpresso(IL.mydataA[,-c(5,8,11,15,29)], normalize.param=list(method='rsn'))## Excluded P30Dev_B, P4Tdrd7_B, E95Dev_C (from p300) and E95Pax6_B  - PCA and PCA-3d performed
# par(cex = 0.75)
# boxplot(exprs(IL.eset.Rank.A),log='y',col = 'green')
# plotPCA(IL.eset.Rank.A)
# plotPCA(IL.eset.Rank.A,pcs=c(1,2,3),plot3d=T)
# sampleNames(IL.eset.Rank.A)

sampleNames(IL.mydataB)
IL.eset.Rank.B <- lumiExpresso(IL.mydataB[,-c(13,14)], normalize.param=list(method='rsn')) ## Not using MAFG WB
# par(cex = 0.75)
# boxplot(exprs(IL.eset.Rank.B),log='y',col = 'green')
# plotPCA(IL.eset.Rank.B)
# plotPCA(IL.eset.Rank.B,pcs=c(1,2,3),plot3d=T)

# sampleNames(IL.mydataC)
# IL.eset.Rank.C <- lumiExpresso(IL.mydataC[,-c(3,17)], normalize.param=list(method='rsn')) ## Exclude E95Dev_C (from p300) and E95Pax6_B - PCA and PCA-3d performed
# par(cex = 0.75)
# boxplot(exprs(IL.eset.Rank.C),log='y',col = 'green')
# plotPCA(IL.eset.Rank.C)
# plotPCA(IL.eset.Rank.C,pcs=c(1,2,3),plot3d=T)
# sampleNames(IL.eset.Rank.C)

#Combine Illumina datasets
IL.eset.Rank.AB <- combine(IL.eset.Rank.A, IL.eset.Rank.B)
sampleNames(IL.eset.Rank.AB)

##P/A filtering ####################################################################################
##Illumina 45 final samples 
presentCount.A <- detectionCall(IL.mydataA[,c(5,8,11,15,29)]) ## Excluded P30Dev_B, P4Tdrd7_B,E95Dev_C (from p300) and E95Pax6_B
presentCount.B <- detectionCall(IL.mydataB[,-c(13,14)])## Without MAFG WB
# presentCount.C <- detectionCall(IL.mydataC[,-c(3,17)])## Exclude E95Dev_C (from p300) and E95Pax6_B

presentCount.AB <- presentCount.A+presentCount.B ## Checked if counts are added for same IDs or not - OK

IL.esetRank.Present <- IL.eset.Rank.AB[presentCount.AB >= 2,] ##Probe is present in at least n samples
IL.esetRank.Present

##PROBE REDUCTION #################################################################################
#Illumina------------------------------------------------------------------------
ID <- featureNames(IL.esetRank.Present)
sampleNames(IL.eset.Rank.AB)
testStat <- rowMedians(IL.esetRank.Present[,-c(36,37)]) ## WB columns not included
idUniq <- findLargest(ID, testStat, data = 'lumiMouseAll.db')
esetLumi.Reduced.AB <- IL.eset.Rank.AB[idUniq,] # Extracts not only expression but imported annotation also
esetLumi.Reduced.AB

# Batch Effect Correction #########################################################################
library(sva)
sampleNames(esetLumi.Reduced.AB)
IL.pheno <- read.table('IL.PhenoData_lumiAll.csv', sep =',', header = T); IL.pheno

batch = IL.pheno$Batch
edata <- exprs(esetLumi.Reduced.AB)
mod = model.matrix(~as.factor(Condition), data=IL.pheno)
combat_edata = ComBat(dat=edata, batch=batch, mod=mod, numCovs=NULL, par.prior=TRUE, prior.plots = TRUE)##Originally used
exprs(esetLumi.Reduced.AB) <- combat_edata ##Fitting back to expression set


##DESIGN AND FIT####################################################################################
## Group<- factor(pData(mydata)[,2 ], levels = levels(pData(mydata)[,2]))
# group<- factor(pData(mydata)[,2 ])## Get levels of column 2 a.k.a stages

##Illumina-------------------------------------------------------------------------
sampleNames(esetLumi.Reduced.AB)
IL.sampletype = IL.pheno$SampleType;IL.sampletype

IL.group <-factor(IL.sampletype)
IL.design<- model.matrix(~0+IL.group)##  required to avoid intercept
IL.design


colnames(IL.design) <- unique(IL.pheno$ColNames)
IL.design
IL.fit <- lmFit(esetLumi.Reduced.AB,IL.design)


## CONTRASTS & FIT #################################################################################
##Illumina----------------------------------------------------------------------------
IL.contrast.matrix<-makeContrasts(E95Dev-WB,E95Dev2-WB,E10Dev-WB,E10Dev2-WB,P4Dev-WB,LensP8-WB,
                                  LensP12-WB,LensP20-WB,P30Dev-WB,LensP42-WB,LensP52-WB,P60Dev-WB,
                                  P4Tdrd7-WB,P30Tdrd7-WB,E95P300-WB,E10P300-WB,E95Pax6-WB,E10Pax6-WB,
                                  P60Mafg-WB,P4Tdrd7-P4Dev,P30Tdrd7-P30Dev,E95P300-E95Dev-WB,E10P300-E10Dev,
                                  E95Pax6-E95Dev2,E10Pax6-E95Dev2,P60Mafg-P60Dev,levels=IL.design)##Old WB and not new
IL.fit2<- contrasts.fit(IL.fit,IL.contrast.matrix)
IL.fit2 <- eBayes(IL.fit2)
IL.Intensities <- (2^exprs(esetLumi.Reduced.AB))
IL.fit2$genes <- data.frame(IL.fit2$genes,IL.Intensities)# Check the intensities with raw file

### RESULTS ########################################################################################
##Illumina-----------------------------------------------------------------------------
topTable(IL.fit2, coef = 8, adjust = "fdr", number = 10,sort.by='logFC')
IL.results.default <- decideTests(IL.fit2, method = 'separate',adjust.method = 'BH')##Default: method = separate
IL.results.global <- decideTests(IL.fit2, method = 'global',adjust.method = 'fdr')
IL.results.nested <- decideTests(IL.fit2, method = 'nestedF',adjust.method = 'BH')

write.fit(IL.fit2,IL.results.default,file='IL_ALL_RESULTS_SEP_v1.txt', sep = '\t')
write.fit(IL.fit2,IL.results.nested,file='IL_ALL_RESULTS_NESTF_v1.txt', sep = '\t')
write.fit(IL.fit2,IL.results.global,file='IL_ALL_RESULTS_GLOB_v1.txt', sep = '\t')
