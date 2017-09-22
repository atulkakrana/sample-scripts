## AFFY Developmental
## Atul Kakrana - 31-Aug

library(affy)
library(limma)
library(annotate)
library(genefilter)
library(GO.db)
library(affycoretools)
library(mouse430a2.db)


## INPUT DATA ################################################################################
getwd()
setwd("/data2/homes/kakrana/3.ProjectR/8.Collab/3.iSyTE2/data_AF/430_2A")
AF.mydata <- ReadAffy()# Read the cel files to create an AffyBatch object
setwd("/data2/homes/kakrana/3.ProjectR/8.Collab/3.iSyTE2")# Main Working directory
sampleNames(AF.mydata)

## NORMALIZATION #############################################################################
AF.esetRMA<-rma(AF.mydata) ##P0C and E19.5C removed
AF.esetRMA

plotPCA(AF.esetRMA)
boxplot(exprs(AF.esetRMA),log='y',col = 'green')

##P/A filtering ####################################################################################
## Affy - Total Samples = 36n - 2n (P0C and E19.5 C) = 34
calls <- mas5calls(AF.mydata) ## Get PMA calls
calls.exprs <- exprs(calls) ## Seprate out P/A flags for each probe for every sample
presentCount <- rowSums(calls.exprs=='P') + rowSums(calls.exprs=='M') ## Tested PresentCount+absentCount = total samples - OK
AF.esetRMA.Present<- AF.esetRMA[presentCount >= 2]##The probe should be Present in atleast n samples
AF.esetRMA.Present

##PROBE REDUCTION #################################################################################
ID <- featureNames(AF.esetRMA.Present) ## Probe sets
sampleNames(AF.esetRMA.Present)
testStat <- rowMedians(exprs(AF.esetRMA.Present))##Exclude Whole body in column 32,33,34
idUniq <- findLargest(ID, testStat, data = 'mouse430a2')
esetRMA.Reduced <- AF.esetRMA.Present[idUniq,] 
esetRMA.Reduced

## ANNOTATIONS #####################################################################################
ID <- featureNames(esetRMA.Reduced)
Symbol <- getSYMBOL(ID, 'mouse430a2.db')
Name <- as.character(lookUp(ID, "mouse430a2.db", "GENENAME"))
tmp <- data.frame(ID=ID, Symbol=Symbol, Name=Name,stringsAsFactors=F)
tmp[tmp=="NA"] <- NA
fData(esetRMA.Reduced) <- tmp
esetRMA.Reduced

##DESIGN AND FIT####################################################################################
AF.sampletype <- c('E15','E15','E15BRG','E15BRG','E15BRG','E15BRG','E15','E15','P0','P0','P0','P0HSF4','P0HSF4','P0HSF4','P0Pax6','P0Pax6','P0Pax6','P02','P02','P02')
AF.group <-factor(AF.sampletype)

# Normal Test
AF.design<- model.matrix(~0+AF.group) ## Required to avoid intercept
colnames(AF.design) <- c('E15','E15BRG','P0','P0HSF4','P0Pax6','P02') ## Personalize the design with col names
AF.fit <- lmFit(esetRMA.Reduced,AF.design)
AF.contrast.matrix<-makeContrasts(E15BRG-E15,P0HSF4-P0,P0Pax6-P02,levels=AF.design)
AF.fit2<- contrasts.fit(AF.fit,AF.contrast.matrix)
AF.fit2 <- eBayes(AF.fit2)
AF.Intensities <- (2^exprs(esetRMA.Reduced))
AF.fit2$genes <- data.frame(AF.fit2$genes,AF.Intensities)# Check the intensities with raw file

# ## Combined T-Test
# AF.group
# AF.comb.grp <- relevel(AF.group, ref="WB")
# AF.comb.grp
# AF.comb.design <- model.matrix(~AF.comb.grp)
# AF.fit <- lmFit(esetRMA.Reduced, AF.comb.design)
# AF.comb.cont <- c(-10,rep(1,10))/10 ## Number of samples
# AF.fit2 <- contrasts.fit(AF.fit, contrast=AF.comb.cont)
# AF.fit2 <- eBayes(AF.fit2)
# AF.Intensities <- (2^exprs(esetRMA.Reduced))
# AF.fit2$genes <- data.frame(AF.fit2$genes,AF.Intensities)# Check the intensities with raw file

### RESULTS ########################################################################################
topTable(AF.fit2, coef = 1, adjust = "fdr", number = 10,sort.by='logFC') 
AF.results.default <- decideTests(AF.fit2, method = 'separate',adjust.method='fdr') ## Default: method = separate
AF.results.global <- decideTests(AF.fit2, method = 'global',adjust.method = 'fdr')
AF.results.nested <- decideTests(AF.fit2, method = 'nestedF',adjust.method = 'fdr')


write.fit(AF.fit2,AF.results.default,file='AF430.2A_RESULTS_COMB_SEP_v1.tsv',adjust = 'fdr', F.adjust = 'fdr',sep = '\t')
write.fit(AF.fit2,AF.results.global,file='AF430.2A_RESULTS_COMB_GLOB_v1.tsv',adjust = 'fdr', F.adjust = 'fdr', sep = '\t')
write.fit(AF.fit2,AF.results.nested,file='AF430.2A_RESULTS_COMB_NESTF_v1.tsv',adjust = 'fdr', F.adjust = 'fdr', sep = '\t')
getwd()

