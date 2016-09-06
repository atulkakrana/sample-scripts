### Agilent two-color analysis
### Written by atulkakrana@udel.edu

## Prepare ####################################################
library(limma)
library(annotate)
library(affy)
library(genefilter)
getwd()
setwd('/data2/homes/kakrana/3.ProjectR/3.Anther/pooled')
###############################################################

### 1. Inputs #################################################
targets <- readTargets("pooled.txt")
targets
rawData = read.maimages(targets$FileName, source="agilent")

## For separate channel analysis
targets2 = targetsA2C(targets)
write.table(targets2,file='targets2',sep='\t')
targets2
###############################################################


### 2. Pre-Process ############################################
## BG correct
# MA.corr = backgroundCorrect(rawData, method="normexp",offset=16)
# MA.corr = backgroundCorrect(rawData, method="subtract")
MA.corr = backgroundCorrect(rawData, method="normexp", offset = 16)
plotDensities(MA.corr)

## Normalize within arrays
MA.with = normalizeWithinArrays(MA.corr, method="loess")
plotDensities(MA.with)
## Normalize between arrays for mixed design only
MA.bet <- normalizeBetweenArrays(MA.with, method="Aquantile")
plotDensities(MA.bet)
###############################################################

### Annotate ##################################################
MA.genes = MA.bet$genes ## Get the gene names from 'genes' slot of MA.bet
sysName = MA.genes$SystematicName ## Get the 
length(sysName)

annoFile = read.delim("Maize4x44_Annotate.txt")
names(annoFile); dim(annoFile)
anno = annoFile[,c(1,3,5,6,7,8,9,10)] ## Extract the columns of interest from annotation file

annoRows = anno[with (anno, match(anno$GeneName, sysName, nomatch = FALSE)),] ## Match the rows between genenames from expressionSet with annotation file
dim(annoRows); names(annoRows)
write.table(MA.bet$genes, "annotated2.txt",sep = '\t')

MA.bet$genes = cbind(MA.bet$genes,annoRows[,c(1,3:8)]) ## Removed geneNames column from annoRows as its same as systematic names in MA.bet$genes
names(MA.bet$genes)
###############################################################


### Filter ####################################################
## Filter out absent genes

## Filter out control genes
MA.genes = MA.bet$genes
names(MA.genes)
controlList = MA.genes$SystematicName[MA.genes$Include == 'N']
length(controlList)

# names(MA.bet)
# MA.bet[-list]
## Reduce probes to genes
###############################################################


### 3. Design #################################################
uniqSamples <- unique(targets2$Target)
uniqSamples
fact <- factor(targets2$Target, levels=uniqSamples)
fact
design <- model.matrix(~0+fact)
design
colnames(design) <- uniqSamples
design
###############################################################


### 4. Something ##############################################
corfit <- intraspotCorrelation(MA.bet, design)
fit <- lmscFit(MA.bet, design, correlation=corfit$consensus)
###############################################################


## 6. SAM ####################################################
library(siggenes)
exprs.mat = exprs.MA(MA.bet)
dim(exprs.mat)
exprs.mat[1:10,1:10]
probe.genes = MA.bet$genes
dim(probe.genes)
length(unique(probe.genes[,5]))

Cy3 = targets$Cy3
Cy5 = targets$Cy5
col.names = c(rbind(Cy3, Cy5))
col.names
colnames(exprs.mat) = col.names
colnames(exprs.mat)
exprs.mat[1:10,1:10]

sam.obj = exprs.mat[,c(113,114,115,116,117,118,119,120,129,130,131,132,133,134,135,136,145,146,147,148,149,150,151,152,57,73,104,106,1,4,5,8,10,11,14,15,17,18,19,20,21,22,23,24, 26,27,29,35,38,39,41,42,43,44,49,50,57,72,82,121,122,123,124,125,126,127,128,137,138,139,140,141,142,143,144,153,154,155,156,157,158,159,160)] ## Select columns for test
dim(sam.obj)
sam.cl = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
sam.out = sam(sam.obj,sam.cl, method = d.stat,rand = 123,gene.names = probe.genes[,5])
summary(sam.out)

sam.diff = summary(sam.out,0.9)
sam.diff
plot(sam.out, 0.9)

## Output results
res = sam.diff@mat.sig
names(res);dim(res)
write.table(res,"SamOutput.csv",sep = ',')


# annoRes = anno[with (anno, match(anno$GeneName, rownames(res), nomatch = FALSE)),] ## Match the rows between genenames from expressionSet with annotation file
# dim(annoRes); names(annoRes)

## Output list of differentially expressed genes
diffGenes = list.siggenes(sam.out, 1.1)
diffGenes
################################################################
# 
# 
# ### 7. LIMMA to SAM ###########################################
# 
# sam.limma = limma2sam(fit2, coef  = 1, moderate = TRUE, sam.control = samControl())
# summary(sam.limma)

###############################################################


### 8. Contrasts ##############################################
cont.matrix <- makeContrasts("mm0.4mac1-mm0.4fs",levels=design)
cont.matrix
fit2 <- contrasts.fit(fit, cont.matrix) 
fit2 <- eBayes(fit2)


### Add intensities to Fit
# Intensities <- (exprs.MA(MA.bet))
# Cy3 = targets$Cy3
# Cy5 = targets$Cy5
# col.names = c(rbind(Cy3, Cy5))
# colnames(Intensities) = col.names
# fit2$genes <- data.frame(MA.bet$genes,Intensities)
###############################################################


### 5. Write Results ##########################################
res = topTable(fit2, adjust="BH",sort.by='logFC',n = 50)
write.table(res,file = 'topTable.csv',sep = ',')
res

length(res$logFC>=2)
summary(dt <- decideTests(fit2))
results.default <- decideTests(fit2, method = 'separate',adjust.method = 'BH') ## Default: method = separate
write.fit(fit2,results.default,file='RESULTS_SEP.tsv', sep = ',')

summary(dt <- decideTests(fit2))

###############################################################
