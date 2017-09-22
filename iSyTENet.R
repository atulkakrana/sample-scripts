## Co-expression Networks for ISyTE 2.0
## Script by Atul Kakrana

## Prepare
library(WGCNA);library(org.Mm.eg.db)
getwd()

## Prepare for analysis - default settings
options(stringsAsFactors = FALSE)
enableWGCNAThreads(18)

### 1. READ DATA - That will be used for network construction #############
setwd("/new_data/data2/homes/kakrana/3.ProjectR/8.Collab/3.iSyTE2/Network/AffyCoex2")
# readData = read.delim("./raw/AF_RESULTS_GLOB_RMA_v2_coef1.6_1N.txt")
readData = read.delim("AF_enrichpval_0.05_2n_con-exnet_AF_AF2_IL_merge_non_red_log.txt")
names(readData);readData[1:5,1:5];dim(readData)
"Sept8" %in% readData$Genes.Symbol ## Check if gene is present with correct name


## Extract just the expression matrix and remove and other auxillary information
datExpr0 = as.data.frame(t(readData[,-c(1)])) ## Extract expression Matrix
names(datExpr0) = readData$Genes.Symbol ## Adding the gene Symbols as Colnames
rownames(datExpr0) = names(readData)[-c(1)] ## Adding the samplenames as rownames i.e. Stage names
rownames(datExpr0)
datExpr0[1:6,1:10];dim(datExpr0)

### 2. Check for missing values and identify outlier samples
gsg = goodSamplesGenes(datExpr0, verbose = 3)
names(gsg)
gsg$allOK ## If result is not OK than consult WGCNA tutorial 1 

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

length(names(datExpr0)) ## Length after genes removed
"Sept8" %in% names(datExpr0) ## Check if gene is present with correct name

## Get Ensembl ids - Not required for co-ex net
ens = select(org.Mm.eg.db,keys = as.character(names(datExpr0)),
             columns = "ENSEMBL",keytype = "SYMBOL")
dim(ens);length(names(datExpr0))

## Get and format string data
string = read.delim("10090.protein.links.v9.1.txt",sep = ' ')
dim(string);string[1:5,]

source.name = as.character(string$protein1)
tar.name = as.character(string$protein2)
score = as.numeric(string$combined_score)
length(source.name);length(tar.name);length(score)

source.splt  = strsplit(source.name, "10090.")
source.splt[[1]]
names.source2 = as.data.frame(matrix(unlist(source.splt), ncol=2, byrow=TRUE))
colnames(names.source2) = c("bogus","source")
names.source2[1:5,];dim(names.source2)
final.source = as.character(names.source2$source)

tar.splt = strsplit(tar.name, "10090.")
names.tar2 = as.data.frame(matrix(unlist(tar.splt), ncol=2, byrow=TRUE))
colnames(names.tar2) = c("bogus","target")
names.tar2[1:5,];dim(names.tar2)
final.tar = as.character(names.tar2$target)

final.string = as.data.frame(cbind(final.source,final.tar,score))
final.string[1:5,]

## Extract string data for enriched genes
ens.keys = ens$ENSEMBL
ens.keys[1:5]
string.res = final.string[final.string$final.source %in% ens.keys,]
dim(string.res);string.res[1:5,]


## Cluster Samples
require (flashClust)
sampleTree = flashClust(dist(datExpr0), method = "complete")
par(cex = 0.7)
plot(sampleTree, main = "Clustering on Samples - Complete Linkage",sub = "",xlab="",
     cex.axis = 1.5, cex.lab = 1, cex.main = 2)

## Remove outlier either manuualy from file or while selecting expression amtrix or here:
# abline(h=120,col = "red")
##cluster under this line
# clust = cutreeStatic(sampleTree,cutHeight = 120,minSize = 2)
# table(clust) ## Samples to keep
# keepSamples = (clust==1)
# keepSamples
# datExpr = datExpr0[keepSamples,]
# nGenes = ncol(datExpr)
# nSamples = nrow(datExpr)
# rownames(datExpr)

## If no outlier needs to be removed - take all - as in case of processed microarray data
datExpr = datExpr0

### 3. Adding the metaData/sampleInfo/traitData
phenoData = read.csv("Pheno_file4.csv")
phenoData;dim(phenoData)

## Prepare to combine with expression matrix
lensSamples = rownames(datExpr);lensSamples
phenoRows = match(lensSamples,phenoData$ColumnName);phenoRows
datPheno = phenoData[phenoRows,-1]
datPheno[1:5,] ## Just to see the values and structure of sample Info
rownames(datPheno) = phenoData[phenoRows,1]
datPheno
collectGarbage()

##Save data - till this point
getwd()
save(datExpr,datPheno, file = "isyte-01-dataInput.RData")
lnames = load(file = "isyte-01-dataInput.RData")

### 6. Network Construction #############################################################

## Choose a set of soft-thersholding power and call network topology analysis function
powers = c(c(1:10),seq(from = 12, to=20,by=2))
sft = pickSoftThreshold(datExpr,powerVector = powers, verbose =5)
names(sft)

### Scale-free topology fit
plot(sft$fitIndices[,1],-sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab = "Soft Threshold (Power)", ylab = "Scale Free Topology Model Fit, signed R^2", type = "n",main = paste("Scale Independence"))
text(sft$fitIndices[,1],-sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels = powers,col = "red")
## This line will correposnd to the cutoff to using an R^2 cut-off of h
abline(h=0.90,col="red")
##Mean connectivity as a function of soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers,col="red")

##Choose power from scale-free topology fit, lowet value at which plot flattens out - need to dig more
pwr = 7 ## Found pwr =10 on Apr10-15
net = blockwiseModules(datExpr, power = pwr, TOMtype = "unassigned", minModuleSize = 30,
                       deepSplit = 4, reaasignThreshold = 0, mergeCutHeight = 0.25,numericLabels=TRUE
                       ,pamRespectDendro = FALSE,saveTOMs = TRUE,saveTOMFileBase = "isyteTOM",verbose =3)
names(net)
## How many modules identifies and  their sizes
table(net$colors)
mergedColors = labels2colors(net$colors) ## Convert labels to colors for plotting
# sizeGrWindow(12,9)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module Colors",
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)

## Save the progress till this point
moduleLabels = net$colors; moduleColors = labels2colors(net$colors); MEs = net$MEs
geneTree = net$dendrograms[[1]]
save(MEs,moduleLabels,moduleColors,geneTree,file = "isyte-02-networkConstruction-auto.RData")

### 7. Exporting Data ##########################################################################
##Load
lnames = load(file = "isyte-02-networkConstruction-auto.RData")

## How many networks
table(moduleLabels); unique(moduleColors) ## All the colors i.e. networks

## Calculate topology overlap
TOM = TOMsimilarityFromExpr(datExpr,power=pwr)
TOM [1:10,1:10]

## Choose module by colors
modules = unique(moduleColors);modules ## All the colors i.e. networks
# modules = c("brown","blue","turquoise","red","green","black","grey","purple","yellow","pink","magenta")
inModule = is.finite(match(moduleColors,modules)) ## Extract modules of interest
probes = names(datExpr) ## All the probes
modProbes = probes[inModule] ## Probes from selected modules
modProbes[1:10]

## Read annotation file
# annot = read.delim("./raw/AF_RESULTS_GLOB_RMA_v2_coef1.6_1N_Anno.txt");dim(annot)

setwd("/data2/homes/kakrana/3.ProjectR/8.Collab/3.iSyTE2/Clust")
affy.data = read.delim("AF_Testpvals_rand1000_track10_int6MB_wt5_09_07_03_56_Final_FC_v2.txt") ## Gene Symbol as rowname
names(affy.data[,1:4])
"Sept8" %in% affy.data$Genes.Symbol ## Check if gene is present with correct name
annot = affy.data[affy.data$Genes.Symbol %in% names(datExpr),1:3]
annot[1:5,];dim(annot);length(names(datExpr))
not.found = names(datExpr)[!names(datExpr) %in% affy.data$Genes.Symbol]
not.found

modGenes = annot$Genes.Symbol[match(modProbes,annot$Genes.Symbol)] ## Annotate by matching modProbes with Probe ID column in annotation file
modGenes[1:10]

##Select the correponding Topological overlap
modTOM = TOM[inModule,inModule]
dimnames(modTOM) = list(modProbes,modProbes)

## Final Export
getwd()
thres = 0.20
cyt = exportNetworkToCytoscape(modTOM,edgeFile = paste("CytoscapeInputv2-edges-",paste(modules,collapse="-"),"-thres-",thres,".txt",sep=""),nodeFile = paste("CytoscapeInputv2-nodes-",paste(modules,collapse="-"),"-thres",thres,".txt",sep=""),weighted = TRUE,threshold = thres,nodeNames = modProbes,altNodeNames = modGenes,nodeAttr = moduleColors[inModule])

## Selecting subnetwork of top hub genes
nTop = 50
IMConn = softConnectivity(datExpr[,modProbes])
top = (rank(-IMConn)<= nTop)
cyt = exportNetworkToCytoscape(modTOM[top,top],edgeFile = paste("CytoscapeInputv2-edges-",paste(modules,collapse="-"),"-thres",thres,"-top.txt",sep=""),nodeFile = paste("CytoscapeInputv2-nodes-",paste(modules,collapse="-"),".-thres",thres,"-top.txt",sep=""),weighted = TRUE,threshold = thres,nodeNames = modProbes,altNodeNames = modGenes,nodeAttr = moduleColors[inModule])


## Visualize gene network - heatmap
dissTOM = 1-TOM
plotTOM = dissTOM^pwr
diag(plotTOM) = NA
TOMplot(plotTOM,geneTree,moduleColors,main = "Network heatmap plot,all genes")
