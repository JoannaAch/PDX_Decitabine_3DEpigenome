library(RUVSeq)
library(EDASeq)

## Read in raw counts and phenotype data.
setwd("/Volumes/Joanna_HD/RNA-seq_Level_2/GAR15-13D/")
save.image("/Volumes/Joanna_HD/RNA-seq_Level_2/GAR15-13D/GAR15-13D-RUV_25032020.RData")

RawCounts <- read.csv("tables/GAR15-13D_count_table.csv", row.names = 1)

head(RawCounts)
tail(RawCounts)

#Pheno <- read.csv("Copy_of_Patient_Information_April_2019_RNA-Seq.csv", row.names = 1)

## Filter out genes that are not expressed. Minimum 5 reads in at least 2 samples to remain
filter <- apply(RawCounts, 1, function(x) length(x[x>5])>=2)
filtered <- RawCounts[filter,]
filtered <- round(filtered,0)

genes <- rownames(filtered)[grep("^ENS", rownames(filtered))]

## Add phenotype data to delineate during the differential expression process
x <- as.factor(rep(c("Comb", "Dec", "Tam", "Veh"), each=4))

set <- newSeqExpressionSet(as.matrix(filtered),phenoData = data.frame(x, row.names=colnames(filtered)))

set
##SeqExpressionSet (storageMode: lockedEnvironment)
#assayData: 18383 features, 16 samples 
#element names: counts, normalizedCounts, offset 
#protocolData: none
#phenoData
#sampleNames: GAR15.13D_Comb_4B_RNA GAR15.13D_Comb_5B_RNA ... GAR15.13D_Vehicle_4_RNA (16 total)
#varLabels: x
#varMetadata: labelDescription
#featureData: none
#experimentData: use 'experimentData(object)'
#Annotation:  

#x <- as.factor(Pheno$Group)
#set <- newSeqExpressionSet(as.matrix(filtered),
#                           phenoData = data.frame(x, row.names=colnames(filtered)))
#set

## Determine if inter-sample correction is necessary by plotting
## log-ratio read count against median read count
library(RColorBrewer)
colors <- brewer.pal(3, "Set2")
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotPCA(set, col=colors[x], cex=1.2)

## Between lane normalisation to correct for the above.
## Uses the EDAseq function to upper-quartile normalise
set <- betweenLaneNormalization(set, which="upper")
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotPCA(set, col=colors[x], cex=1.2)

## Create an original design matrix to undergo differential expression
## This is to create a set of non-differentially expressed genes to use as an empirical control
design <- model.matrix(~x, data=pData(set))
y <- DGEList(counts=counts(set), group=x)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)
top <- topTags(lrt, n=nrow(set))$table
empirical <- rownames(set)[which(!(rownames(set) %in% rownames(top)[1:5000]))]



##RUVg -----> MAY NOT BE CORRECT ASSUMPTION AS ALL GENES WOULD BE CHANGING WITH LOSS OF DNA METHYLATION?? (Jo)
#To estimate the factors of unwanted variation, we need a set ofnegative control genes,i.e.,
#genes that can be assumed not to be influenced by the covariates of interes

##Empirical control genes:
#If no genes are knowna priorinot to be influenced by the covariates of interest, one canobtain a set of “in-silico empirical” negative controls,
#e.g., least significantly DE genesbased on a first-pass DE analysis performed prior to RUVg normalization.
## All but the top 5000 genes are considered empirical controls
set2 <- RUVg(set, empirical, k=1)
pData(set2)
plotRLE(set2, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotPCA(set2, col=colors[x], cex=1.2)

## Extract normalised counts data for alternative use
normCounts <- normCounts(set2)

## Differential analysis is undertaken again with empirical controls used if so desired
design <- model.matrix(~x + W_1, data=pData(set2))
y <- DGEList(counts=counts(set2), group=x)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)
topTags(lrt)



#RUVr
#Finally,  a  third  approach  is  to  consider  the  residuals  (e.g.,  deviance  residuals)  from 
#a first-pass GLM regression of the counts on the covariates of interest.
#This can beachieved with theRUVr method

design <- model.matrix(~x, data=pData(set))
y <- DGEList(counts=counts(set), group=x)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)


fit <- glmFit(y, design)
res <- residuals(fit, type="deviance")

set4 <- RUVr(set, genes, k=1, res)
pData(set4)
plotRLE(set4, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotPCA(set4, col=colors[x], cex=1.2)




##RUVs
#RUVsmethod to estimate the factors of unwanted variation using replicate/negative control samples for which the covariates of interest are constant
differences <- makeGroups(x)
differences
#     [,1] [,2] [,3] [,4]
#[1,]    1    2    3    4
#[2,]    5    6    7    8
#[3,]    9   10   11   12
#[4,]   13   14   15   16

set3 <- RUVs(set, genes, k=1, differences)
pData(set3)
plotRLE(set3, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotPCA(set3, col=colors[x], cex=1.2)



###GOExpress using RUVr - set4
library(GOexpress)

RawCounts <- read.csv("tables/GAR15-13D_count_table.csv", row.names = 1)

head(RawCounts)
tail(RawCounts)

dataDirectory <- "/Volumes/Joanna_HD/RNA-seq_Level_2/GAR15-13D/"

pDataFile <- file.path(dataDirectory, "pData.txt")
pData <- read.table(pDataFile,row.names=1, header=TRUE, sep="\t")

summary(pData)

names(pData)

metadata <- data.frame(labelDescription=c("Treatment"), row.names=c("Treatment"))


phenoData <- new("AnnotatedDataFrame", data=pData, varMetadata=metadata)
#An object of class 'AnnotatedDataFrame'
#rowNames: GAR15.13D_Comb_4B_RNA GAR15.13D_Comb_5B_RNA ...
#GAR15.13D_Vehicle_4_RNA (16 total)
#varLabels: Treatment
#varMetadata: labelDescription
#row.names(pData)

set4Set <- ExpressionSet(assayData=counts(set4), phenoData=phenoData)
#ExpressionSet (storageMode: lockedEnvironment)
#assayData: 18383 features, 16 samples 
#element names: exprs 
#protocolData: none
#phenoData
#sampleNames: GAR15.13D_Comb_4B_RNA GAR15.13D_Comb_5B_RNA ...
#GAR15.13D_Vehicle_4_RNA (16 total)
#varLabels: Treatment
#varMetadata: labelDescription
#featureData: none
#experimentData: use 'experimentData(object)'
#Annotation: 
  
  
#Go Analyse without provided GO annotations

set4Set_results <- GO_analyse(eSet = set4Set, f = "Treatment")

#0 features from ExpressionSet found in the mapping table.  
  
set4_results.pVal = pValue_GO(result=set4Set_results, N=100)
  
####### fGSEA

library(fgsea)





