#DiffBind analysis of ER ChIP-Seq
library(DiffBind)
library(rgl)
library(DESeq2)
library(dplyr)
library(tidyverse)
library(rtracklayer)

##Load in sampleSheet
setwd("/Library/Frameworks/R.framework/Versions/4.0/Resources/library/DiffBind/extra")
samples <- read.csv("ER_SampleSheet.csv")
names(samples)
samples

##Create DBA from files and sample sheet 
basedir <- system.file("extra", package="DiffBind")
ER_ChIP <- dba(sampleSheet = samples, dir=basedir)

##Correlation heatmap, using occupancy (peak caller score) data
##Heatmap indicates correlation between the location of peaks between samples
par(mar = rep(2, 4))
plot(ER_ChIP)

#blacklist regions hg38
ER_ChIP_Hg38 <- dba.blacklist(ER_ChIP, blacklist = TRUE, greylist = FALSE)
#Genome detected: Hsapiens.UCSC.hg38
#Applying blacklist...
#Removed: 43 of 222560 intervals.
#Removed: 15 merged (of 58199) and 9 (of 39857) consensus

#Count reads in peaks
#Create correlation heatmap, using affinity (read count) data
#Heatmap indicates correlation between the location and signal strength of peaks (number of reads/read depth) between samples
#Defult count method is DESeq2
ER_ChIP_Counts <- dba.count(ER_ChIP)
ER_ChIP_Counts
plot(ER_ChIP_Counts)

#DESeq2 RLE normalisation
ER_ChIP_Norm <-  dba.normalize(ER_ChIP_Counts, normalize=DBA_NORM_NATIVE)

#Create Diferential analysis by defining contrast design for treatments
ER_ChIP_Contrast <- dba.contrast(ER_ChIP_Norm, categories = c(DBA_TREATMENT, DBA_REPLICATE))
ER_ChIP_PDX_DBA <- dba.analyze(ER_ChIP_Contrast)
#Applying Blacklist/Greylists...
#Genome detected: Hsapiens.UCSC.hg38
#Applying blacklist...
#Removed: 8 of 39693 intervals.
#Counting control reads for greylist...
#Building greylist: /Library/Frameworks/R.framework/Versions/4.0/Resources/library/DiffBind/extra/reads/IVC018_Gar15-13_Input.asd.bam
#coverage: 12271616 bp (0.40%)
#Master greylist: 1136 ranges, 12271616 bases
#Removed: 1842 of 39685 intervals.
#Re-normalizing...
#Removed 1850 (of 39693) consensus peaks.
#Analyzing...
#gene-wise dispersion estimates
#mean-dispersion relationship
#final dispersion estimates


op <- par(oma=c(5,7,1,1))
dev.off()
ER_ChIP_PDX_DBA
par(mar = rep(2, 4))
#contrast 1 is Vehicle vs Decitabine
plot(ER_ChIP_PDX_DBA, contrast=1)


#Reterieve differentially bound sites
ER_VEHvsDEC_DB <- dba.report(ER_ChIP_PDX_DBA, contrast = 1)
ER_VEHvsDEC_DB

#save as bed files
write.table(ER_VEHvsDEC_DB, file="ER_VEHvsDEC_DB.bed", quote=F, sep="\t", row.names=F, col.names=F)
#Report Lost peaks
DEC_Lost_ER_Peaks <-ER_VEHvsDEC_DB%>% filter(V9 > 0)
head(DEC_Lost_ER_Peaks)
write.table(DEC_Lost_ER_Peaks, file="ER_Dec_Lost_DBA.bed", quote=F, sep="\t", row.names=F, col.names=F)

#Report gained peaks
DEC_Gained_ER_Peaks <- ER_VEHvsDEC_DB%>% filter(V9 < 0)
head(DEC_Gained_ER_Peaks)
write.table(DEC_Gained_ER_Peaks, file="ER_Dec_Gained_DBA.bed", quote=F, sep="\t", row.names=F, col.names=F)

#VennDiagram
dba.plotVenn(ER_ChIP_PDX_DBA, contrast=1, bDB=TRUE, bGain=TRUE, bLoss=TRUE, bAll=FALSE)

#PCA plots
dba.plotPCA(ER_ChIP_PDX_DBA,DBA_TREATMENT,label=DBA_TREATMENT)
dba.plotPCA(ER_ChIP_PDX_DBA, contrast=1, label=DBA_ID)


#Overlaps
olap.rate <- dba.overlap(ER_ChIP_Contrast,mode=DBA_OLAP_RATE)
plot(olap.rate,type='b',ylab='# peaks',xlab ='Overlap at least this many peaksets')


write.table(consensus_peaks, file="consensus_peaks.bed", quote=F, sep="\t", row.names=F, col.names=F)
