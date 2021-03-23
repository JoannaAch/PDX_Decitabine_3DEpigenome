library(edgeR)
library(ggplot2)
library(gridExtra)
library(data.table)
library(goseq)
library(GO.db)
library(dplyr)


setwd("/Volumes/Hippo/NCI_Epigenetics/RNA-Seq_Level_2/hg38/PDX_Decitabine/GAR15-13D/")
load("GAR15-13D_edgeR_24032020.RData")
save.image("GAR15-13D_edgeR_24032020.RData")


SAMPLES.NAME <- "GAR15-13D"

OUTPUT.FOLDER <- "DGE"

TABLES.FOLDER <- "tables/"

raw.read.counts <- read.csv(paste0(TABLES.FOLDER, SAMPLES.NAME, "_count_table.csv"), row.names=1)
tpm.data <- read.csv(paste0(TABLES.FOLDER, SAMPLES.NAME, "_TPM_table.csv"), row.names=1)

# round the rsem gene expected counts values to the nearest integer to input into edgeR
raw.read.counts <- round(raw.read.counts)
#counts<-counts[, c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17)]
# matrix of counts with ENSGxxxxxxxx tags


counts <- data.matrix(raw.read.counts) 
colnames(counts)<-c("GAR15.13D_Comb_rep1", "GAR15.13D_Comb_rep2", "GAR15.13D_Comb_rep3", "GAR15.13D_Comb_rep4", "GAR15.13D_Dec_rep1", "GAR15.13D_Dec_rep2", "GAR15.13D_Dec_rep3", "GAR15.13D_Dec_rep4", "GAR15.13D_Tam_rep1", "GAR15.13D_Tam_rep2", "GAR15.13D_Tam_rep3", "GAR15.13D_Tam_rep4", "GAR15.13D_Vehicle_rep1", "GAR15.13D_Vehicle_rep2", "GAR15.13D_Vehicle_rep3", "GAR15.13D_Vehicle_rep4")
counts<-counts[, c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)]

colnames(tpm.data)<-c("GAR15.13D_Comb_rep1", "GAR15.13D_Comb_rep2", "GAR15.13D_Comb_rep3", "GAR15.13D_Comb_rep4", "GAR15.13D_Dec_rep1", "GAR15.13D_Dec_rep2", "GAR15.13D_Dec_rep3", "GAR15.13D_Dec_rep4", "GAR15.13D_Tam_rep1", "GAR15.13D_Tam_rep2", "GAR15.13D_Tam_rep3", "GAR15.13D_Tam_rep4", "GAR15.13D_Vehicle_rep1", "GAR15.13D_Vehicle_rep2", "GAR15.13D_Vehicle_rep3", "GAR15.13D_Vehicle_rep4")
tpm.data<-tpm.data[, c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)]
tpm.data.matrix <- data.matrix(tpm.data)

#sample.names <- colnames(raw.read.counts)
sample.names <- colnames(c("GAR15.13D_Comb_rep1", "GAR15.13D_Comb_rep2", "GAR15.13D_Comb_rep3", "GAR15.13D_Comb_rep4", "GAR15.13D_Dec_rep1", "GAR15.13D_Dec_rep2", "GAR15.13D_Dec_rep3", "GAR15.13D_Dec_rep4", "GAR15.13D_Tam_rep1", "GAR15.13D_Tam_rep2", "GAR15.13D_Tam_rep3", "GAR15.13D_Tam_rep4", "GAR15.13D_Vehicle_rep1", "GAR15.13D_Vehicle_rep2", "GAR15.13D_Vehicle_rep3", "GAR15.13D_Vehicle_rep4"))
# range of library size/sequencing depth
library.size <- round(colSums(counts)/1e6, 1)
library.size
# GAR15.13D_Comb_rep4b   GAR15.13D_Comb_rep5b    GAR15.13D_Comb_rep7    GAR15.13D_Comb_rep8    GAR15.13D_Dec_rep10     GAR15.13D_Dec_rep4     GAR15.13D_Dec_rep5     GAR15.13D_Dec_rep9     GAR15.13D_Tam_rep3     GAR15.13D_Tam_rep4     GAR15.13D_Tam_rep5 
#27.0                   27.3                   25.8                   23.1                   27.8                   32.5                   25.2                   27.4                   27.7                   26.7                   24.7 
#GAR15.13D_Tam_rep6 GAR15.13D_Vehicle_rep1 GAR15.13D_Vehicle_rep2 GAR15.13D_Vehicle_rep3 GAR15.13D_Vehicle_rep4 
#24.1                   28.2                   27.2                   27.1                   22.1 



# perform PCA on log transformed count values
pca <- prcomp(t(log(counts+1)), center=TRUE)
# also see summary(pca)
pc1.var <- round(pca$sdev[1]^2/sum(pca$sdev^2)*100,2) 
pc2.var <- round(pca$sdev[2]^2/sum(pca$sdev^2)*100,2)
pc.data.frame <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], Name=colnames(counts), stringsAsFactors=F)

makeLab <- function(x,pc) paste0("PC",pc,": ", x, "% variance")

ggplot(pc.data.frame, aes(x=PC1, y=PC2, label=Name))+
  geom_point() + geom_text(aes(colour = factor(gsub("rep1|rep2|rep3|rep4", "", pc.data.frame$Name)))) +
  xlab(makeLab(pc1.var,1)) + ylab(makeLab(pc2.var,2)) +
  ggtitle("PC1 vs PC2 for log(count) of GAR15-13D") +
  expand_limits(x=c(-100, 100)) +
  theme(legend.position = "none")

# perform PCA on TPM values (use log(TPM) as TPM can be quite swayed by 
# differences in sample library size)
pca <- prcomp(t(log(tpm.data.matrix+1)), center=TRUE)
# also see summary(pca)
pc1.var <- round(pca$sdev[1]^2/sum(pca$sdev^2)*100,2) 
pc2.var <- round(pca$sdev[2]^2/sum(pca$sdev^2)*100,2)
pc.data.frame <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], Name=colnames(tpm.data.matrix))

ggplot(pc.data.frame, aes(x=PC1, y=PC2, label=Name))+
  geom_point() + geom_text(aes(colour=factor(gsub("rep1|rep2|rep3|rep4", "", pc.data.frame$Name)))) +
  xlab(makeLab(pc1.var,1)) + ylab(makeLab(pc2.var,2)) +
  ggtitle("PC1 vs PC2 for TPM values of GAR15-13D") +
  expand_limits(x=c(-50,50)) +
  theme(legend.position = "none")




# hierarchical clustering for log(TPM) vlaues
hc <- hclust(dist(t(log(tpm.data.matrix+1))))       
plot(hc, labels=colnames(tpm.data.matrix), main="Clustering samples on log(TPM+1)")

hc <- hclust(dist(t(log(counts+1))))       
plot(hc, labels=colnames(counts), main="Clustering samples on log(counts+1)")

# Filter out ENSGxxxx tags whose coverage is so low that any group differences 
# aren't truly "real". 
# filter out tags whose rowcount <= degrees of freedom.
counts <- counts[rowSums(counts) >= 3,]
tpm.data.matrix <- tpm.data.matrix[rowSums(tpm.data.matrix) >= 3, ]

# set up design matrix
group.types <- gsub("_rep1|_rep2|_rep3|_rep4", "", colnames(counts))

group <- factor(group.types)

design <- model.matrix(~0+group) 
colnames(design) <- gsub("group", "", colnames(design))

#   GAR15.13D_Comb GAR15.13D_Dec GAR15.13D_Tam GAR15.13D_Vehicle
#1               1             0             0                 0
#2               1             0             0                 0
#3               1             0             0                 0
#4               1             0             0                 0
#5               0             1             0                 0
#6               0             1             0                 0
#7               0             1             0                 0
#8               0             1             0                 0
#9               0             0             1                 0
#10              0             0             1                 0
#11              0             0             1                 0
#12              0             0             1                 0
#13              0             0             0                 1
#14              0             0             0                 1
#15              0             0             0                 1
#16              0             0             0                 1
#attr(,"assign")
#[1] 1 1 1 1
#attr(,"contrasts")
#attr(,"contrasts")$group
#[1] "contr.treatment"



dge_gene <- DGEList(counts=counts, group=group)
dge_gene$samples


###Filter out lowly expressed genes
cpm.cutoff <- 10/(min(dge_gene$samples$lib.size)/1000000)
cpm.cutoff
#0.4531546'
##Use 0.6 for simplicity
keep <- rowSums(cpm(dge_gene) > 0.6) >= 4 # 4 is used here as each group as 4 replicates
table(keep)

#keep
#FALSE  TRUE 
#10811 15887 


# TMM normalization
dge_gene_norm1 <- calcNormFactors(dge_gene, method="TMM")
dge_gene_norm1$samples


# perform PCA on log(CPMS)
cpms <- cpm(dge_gene_norm1)
log.cpms <- log(cpms + 1)
# perform PCA on CPM values (use log(CPM) as CPM can be quite affected by differences in 
# sample library size)
pca <- prcomp(t(log.cpms), center=TRUE)
# also see summary(pca)
pc1.var <- round(pca$sdev[1]^2/sum(pca$sdev^2)*100,2) 
pc2.var <- round(pca$sdev[2]^2/sum(pca$sdev^2)*100,2)
pc.data.frame <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], Name=colnames(cpms))






# The calcNormFactors() function normalizes for RNA composition by finding 
# a set of scaling factors for the library sizes that minimize the log-fold 
# changes between the samples for most genes. 
# The default method for computing these scale factors uses a trimmed mean 
# of M-values (TMM) between each pair of samples. It is based on the hypothesis 
# that most genes are not DE.(However, by switching off acetylation, you might
# expect to get a lot of down-regulated genes so TMM may not be the most
# appropriate norm method)



# perform PCA on log(CPMS)
cpms <- cpm(y)
log.cpms <- log(cpms + 1)

# perform PCA on CPM values (use log(CPM) as CPM can be quite affected by differences in 
# sample library size)

pca <- prcomp(t(log.cpms), center=TRUE)
# also see summary(pca)
pc1.var <- round(pca$sdev[1]^2/sum(pca$sdev^2)*100,2) 
pc2.var <- round(pca$sdev[2]^2/sum(pca$sdev^2)*100,2)
pc.data.frame <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], Name=colnames(cpms))


# GLM estimates of variance (dispersion)

# Dispersion
dge_gene_norm1 <- estimateDisp(dge_gene_norm1, design)
dge_gene_norm1$common.dispersion




# Fitting a model in edgeR takes several steps. 
# First, you must fit the common dispersion. Then you need to fit a trended model 
# (if you do not fit a trend, the default is to use the common dispersion as a trend). 
# Then you can fit the tagwise dispersion which is a function of this model.
dge_gene_norm1 <- estimateGLMCommonDisp(dge_gene_norm1, design, verbose=TRUE)
#Disp = 0.06275 , BCV = 0.2505 
dge_gene_norm1 <- estimateGLMTrendedDisp(dge_gene_norm1, design)
dge_gene_norm1 <- estimateGLMTagwiseDisp(dge_gene_norm1, design)

# plot the genewise biological coefficient of variation (BCV) against gene abundance 
# (in log2 counts per million). 
# It displays the common, trended and tagwise BCV estimates.
plotBCV(dge_gene_norm1)

# Fitting linear model
fit <- glmFit(dge_gene_norm1, design)

# find the tags that are interesting by using a LRT (Likelihood Ratio Test)
# alternative to use ? makeContrasts (limma)

# absolute differences between various pairs of conditions

##Vehicle is -1
lrt.Veh.vs.Dec <- glmLRT(fit, contrast=c(0, 1, 0, -1))

lrt.Veh.vs.Tam  <- glmLRT(fit, contrast=c(0, 0, 1, -1))

lrt.Veh.vs.Comb  <- glmLRT(fit, contrast=c(1, 0, 0, -1))


##GLM
fit_glm <- glmQLFit(dge_gene_norm1, design)

glm.Veh.vs.Dec <- glmQLFTest(fit_glm, contrast=c(0, 1, 0, -1))
glm.Veh.vs.Tam <- glmQLFTest(fit_glm, contrast=c(0, 0, 1, -1))
glm.Veh.vs.Comb <- glmQLFTest(fit_glm, contrast=c(1, 0, 0, -1))



# add all contrasts to a list for subsequent processing
lrt.list <- list(lrt.Veh.vs.Dec=lrt.Veh.vs.Dec, 
                 lrt.Veh.vs.Tam=lrt.Veh.vs.Tam,
                 lrt.Veh.vs.Comb=lrt.Veh.vs.Comb)

glm.list <- list(glm.Veh.vs.Dec=glm.Veh.vs.Dec, 
                 glm.Veh.vs.Tam=glm.Veh.vs.Tam,
                 glm.Veh.vs.Comb=glm.Veh.vs.Comb)


# annotation
# FROM GTF FILE
gtf.annotation.file <- "data/gene_annotation.tsv"
gtf.anno <- read.table(gtf.annotation.file, sep= "\t", header=TRUE, stringsAsFactors=F)

# annotation file downloaded from http://www.ensembl.org/biomart/martview/ (July2017)
# has following columns: ensembl.gene.ID, chr, start, end, strand, 
# description, HGNC.symbol, entrez.gene.ID
annotation.file <- "data/gene_annotation_mart_export.tsv"
if (!(file.exists(annotation.file))){
  untar(paste0("input/", annotation.file, ".tgz"))
}

DT <- fread(annotation.file)
setnames(DT, gsub(" ", ".", colnames(DT)))
setkey(DT, HGNC.symbol)

# implement once here for use with repeated GOseq analysis within lapply
# GOseq
lengths.gene <- read.csv(paste0(TABLES.FOLDER, SAMPLES.NAME, "_effective_length_table.csv"), row.names=1)
# get the average gene length for each row
lengths <- apply(lengths.gene, 1, mean)

getgoresults <- function(DMgenes, bias.data, output.folder){
  
  # fitting the Probability Weighting Function
  # PWF quantifies how the probability of a gene selected as DE changes as 
  # a function of its transcript length
  pwf <- nullp(DEgenes=DMgenes, genome="hg38", id="geneSymbol", 
               bias.data=bias.data, plot.fit=FALSE)
  pdf(paste0(output.folder, "pwf.goodness.of.fit.plot.pdf"))
  plotPWF(pwf)
  dev.off()
  # calculate  the  over  and  under  expressed  GO
  # categories among DE genes
  # goseqres ordered by GO category over representation amongst DE genes.
  goseqres <- goseq(pwf, "hg38", "geneSymbol")
  # multiple correction
  goseqres$over_fdr <- p.adjust(goseqres$over_represented_pvalue, method="BH")
  goseqres$under_fdr <- p.adjust(goseqres$under_represented_pvalue, method="BH")
  
  over <- goseqres[order(goseqres$over_fdr),]
  
  go <- getgo(names(DMgenes),"hg38","geneSymbol")
  go <- go[DMgenes==1]
  go <- go[!is.na(names(go))]
  # create array of gene name-GO term
  gotable <- array(0, c(0, 2))
  for (i in 1:length(go)){
    gotable <- rbind(gotable, cbind(rep(names(go)[[i]], length(go[[i]])), unlist(go[[i]])))
  }
  # split by go term into genes
  gomap <- split(gotable[,1], gotable[,2])
  # add to goseq results table
  m <- match(goseqres$category, names(gomap))
  over$genes <- sapply(m, function(x) paste(unlist(gomap[x]), collapse=', '))
  over
}

FDR.CUTOFF <- 0.1
LOG.FC.CUTOFF <- 1.5

res <- lapply(names(lrt.list), function(contrast){
  
  print (contrast)
  comparison <- gsub("lrt.", "", contrast)
  
  partA <- strsplit(comparison, ".vs.")[[1]][1]
  partB <- strsplit(comparison, ".vs.")[[1]][2]
  
  SUB.FOLDER <- paste0(OUTPUT.FOLDER, gsub("lrt.", "", contrast), "/")
  
  if (!(file.exists(SUB.FOLDER))){
    dir.create(SUB.FOLDER)
  }
  # Top table
  # the default method used to adjust p-values for multiple testing is BH.
  tt <- topTags(lrt.list[[contrast]], n=nrow(counts))$table
  
  m <- match(rownames(tt), gtf.anno$gene.id)
  
  # assign the rownames(tt) as the gene_id as more specific with version 
  # number at end as originates from original gtf
  tt$gene.id <- rownames(tt)
  tt$gene.symbol <- gtf.anno$gene.name[m]
  tt$chr <- gtf.anno$chr[m]
  tt$start <- gtf.anno$start[m]
  tt$end <- gtf.anno$end[m]
  tt$strand <- gtf.anno$strand[m]
  tt$gene.type <- gtf.anno$gene.type[m]
  
  m <- match(gsub("\\.[0-9]*", "", rownames(tt)), DT$Gene.stable.ID)
  
  tt$description <- DT$Gene.description[m]
  tt$entrez.gene.id <- DT$Gene.name[m]
  
  # only keep chromosome names beginning with chr1..22, X, Y; 
  # remove patch chromosome assignments 
  # like JH806587.1, JH806587.1 etc
  tt <- tt[grep("chr*", tt$chr),]
  
  #Volcano plot
  #  plot(tt$logFC, -log10(tt$PValue), type="n", xlab="BRG1 KD <- -> scrambled logFC", ylab="-log10(p.value)", main="Volcano plot of Scrambled vs BRG1 KD")
  #plot(tt$logFC, -log10(tt$PValue), type="n", xlab=paste0(partB, " <- -> ", partA, " logFC"), ylab="-log10(p.value)", main=paste0("Volcano plot of ", partA, " vs ", partB))
  #text(tt$logFC, -log10(tt$PValue), labels = tt$gene.symbol, cex=0.5)
  #abline(h=-log10(tt$PValue[sum(tt$FDR < 0.05)]), col="red")
  
  vp.data <- tt[c("gene.symbol", "logFC", "PValue", "FDR")]
  vp.data = mutate(vp.data, sig=ifelse(((vp.data$FDR<FDR.CUTOFF)&(abs(vp.data$logFC)>LOG.FC.CUTOFF)),
                                       paste0("FDR<", FDR.CUTOFF), "Not Sig"))
  
  p <- ggplot(vp.data, aes(logFC, -log10(PValue))) +
    geom_point(aes(col=sig)) +
    scale_color_manual(values=c("red", "black")) +
    geom_text(data=filter(vp.data, ((FDR<FDR.CUTOFF)&(abs(logFC)>LOG.FC.CUTOFF))), aes(label=gene.symbol)) +
    ggtitle(paste0(comparison, " volcano plot")) +
    # edit the x label to signify contrast direction
    xlab(paste0(partB, " <- ", "logFC", " -> ", partA))
  
  print (p)
  
  # The function plotSmear generates a plot of the tagwise log-fold-changes against 
  # log-cpm (analogous to an MA-plot for microarray data). 
  # DE tags are highlighted on the plot
  de2 <- decideTestsDGE(lrt.list[[contrast]], p.value=FDR.CUTOFF, lfc=LOG.FC.CUTOFF)
  de2tags <- rownames(y)[as.logical(de2)]
  plotSmear(lrt.list[[contrast]], de.tags=de2tags, 
            main=paste0("smear plot with p.value < ", FDR.CUTOFF, " and LFC=", 
                        LOG.FC.CUTOFF, " cutoffs"))
  abline(h = c(-2, 2), col = "blue")
  
  # defining significant as FDR < FDR.CUTOFF and abs(logFC) > LOG.FC.CUTOFF 
  # add DGE.status column
  # UP, DOWN, NC
  tt$DGE.status <- "NC"
  # defining significant as FDR < FDR.CUTOFF and abs(logFC) > LOG.FC.CUTOFF 
  # for UP/DOWN
  if (nrow(tt[((tt$FDR < FDR.CUTOFF)&(tt$logFC > LOG.FC.CUTOFF)),]) > 0){
    tt[((tt$FDR < FDR.CUTOFF)&(tt$logFC > LOG.FC.CUTOFF)),]$DGE.status <- "UP"
  }
  if (nrow(tt[((tt$FDR < FDR.CUTOFF)&(tt$logFC < -LOG.FC.CUTOFF)),]) > 0){
    tt[((tt$FDR < FDR.CUTOFF)&(tt$logFC < -LOG.FC.CUTOFF)),]$DGE.status <- "DOWN"
  }
  
  # filtering/re-ordering
  column.order <- c(6:13, 1:5)
  tt <- tt[column.order]
  tt[is.na(tt)] <- ""
  
  write.table(tt, paste0(SUB.FOLDER, comparison, ".DGE.tsv"), sep="\t", quote=F, row.names=F)
  
  print (nrow(tt[((tt$FDR < FDR.CUTOFF)&(abs(tt$logFC) > LOG.FC.CUTOFF)),]))
  print (paste0("UP: ", nrow(tt[(tt$DGE.status=="UP"),])))
  print (paste0("DOWN: ", nrow(tt[(tt$DGE.status=="DOWN"),])))
  
  # defining significant as FDR < FDR.CUTOFF and abs(logFC) > LOG.FC.CUTOFF
  sigtt <- tt[((tt$FDR < FDR.CUTOFF)&(abs(tt$logFC) > LOG.FC.CUTOFF)),]
  
})

#[1] "lrt.Veh.vs.Dec"
#[1] 527
#[1] "UP: 334"  ####FC > 0
#[1] "DOWN: 193"
#[1] "lrt.Veh.vs.Tam"
#[1] 6151
#[1] "UP: 4729"
#[1] "DOWN: 1422"
#[1] "lrt.Veh.vs.Comb"
#[1] 5809
#[1] "UP: 4311"
#[1] "DOWN: 1498"

save.image("GAR15-13D_edgeR_LTR_GLM_20102020.RData")

##GLM
res <- lapply(names(glm.list), function(contrast){
  
  print (contrast)
  comparison <- gsub("glm.", "", contrast)
  
  partA <- strsplit(comparison, ".vs.")[[1]][1]
  partB <- strsplit(comparison, ".vs.")[[1]][2]
  
  SUB.FOLDER <- paste0(OUTPUT.FOLDER, gsub("glm.", "", contrast), "/")
  
  if (!(file.exists(SUB.FOLDER))){
    dir.create(SUB.FOLDER)
  }
  # Top table
  # the default method used to adjust p-values for multiple testing is BH.
  tt <- topTags(glm.list[[contrast]], n=nrow(counts))$table
  
  m <- match(rownames(tt), gtf.anno$gene.id)
  
  # assign the rownames(tt) as the gene_id as more specific with version 
  # number at end as originates from original gtf
  tt$gene.id <- rownames(tt)
  tt$gene.symbol <- gtf.anno$gene.name[m]
  tt$chr <- gtf.anno$chr[m]
  tt$start <- gtf.anno$start[m]
  tt$end <- gtf.anno$end[m]
  tt$strand <- gtf.anno$strand[m]
  tt$gene.type <- gtf.anno$gene.type[m]
  
  m <- match(gsub("\\.[0-9]*", "", rownames(tt)), DT$Gene.stable.ID)
  
  tt$description <- DT$Gene.description[m]
  tt$entrez.gene.id <- DT$Gene.name[m]
  
  # only keep chromosome names beginning with chr1..22, X, Y; 
  # remove patch chromosome assignments 
  # like JH806587.1, JH806587.1 etc
  tt <- tt[grep("chr*", tt$chr),]
  
  #Volcano plot
  #  plot(tt$logFC, -log10(tt$PValue), type="n", xlab="BRG1 KD <- -> scrambled logFC", ylab="-log10(p.value)", main="Volcano plot of Scrambled vs BRG1 KD")
  #plot(tt$logFC, -log10(tt$PValue), type="n", xlab=paste0(partB, " <- -> ", partA, " logFC"), ylab="-log10(p.value)", main=paste0("Volcano plot of ", partA, " vs ", partB))
  #text(tt$logFC, -log10(tt$PValue), labels = tt$gene.symbol, cex=0.5)
  #abline(h=-log10(tt$PValue[sum(tt$FDR < 0.05)]), col="red")
  
  vp.data <- tt[c("gene.symbol", "logFC", "PValue", "FDR")]
  vp.data = mutate(vp.data, sig=ifelse(((vp.data$FDR<FDR.CUTOFF)&(abs(vp.data$logFC)>LOG.FC.CUTOFF)),
                                       paste0("FDR<", FDR.CUTOFF), "Not Sig"))
  
  p <- ggplot(vp.data, aes(logFC, -log10(PValue))) +
    geom_point(aes(col=sig)) +
    scale_color_manual(values=c("red", "black")) +
    geom_text(data=filter(vp.data, ((FDR<FDR.CUTOFF)&(abs(logFC)>LOG.FC.CUTOFF))), aes(label=gene.symbol)) +
    ggtitle(paste0(comparison, " volcano plot")) +
    # edit the x label to signify contrast direction
    xlab(paste0(partB, " <- ", "logFC", " -> ", partA))
  
  print (p)
  
  # The function plotSmear generates a plot of the tagwise log-fold-changes against 
  # log-cpm (analogous to an MA-plot for microarray data). 
  # DE tags are highlighted on the plot
  de2 <- decideTestsDGE(glm.list[[contrast]], p.value=FDR.CUTOFF, lfc=LOG.FC.CUTOFF)
  de2tags <- rownames(y)[as.logical(de2)]
  plotSmear(glm.list[[contrast]], de.tags=de2tags, 
            main=paste0("smear plot with GLM p.value < ", FDR.CUTOFF, " and LFC=", 
                        LOG.FC.CUTOFF, " cutoffs"))
  abline(h = c(-2, 2), col = "blue")
  
  # defining significant as FDR < FDR.CUTOFF and abs(logFC) > LOG.FC.CUTOFF 
  # add DGE.status column
  # UP, DOWN, NC
  tt$DGE.status <- "NC"
  # defining significant as FDR < FDR.CUTOFF and abs(logFC) > LOG.FC.CUTOFF 
  # for UP/DOWN
  if (nrow(tt[((tt$FDR < FDR.CUTOFF)&(tt$logFC > LOG.FC.CUTOFF)),]) > 0){
    tt[((tt$FDR < FDR.CUTOFF)&(tt$logFC > LOG.FC.CUTOFF)),]$DGE.status <- "UP"
  }
  if (nrow(tt[((tt$FDR < FDR.CUTOFF)&(tt$logFC < -LOG.FC.CUTOFF)),]) > 0){
    tt[((tt$FDR < FDR.CUTOFF)&(tt$logFC < -LOG.FC.CUTOFF)),]$DGE.status <- "DOWN"
  }
  
  # filtering/re-ordering
  column.order <- c(6:13, 1:5)
  tt <- tt[column.order]
  tt[is.na(tt)] <- ""
  
  write.table(tt, paste0(SUB.FOLDER, comparison, ".GLM.DGE.tsv"), sep="\t", quote=F, row.names=F)
  
  print (nrow(tt[((tt$FDR < FDR.CUTOFF)&(abs(tt$logFC) > LOG.FC.CUTOFF)),]))
  print (paste0("UP: ", nrow(tt[(tt$DGE.status=="UP"),])))
  print (paste0("DOWN: ", nrow(tt[(tt$DGE.status=="DOWN"),])))
  
  # defining significant as FDR < FDR.CUTOFF and abs(logFC) > LOG.FC.CUTOFF
  sigtt <- tt[((tt$FDR < FDR.CUTOFF)&(abs(tt$logFC) > LOG.FC.CUTOFF)),]
  
})

#[1] "glm.Veh.vs.Dec"
#[1] 306
#[1] "UP: 193"
#[1] "DOWN: 113"
#[1] "glm.Veh.vs.Tam"
#[1] 6049
#[1] "UP: 4638"
#[1] "DOWN: 1411"
#[1] "glm.Veh.vs.Comb"
#[1] 5708
#[1] "UP: 4222"
#[1] "DOWN: 1486"

# add annotation to TPM data and output
m <- match(rownames(tpm.data), gtf.anno$gene.id)

# assign the rownames(tt) as the gene_id as more specific with version number at end
# as originates from original gtf
tpm.data$gene.id <- rownames(tpm.data)
tpm.data$gene.symbol <- gtf.anno$gene.name[m]
tpm.data$chr <- gtf.anno$chr[m]
tpm.data$start <- gtf.anno$start[m]
tpm.data$end <- gtf.anno$end[m]
tpm.data$strand <- gtf.anno$strand[m]
tpm.data$gene.type <- gtf.anno$gene.type[m]

m <- match(gsub("\\.[0-9]*", "", rownames(tpm.data)), DT$ensembl.gene.ID)

tpm.data$description <- DT$description[m]
tpm.data$entrez.gene.id <- DT$entrez.gene.ID[m]

# get averages for each set of duplicates/triplicates
tpm.data$mean_Veh <- 
  round(rowMeans(tpm.data[,grep("GAR15.13D_Vehicle", colnames(tpm.data))]), 2)
tpm.data$mean_Tam <- 
  round(rowMeans(tpm.data[,grep("GAR15.13D_Tam", colnames(tpm.data))]), 2)
tpm.data$mean_Dec <- 
  round(rowMeans(tpm.data[,grep("GAR15.13D_Dec", colnames(tpm.data))]), 2)
tpm.data$mean_Comb <- 
  round(rowMeans(tpm.data[,grep("GAR15.13D_Comb", colnames(tpm.data))]), 2)


column.order <- c(17:27, 1:16)
tpm.data <- tpm.data[column.order]
tpm.data[is.na(tpm.data)] <- ""

write.table(tpm.data, paste0("~/Desktop/PROJECTS/PDX_Decitabine/GAR15-13D/hg38/RNA-seq_hg38/annotated_GAR15-13D_TPM_table_GLM.tsv"), 
            sep="\t", quote=F, row.names=F)

#### NEXT

###Create top tables per comparison

# Top table Veh vs Dec
# the default method used to adjust p-values for multiple testing is BH.
tt1 <- topTags(lrt.Veh.vs.Dec, n=nrow(counts))$table

m <- match(rownames(tt1), gtf.anno$gene.id)

# assign the rownames(tt) as the gene_id as more specific with version number at end
# as originates from original gtf
tt1$gene.id <- rownames(tt1)
tt1$gene.symbol <- gtf.anno$gene.name[m]
tt1$chr <- gtf.anno$chr[m]
tt1$start <- gtf.anno$start[m]
tt1$end <- gtf.anno$end[m]
tt1$strand <- gtf.anno$strand[m]
tt1$gene.type <- gtf.anno$gene.type[m]

m <- match(gsub("\\.[0-9]*", "", rownames(tt1)), DT$Ensembl.Gene.ID)

tt1$description <- DT$description[m]
tt1$entrez.gene.id <- DT$entrez.gene.ID[m]

# only keep chromosome names beginning with chr1..22, X, Y; remove patch chromosome assignments 
# like JH806587.1, JH806587.1 etc
tt1 <- tt1[grep("chr*", tt1$chr),]

# GOseq Veh vs Dec
lengths.gene <- read.csv(paste0("~/Desktop/PROJECTS/PDX_Decitabine/GAR15-13D/hg38/RNA-seq_hg38/tables/GAR15-13D_effective_length_table.csv"), row.names=1)
# get the average gene length for each row
lengths <- apply(lengths.gene, 1, mean)
bias.data <- lengths[rownames(tt1)]
names(bias.data) <- tt1$gene.symbol
bias.data <- bias.data[!duplicated(names(bias.data))]
if (length(names(bias.data[(names(bias.data) == "")])) > 0){
  bias.data <- bias.data[-which(names(bias.data)=="")]
}
bias.data <- bias.data[-which(bias.data==0)]
if (length(names(bias.data[(is.na(names(bias.data)))])) > 0){
  bias.data <- bias.data[-which(is.na(names(bias.data)))]
}
sigtt1 <- tt1[((tt1$FDR < FDR.CUTOFF)&(abs(tt1$logFC) > LOG.FC.CUTOFF)),]
comparison.UP <- sigtt1$gene.symbol[sigtt1$logFC > 0]
comparison.DOWN <- sigtt1$gene.symbol[sigtt1$logFC < 0]

comparison.UP.DE <- sapply(names(bias.data), function (x) as.numeric(x %in% comparison.UP))
comparison.DOWN.DE <- sapply(names(bias.data), function (x) as.numeric(x %in% comparison.DOWN))

write.table(tt1, "~/Desktop/PROJECTS/PDX_Decitabine/GAR15-13D/hg38/RNA-seq_hg38/DGE_output_defaultQL/GAR15-13D_Dec.vs.Veh_DGE.tsv", sep = "\t", quote = F)

table(comparison.UP.DE)
#comparison.UP.DE
#0     1 
#24527   170

table(comparison.DOWN.DE)

#comparison.DOWN.DE
#0     1 
#24647    50 
library(goseq)
library(magrittr)
library(dplyr)
library(ggplot2)
pwf1_up <- nullp(comparison.UP.DE,"hg38","geneSymbol")
pwf1_down <- nullp(comparison.DOWN.DE,"h38","geneSymbol")

#Using the Wallenius approximation
GO.wall_UP <- goseq(pwf1_up,"hg19","geneSymbol")


head(GO.wall_UP)
write.table(GO.wall_UP, paste0("GAR15-13D_Veh.vs.Dec.GOterms_UP.tsv"), sep="\t", quote=F, row.names=F)
##visualise Up Dec.vs.Veh
##plot TOP 10 UP
GO.wall_UP %>% 
  top_n(10, wt=-over_represented_pvalue) %>% 
  mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
  ggplot(aes(x=hitsPerc, 
             y=term, 
             colour=over_represented_pvalue, 
             size=numDEInCat)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="Hits (%)", y="GO term", colour="p value", size="Count")

##BP only
GO.wall_UP_BP <- goseq(pwf1_up,"hg19","geneSymbol",test.cats=c("GO:BP"))

GO.wall_UP_BP %>% 
  top_n(10, wt=-over_represented_pvalue) %>% 
  mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
  ggplot(aes(x=hitsPerc, 
             y=term, 
             colour=over_represented_pvalue, 
             size=numDEInCat)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="Hits (%)", y="GO term", colour="p value", size="Count", title = "Dec_GO.wall_UP_BP")

GO.wall_DOWN_BP <- goseq(pwf1_down,"hg19","geneSymbol",test.cats=c("GO:BP"))
GO.wall_DOWN_BP %>% 
  top_n(10, wt=-over_represented_pvalue) %>% 
  mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
  ggplot(aes(x=hitsPerc, 
             y=term, 
             colour=over_represented_pvalue, 
             size=numDEInCat)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="Hits (%)", y="GO term", colour="p value", size="Count", title = "Dec_GO.wall_DOWN_BP")



GO.wall_DOWN <- goseq(pwf1_down,"hg19","geneSymbol")
head(GO.wall_DOWN)
write.table(GO.wall_DOWN, paste0("GAR15-13D_Veh.vs.Dec.GOterms_DOWN.tsv"), sep="\t", quote=F, row.names=F)


##Visualise top 10 down
GO.wall_DOWN %>% 
  top_n(10, wt=-over_represented_pvalue) %>% 
  mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
  ggplot(aes(x=hitsPerc, 
             y=term, 
             colour=over_represented_pvalue, 
             size=numDEInCat)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="Hits (%)", y="GO term", colour="p value", size="Count")




#Using random sampling
GO.sampUP <- goseq(pwf1_up,"hg19","geneSymbol",method="Sampling",repcnt=1000)


#Compare both methods
plot(log10(GO.wall_UP[,2]), log10(GO.sampUP[match(GO.sampUP[,1],GO.wall_UP[,1]),2]),
     xlab="log10(Wallenius p-values)",ylab="log10(Sampling p-values)",
     xlim=c(-3,0))
abline(0,1,col=3,lty=2)

GO.sampDOWN <- goseq(pwf1_down,"hg19","geneSymbol",method="Sampling",repcnt=1000)







##Enriched GO terms
enriched.GO_UP <- GO.wall_UP$category[p.adjust(GO.wall_UP$over_represented_pvalue, method="BH")<.05]


##Heatmap - DEG in Veh vs Dec ordered by PValue Veh vs Dec (top 230)
logCPM <- cpm(y,prior.count = 2,log = TRUE) 
o <- order(lrt.Veh.vs.Dec$table$PValue)
logCPM <- logCPM[o[1:230],]
logCPM <- t(scale(t(logCPM)))
library(gplots)
col.pan <- colorpanel(100,"blue","white","red")
heatmap.2(logCPM,col=col.pan,Rowv = TRUE,scale = "none",trace = "none",dendrogram = "both",cexRow = 1, cexCol = 1.4,margin = c(10,9), lhei = c(2,10),lwid = c(2,6))



###Heatmap - DEG in Veh vs Tam ordered by PValue Veh vs Tam (top 300)
logCPM <- cpm(y,prior.count = 2,log = TRUE) 
o <- order(lrt.Veh.vs.Tam$table$PValue)
logCPM <- logCPM[o[1:300],]
logCPM <- t(scale(t(logCPM)))
library(gplots)
col.pan <- colorpanel(100,"blue","white","red")
heatmap.2(logCPM,col=col.pan,Rowv = TRUE,scale = "none",trace = "none",dendrogram = "both",cexRow = 1, cexCol = 1.4,margin = c(10,9), lhei = c(2,10),lwid = c(2,6))

###Heatmap - DEG in Veh vs Comb ordered by PValue Veh vs Comb (top 300)
logCPM <- cpm(y,prior.count = 2,log = TRUE) 
o <- order(lrt.Veh.vs.Comb$table$PValue)
logCPM <- logCPM[o[1:30],]
logCPM <- t(scale(t(logCPM)))
library(gplots)
col.pan <- colorpanel(100,"blue","white","red")
heatmap.2(logCPM,col=col.pan,Rowv = TRUE,scale = "none",trace = "none",dendrogram = "both",cexRow = 1, cexCol = 1.4,margin = c(10,9), lhei = c(2,10),lwid = c(2,6))


#### Top table Veh vs Tam
# the default method used to adjust p-values for multiple testing is BH.
tt2 <- topTags(lrt.Veh.vs.Tam, n=nrow(counts))$table

m <- match(rownames(tt2), gtf.anno$gene.id)

# assign the rownames(tt) as the gene_id as more specific with version number at end
# as originates from original gtf
tt2$gene.id <- rownames(tt2)
tt2$gene.symbol <- gtf.anno$gene.name[m]
tt2$chr <- gtf.anno$chr[m]
tt2$start <- gtf.anno$start[m]
tt2$end <- gtf.anno$end[m]
tt2$strand <- gtf.anno$strand[m]
tt2$gene.type <- gtf.anno$gene.type[m]

m <- match(gsub("\\.[0-9]*", "", rownames(tt2)), DT$Ensembl.Gene.ID)

tt2$description <- DT$description[m]
tt2$entrez.gene.id <- DT$entrez.gene.ID[m]

# only keep chromosome names beginning with chr1..22, X, Y; remove patch chromosome assignments 
# like JH806587.1, JH806587.1 etc
tt2 <- tt2[grep("chr*", tt1$chr),]

# GOseq Veh vs Dec
lengths.gene <- read.csv(paste0("/Volumes/Joanna_HD/RNA-seq_Level_2/GAR15-13D/tables/GAR15-13D_effective_length_table.csv"), row.names=1)
# get the average gene length for each row
lengths <- apply(lengths.gene, 1, mean)
bias.data <- lengths[rownames(tt2)]
names(bias.data) <- tt2$gene.symbol
bias.data <- bias.data[!duplicated(names(bias.data))]
if (length(names(bias.data[(names(bias.data) == "")])) > 0){
  bias.data <- bias.data[-which(names(bias.data)=="")]
}
bias.data <- bias.data[-which(bias.data==0)]
if (length(names(bias.data[(is.na(names(bias.data)))])) > 0){
  bias.data <- bias.data[-which(is.na(names(bias.data)))]
}
sigtt2 <- tt2[((tt2$FDR < FDR.CUTOFF)&(abs(tt2$logFC) > LOG.FC.CUTOFF)),]
comparison.UP <- sigtt2$gene.symbol[sigtt2$logFC > 0]
comparison.DOWN <- sigtt2$gene.symbol[sigtt2$logFC < 0]

comparison.UP.DE <- sapply(names(bias.data), function (x) as.numeric(x %in% comparison.UP))
comparison.DOWN.DE <- sapply(names(bias.data), function (x) as.numeric(x %in% comparison.DOWN))

table(comparison.UP.DE)
#comparison.UP.DE
#0     1 
#21807  2891 

table(comparison.DOWN.DE)

#comparison.DOWN.DE
#0     1 
#23953   745 

pwf1_up <- nullp(comparison.UP.DE,"hg19","geneSymbol")
pwf1_down <- nullp(comparison.DOWN.DE,"hg19","geneSymbol")

#Using the Wallenius approximation
GO.wall_UP <- goseq(pwf1_up,"hg19","geneSymbol")
head(GO.wall_UP)
write.table(GO.wall_UP, paste0("GAR15-13D_Veh.vs.Tam.GOterms_UP.tsv"), sep="\t", quote=F, row.names=F)

GO.wall_DOWN <- goseq(pwf1_down,"hg19","geneSymbol")
head(GO.wall_DOWN)
write.table(GO.wall_DOWN, paste0("GAR15-13D_Veh.vs.Tam.GOterms_DOWN.tsv"), sep="\t", quote=F, row.names=F)

#Using random sampling
GO.sampUP <- goseq(pwf1_up,"hg19","geneSymbol",method="Sampling",repcnt=1000)

#Compare both methods
plot(log10(GO.wall_UP[,2]), log10(GO.sampUP[match(GO.sampUP[,1],GO.wall_UP[,1]),2]),
     xlab="log10(Wallenius p-values)",ylab="log10(Sampling p-values)",
     xlim=c(-3,0))
abline(0,1,col=3,lty=2)

GO.sampDOWN <- goseq(pwf1_down,"hg19","geneSymbol",method="Sampling",repcnt=1000)

#Compare both methods
plot(log10(GO.wall_DOWN[,2]), log10(GO.sampDOWN[match(GO.sampDOWN[,1],GO.wall_DOWN[,1]),2]),
     xlab="log10(Wallenius p-values)",ylab="log10(Sampling p-values)",
     xlim=c(-3,0))
abline(0,1,col=3,lty=2)






## All genes from exp per cell line output

cpms <- cpm(y)
log.cpms <- log(cpms + 1)
col.order <- c(3,4,5,6,1,2)
#column.order <- c(1:15, 1:5)
gtf.anno2 <- gtf.anno[col.order]
write.table(gtf.anno2, paste0("GTF_anno_bed.bed"), sep="\t", quote=F, row.names=F)

exp_file <- "/Volumes/Joanna_HD/Caldon_Lab/RNAseq_MCF7_Garvan/LizCaldon_MCF7TAMRFASR_temporal_RNAseq/expr_per_cell_line_output/GENE_Liz_TAMR_expression_data.tsv"
temp <- read.table(exp_file, sep="\t", header=TRUE)
col.order <- c(3,4,5,6,1,2,9:11,13:18)
temp <- temp[col.order]
temp <- temp[grep("chr*", temp$chr),]
temp[is.na(temp)] <- ""
write.table(temp, paste0("ex_per_cell_anno_bed.bed"), sep="\t", quote=F, row.names=F)

test <- read.table("/Users/joanna/Desktop/PROJECTS/MCF7_Garvan_Liz_temporal/RNAseq/DMRs/exp_TAMR_hyper_heatmap_uniq.txt", header = F, row.names=1, col.names=0)
logdata <- log1p(test)

#### Heatmap genes at DMRs
##Heatmap
test <- read.table("Gene_Symbol_DMRs", header=TRUE)
#TPM <- read.table("/Volumes/Joanna_HD/Caldon_Lab/RNAseq_MCF7_Garvan/LizCaldon_MCF7TAMRFASR_temporal_RNAseq/expr_per_cell_line_output/GENE_Liz_TAMRF_expression_data.tsv", row.names=1, header=TRUE)

m <- match(rownames(logCPM), gtf.anno$gene.id)

# assign the rownames
logCPM$gene.id <- rownames(logCPM)
logCPM$gene.symbol <- gtf.anno$gene.name[m]
logCPM$chr <- gtf.anno$chr[m]
logCPM$start <- gtf.anno$start[m]
logCPM$end <- gtf.anno$end[m]
logCPM$strand <- gtf.anno$strand[m]
logCPM$gene.type <- gtf.anno$gene.type[m]


################# Continue GoSeq

# Top table Veh vs Tam
# the default method used to adjust p-values for multiple testing is BH.
tt2 <- topTags(lrt.Veh.vs.Tam, n=nrow(counts))$table

m <- match(rownames(tt2), gtf.anno$gene.id)

# assign the rownames(tt) as the gene_id as more specific with version number at end
# as originates from original gtf
tt2$gene.id <- rownames(tt2)
tt2$gene.symbol <- gtf.anno$gene.name[m]
tt2$chr <- gtf.anno$chr[m]
tt2$start <- gtf.anno$start[m]
tt2$end <- gtf.anno$end[m]
tt2$strand <- gtf.anno$strand[m]
tt2$gene.type <- gtf.anno$gene.type[m]

m <- match(gsub("\\.[0-9]*", "", rownames(tt2)), DT$Ensembl.Gene.ID)

tt2$description <- DT$description[m]
tt2$entrez.gene.id <- DT$entrez.gene.ID[m]

# only keep chromosome names beginning with chr1..22, X, Y; remove patch chromosome assignments 
# like JH806587.1, JH806587.1 etc
tt2 <- tt2[grep("chr*", tt2$chr),]

# GOseq Veh vs Tam
lengths.gene <- read.csv(paste0("/Volumes/Joanna_HD/RNA-seq_Level_2/GAR15-13D/tables/GAR15-13D_effective_length_table.csv"), row.names=1)
# get the average gene length for each row
lengths <- apply(lengths.gene, 1, mean)
bias.data <- lengths[rownames(tt2)]
names(bias.data) <- tt2$gene.symbol
bias.data <- bias.data[!duplicated(names(bias.data))]
if (length(names(bias.data[(names(bias.data) == "")])) > 0){
  bias.data <- bias.data[-which(names(bias.data)=="")]
}
bias.data <- bias.data[-which(bias.data==0)]
if (length(names(bias.data[(is.na(names(bias.data)))])) > 0){
  bias.data <- bias.data[-which(is.na(names(bias.data)))]
}
sigtt2 <- tt2[((tt2$FDR < FDR.CUTOFF)&(abs(tt2$logFC) > LOG.FC.CUTOFF)),]
comparison.UP <- sigtt2$gene.symbol[sigtt2$logFC > 0]
comparison.DOWN <- sigtt2$gene.symbol[sigtt2$logFC < 0]

comparison.UP.DE <- sapply(names(bias.data), function (x) as.numeric(x %in% comparison.UP))
comparison.DOWN.DE <- sapply(names(bias.data), function (x) as.numeric(x %in% comparison.DOWN))

table(comparison.UP.DE)
#comparison.UP.DE
#0     1 
#21807  2891 

table(comparison.DOWN.DE)

#comparison.DOWN.DE
#0     1 
#23953   745

pwf1_up <- nullp(comparison.UP.DE,"hg19","geneSymbol")
pwf1_down <- nullp(comparison.DOWN.DE,"hg19","geneSymbol")

#Using the Wallenius approximation
GO.wall_UP <- goseq(pwf1_up,"hg19","geneSymbol")

##visualise Up Tam.vs.Veh
##plot TOP 10 UP
GO.wall_UP %>% 
  top_n(10, wt=-over_represented_pvalue) %>% 
  mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
  ggplot(aes(x=hitsPerc, 
             y=term, 
             colour=over_represented_pvalue, 
             size=numDEInCat)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="Hits (%)", y="GO term", colour="p value", size="Count")



GO.wall_DOWN <- goseq(pwf1_down,"hg19","geneSymbol")


##Visualise top 10 down
GO.wall_DOWN %>% 
  top_n(10, wt=-over_represented_pvalue) %>% 
  mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
  ggplot(aes(x=hitsPerc, 
             y=term, 
             colour=over_represented_pvalue, 
             size=numDEInCat)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="Hits (%)", y="GO term", colour="p value", size="Count")

#
##BP only
GO.wall_UP_BP <- goseq(pwf1_up,"hg19","geneSymbol",test.cats=c("GO:BP"))

GO.wall_UP_BP %>% 
  top_n(10, wt=-over_represented_pvalue) %>% 
  mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
  ggplot(aes(x=hitsPerc, 
             y=term, 
             colour=over_represented_pvalue, 
             size=numDEInCat)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="Hits (%)", y="GO term", colour="p value", size="Count", title = "Tam_GO.wall_UP_BP")

GO.wall_DOWN_BP <- goseq(pwf1_down,"hg19","geneSymbol",test.cats=c("GO:BP"))
GO.wall_DOWN_BP %>% 
  top_n(10, wt=-over_represented_pvalue) %>% 
  mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
  ggplot(aes(x=hitsPerc, 
             y=term, 
             colour=over_represented_pvalue, 
             size=numDEInCat)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="Hits (%)", y="GO term", colour="p value", size="Count", title = "Tam_GO.wall_DOWN_BP")






# Top table Veh vs Comb
# the default method used to adjust p-values for multiple testing is BH.
tt3 <- topTags(lrt.Veh.vs.Comb, n=nrow(counts))$table

m <- match(rownames(tt3), gtf.anno$gene.id)

# assign the rownames(tt) as the gene_id as more specific with version number at end
# as originates from original gtf
tt3$gene.id <- rownames(tt3)
tt3$gene.symbol <- gtf.anno$gene.name[m]
tt3$chr <- gtf.anno$chr[m]
tt3$start <- gtf.anno$start[m]
tt3$end <- gtf.anno$end[m]
tt3$strand <- gtf.anno$strand[m]
tt3$gene.type <- gtf.anno$gene.type[m]

m <- match(gsub("\\.[0-9]*", "", rownames(tt3)), DT$Ensembl.Gene.ID)

tt3$description <- DT$description[m]
tt3$entrez.gene.id <- DT$entrez.gene.ID[m]

# only keep chromosome names beginning with chr1..22, X, Y; remove patch chromosome assignments 
# like JH806587.1, JH806587.1 etc
tt3 <- tt3[grep("chr*", tt3$chr),]

# GOseq Veh vs Comb
lengths.gene <- read.csv(paste0("/Volumes/Joanna_HD/RNA-seq_Level_2/GAR15-13D/tables/GAR15-13D_effective_length_table.csv"), row.names=1)
# get the average gene length for each row
lengths <- apply(lengths.gene, 1, mean)
bias.data <- lengths[rownames(tt3)]
names(bias.data) <- tt3$gene.symbol
bias.data <- bias.data[!duplicated(names(bias.data))]
if (length(names(bias.data[(names(bias.data) == "")])) > 0){
  bias.data <- bias.data[-which(names(bias.data)=="")]
}
bias.data <- bias.data[-which(bias.data==0)]
if (length(names(bias.data[(is.na(names(bias.data)))])) > 0){
  bias.data <- bias.data[-which(is.na(names(bias.data)))]
}
sigtt3 <- tt3[((tt3$FDR < FDR.CUTOFF)&(abs(tt3$logFC) > LOG.FC.CUTOFF)),]
comparison.UP <- sigtt3$gene.symbol[sigtt3$logFC > 0]
comparison.DOWN <- sigtt3$gene.symbol[sigtt3$logFC < 0]

comparison.UP.DE <- sapply(names(bias.data), function (x) as.numeric(x %in% comparison.UP))
comparison.DOWN.DE <- sapply(names(bias.data), function (x) as.numeric(x %in% comparison.DOWN))

table(comparison.UP.DE)
#comparison.UP.DE
#0     1 
#22186  2511  

table(comparison.DOWN.DE)

#comparison.DOWN.DE
#0     1 
#23861   836 


pwf1_up <- nullp(comparison.UP.DE,"hg19","geneSymbol")
pwf1_down <- nullp(comparison.DOWN.DE,"hg19","geneSymbol")

#Using the Wallenius approximation
GO.wall_UP <- goseq(pwf1_up,"hg19","geneSymbol")

##visualise Up Comb.vs.Veh
##plot TOP 10 UP
GO.wall_UP %>% 
  top_n(10, wt=-over_represented_pvalue) %>% 
  mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
  ggplot(aes(x=hitsPerc, 
             y=term, 
             colour=over_represented_pvalue, 
             size=numDEInCat)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="Hits (%)", y="GO term", colour="p value", size="Count")



GO.wall_DOWN <- goseq(pwf1_down,"hg19","geneSymbol")


##Visualise top 10 down
GO.wall_DOWN %>% 
  top_n(10, wt=-over_represented_pvalue) %>% 
  mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
  ggplot(aes(x=hitsPerc, 
             y=term, 
             colour=over_represented_pvalue, 
             size=numDEInCat)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="Hits (%)", y="GO term", colour="p value", size="Count")

####BP only
GO.wall_UP_BP <- goseq(pwf1_up,"hg19","geneSymbol",test.cats=c("GO:BP"))

GO.wall_UP_BP %>% 
  top_n(10, wt=-over_represented_pvalue) %>% 
  mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
  ggplot(aes(x=hitsPerc, 
             y=term, 
             colour=over_represented_pvalue, 
             size=numDEInCat)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="Hits (%)", y="GO term", colour="p value", size="Count", title = "Comb_GO.wall_UP_BP")

GO.wall_DOWN_BP <- goseq(pwf1_down,"hg19","geneSymbol",test.cats=c("GO:BP"))
GO.wall_DOWN_BP %>% 
  top_n(10, wt=-over_represented_pvalue) %>% 
  mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
  ggplot(aes(x=hitsPerc, 
             y=term, 
             colour=over_represented_pvalue, 
             size=numDEInCat)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="Hits (%)", y="GO term", colour="p value", size="Count", title = "Comb_GO.wall_DOWN_BP")

