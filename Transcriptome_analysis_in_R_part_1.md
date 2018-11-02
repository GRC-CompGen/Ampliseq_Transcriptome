---
title: "Transcriptome Analysis"
author: "Nicholas Harvey"
date: "November 2, 2018"
output:
  html_document: default
  pdf_document: default
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

```{R} 
require("BioCinstaller")
```

###Downloading relevant Packages

This script requires R V3.5.1
BiocInstaller::biocLite(c("tximport", "readr", "edgeR", "Homo.sapiens", "biomaRt", "qusage", "limma", "Glimma", "ChAMP"))


##info

This document outlines the necessary coding and steps to upload transcriptome Salmon Quant files into R using the TxImport package
note: this setup requires that all quant files are in directory stucture `main_folder/2-Input/quants/sample_name_quant/quant.sf`


#### Set working directory 
```{r working directory, include=FALSE}
setwd("H:/Transcriptome_R")
dir.create("3-QC")
dir.create("4-Output")
options(digits=3)
getwd()
```

####load required packages
```{r, include=FALSE}
library(tximport)
library(readr)
```

##Using the tximport package to load the salmon quant files into R
####tximport of pre-constructed tx2gene table - created from bioMart Ensembl site
```{R transcript to gene file, include=FALSE}
tx2gene <- read.csv("H:/Transcriptome_R/2-Input/Homo_sapiens.GRCh38.91_tx2gene.csv")
head(tx2gene, 5)
```

####create links to salmon quant folders for import
```{R}
folder <- c("H:/Transcriptome_R/2-Input/quants")
salmon.dir <- as.matrix(read.csv(file="H:/Transcriptome_R/2-Input/quant_filenames.csv", sep=",", header=F))

salmon.files <- file.path(folder, salmon.dir, "quant.sf")
names(salmon.files) <- as.matrix(read.csv(file="H:/Transcriptome_R/2-Input/names.csv", sep=",", header=F))
file.exists(salmon.files)
all(file.exists(salmon.files))
```


####tximport
countsFromAbundance default is "no" and just returns counts not Transcript Per Million (TPM). scaledTPM is TPM scaled up to library size while lengthScaledTPM first multiplies TPM by feature length and then scales up to library size. Both are then quantities that are on the same scale as original counts,except no longer correlated with feature length across samples.

```{R}
require(tximport)
txi <- tximport(salmon.files, type="salmon", tx2gene=tx2gene, countsFromAbundance="lengthScaledTPM") 
names(txi)
head(txi$counts, 3)
write.csv(txi, file = "H:/Transcriptome_R/4-Output/geneSMART_tximport_matrix.csv")
```

#### prepare data
EdgeR is used to read the count data into R
```{R}
library(edgeR)
library(Homo.sapiens)

y <- DGEList(txi$counts)
dim(y)
```

#### Read in sample information table

In this case sample information is coming from the GeneSMART phenotypes table 
Ideally this table should contain the necessary variables to both perform a reasonable PCA plot and account for the most variance in the data set
```{R}
csvfile <- file.path("H:/Transcriptome_R/2-Input/geneSMART_sample_table.csv")
sampleTable <- read.csv(csvfile, row.names=1)
y$samples$name <- sampleTable$name
y$samples$filename <- sampleTable$filename
y$samples$sampleID <- sampleTable$sampleID
y$samples$subject <- sampleTable$subject
y$samples$timepoint <- sampleTable$timepoint
y$samples$id <- sampleTable$id
y$samples$run_date <- sampleTable$run_date
y$samples
```


#### Add gene annotation from biomaRt package
```{R}
geneid <- rownames(y)

library("biomaRt")
geneid <- rownames(y)
ensembl91 <- useMart(host="dec2017.archive.ensembl.org", 
                     biomart="ENSEMBL_MART_ENSEMBL", 
                     dataset="hsapiens_gene_ensembl")
attributes <- listAttributes(ensembl91)
attributes[1:20,]

genes <- select(ensembl91, keys=geneid, keytype="ensembl_gene_id",
                columns=c("ensembl_gene_id", "external_gene_name",
                          "description","entrezgene","chromosome_name","gene_biotype"))

(colnames(genes) <- c("ENSEMBL","SYMBOL","GENENAME","GENEID","TXCHROM","BIOTYPE"))
```

#### Remove duplicated genes
```{R}
genes <- genes[!duplicated(genes$ENSEMBL),]
y$genes <- genes
head(genes,5)
dim(genes)
```

#### Create barplot of library sizes
creates barplot of library sizes to qualitatively assess outliers (the abline functions are for Ampliseq Transcriptome library sizes and should be adjusted for other experiments)
```{R}
png("H:/Transcriptome_R/3-QC/Barplot of library sizes.png", width = 90, height = 30, units = 'cm', res = 300)
col <- as.numeric(y$sample$timepoint)
dt <- colSums((y$counts)*1e-6)
barplot(dt, names=colnames(dt), col=c("lightsalmon","lightcoral","steelblue1","cornsilk")[col], las=2, cex.names=0.8)
abline(h=5,col="black",lty=3)
abline(h=10,col="black",lty=3)
abline(h=15,col="black",lty=3)
title(main="Barplot of library sizes",ylab="Library size (millions)")
dev.off()
```

#### Summary of total mapped read counts
```{R} 
summary(dt)
```
#### Add all reads and average count for each gene in all samples, sorted from highest to lowest
```{R}
CountMeans <- rowMeans(y$counts)
TotalCounts <- sum(CountMeans)


PercentReads <- (CountMeans/TotalCounts)*100
AvgCounts <- rowMeans(y$counts)
Symbol <- y$genes$SYMBOL
GeneName <- y$genes$GENENAME
AvgCountsTable <- (data.frame(Symbol,GeneName,AvgCounts,PercentReads))
AvgCountsTable <- AvgCountsTable[order(-AvgCountsTable$AvgCounts),]
head(AvgCountsTable,20)
write.csv(AvgCountsTable, file ="H:/Transcriptome_R/4-Output/Average_Counts.csv")

top5 <- ((sum(head(AvgCountsTable, 5)$AvgCounts))/TotalCounts)*100
top10 <- ((sum(head(AvgCountsTable, 10)$AvgCounts))/TotalCounts)*100
top25 <- ((sum(head(AvgCountsTable, 25)$AvgCounts))/TotalCounts)*100

top5; top10; top25
```

#### Summary of the lengthscaledTPM count data
```{R}
dim(y)
summary(rowMeans(y$counts))
```
#### Create histogram of count distributions across samples
 Note that here we have added a pseudocount of 0.25 to all values to prevent logging zero values and therefore O CPM = lcpm -6.9

```{R}
lcpm.AvgCounts <- cpm(AvgCounts, log = TRUE)
png("H:/Transcriptome_R/3-QC/Read Count Distribution-B4.png", width = 45, height = 25, units = 'cm', res = 600)
par(mar=c(5,6,4,1)+.1)
hist(lcpm.AvgCounts, col="salmon", border="salmon",
     cex.lab=2.5, cex.main=3, cex.axis=2,
     xlab="Median log2-CPM", ylab="No. of Transcripts",
     breaks=100, xlim=c(-10,20), main ="Read Count Distribution (before)")
dev.off()

table(AvgCounts>=1)
table(AvgCounts>=10)
table(AvgCounts>=100)
table(AvgCounts>=1000)
table(AvgCounts>=10000)
table(AvgCounts>=100000)
```


Filter unexpressed and very low expressed genes and Count genes with zero counts across all 211 samples
```{R}
table(rowSums(y$counts==0)==211)
```
## FILTERING Based on median log2 CPM
We expect approximately 9K-12K genes to remain after filtering https://www.biostars.org/p/211954/
```{R}
median_cpm <- apply(cpm(y), 1, median)
expr_cutoff <- 0.5 # in cpm
sum(median_cpm > expr_cutoff)
y.Filt <- y[median_cpm > expr_cutoff, ]
lcpm.Filt <- cpm(y.Filt$counts, log = TRUE)

log.cutoff <- log2(expr_cutoff)

summary(rowMeans(y.Filt$counts))
```


#### Density of count values (code modified from 'RNAseq 1-2-3')
```{R}
png("H:/Transcriptome_R/3-QC/Density of count values.png", width = 10, height = 20, units = 'cm', res = 600)
nsamples <- ncol(y)
col <- rainbow(nsamples)
par(mfrow=c(2,1))
lcpm.Raw <- cpm(y$counts, log=TRUE)
plot(density(lcpm.Raw[,1]), col=col[1], xlim=c(-10,20), ylim=c(0,0.3), main="", xlab="")
for (i in 2:nsamples){
  den <-density(lcpm.Raw[,i])
  lines(den$x, den$y, col=col[i])
}
abline(v=log.cutoff, col="red", lwd=1, lty=2, main="")
title("Raw data",xlab="log2-CPM")

lcpm.Filt <- cpm(y.Filt$counts, log=TRUE)
plot(density(lcpm.Filt[,1]), col=col[1], xlim=c(-10,20), ylim=c(0,0.3), main="", xlab="")
for (i in 2:nsamples){
  den <-density(lcpm.Filt[,i])
  lines(den$x, den$y, col=col[i])
}
abline(v=log.cutoff, col="red", lwd=1, lty=2, main="")
title("Filtered data (median CPM > 0.5)",xlab="log2-CPM")
dev.off()
```

## Histogram of count distribution
pseudocount of 0.25 added to values to prevent logging zero values therefore O CPM = lcpm -6.9
```{R}
AvgCounts <- rowMeans(y.Filt$counts)
lcpm.AvgCounts <- cpm(AvgCounts, log = TRUE)
png("H:/Transcriptome_R/3-QC/Read Count Distribution-AFTER.png", width = 45, height = 25, units = 'cm', res = 600)
par(mar=c(5,6,4,1)+.1)
hist(lcpm.AvgCounts, col="salmon", border="salmon",
     cex.lab=2.5, cex.main=3, cex.axis=2,
     xlab="Median log2-CPM", ylab="No. of Transcripts",
     breaks=100, xlim=c(-10,20), main ="Read Count Distribution (after)")
dev.off()
```

####Heatmap of samples
```{R}
png("3-QC/Sample heatmap.png", width = 60, height = 60, units = 'cm', res = 300)
par(mfrow=c(1,2))
lcpm.Raw <- cpm(y$counts, log = TRUE)
heatmap(cor(lcpm.Raw))
title("Raw data")
lcpm.Filt <- cpm(y.Filt$counts, log = TRUE)
heatmap(cor(lcpm.Filt))
title("Filtered data")
dev.off()
```

##  limma

the Limma package is used firstly for normalisation of the count data and then comparison of upregulated/downregulated gene expression
```{R}
library(limma)

y.Norm <- calcNormFactors(y.Filt, method="TMM")
y.Norm$samples$norm.factors
summary(y.Norm$samples$norm.factors)
```

####boxplot BEFORE normalisation
```{R}
png("3-QC/Boxplots.png", width = 90, height = 45, units = 'cm', res = 300)
par(mfrow=c(2,1))
lcpm.Filt <- cpm(y.Filt, log=TRUE)
col <- as.numeric(y$sample$timepoint)
boxplot(lcpm.Filt,las=2,
        col=c("lightsalmon","lightcoral","steelblue1","cornsilk")[col],main="")
abline(h=median(lcpm.Filt),col="black",lty=3)
title(main="Unnormalised data",ylab="log-counts")
```
####boxplot AFTER normalisation
```{R}
lcpm.Norm <- cpm(y.Norm, log=TRUE)
boxplot(lcpm.Norm,las=2,
        col=c("lightsalmon","lightcoral","steelblue1","cornsilk")[col],main="")
abline(h=median(lcpm.Norm),col="black",lty=3)
title(main="TMM-Normalised data",ylab="log-counts")
dev.off()
```


## UNSUPERVISED CLUSTERING OF SAMPLES
Creation of MDS plots for visualisation of data pre and post normalisation
```{R}
dir.create("5-Glimma")
library(Glimma)
glMDSPlot(y.Norm, top=500, labels=y.Norm$samples$subject,
          groups=y.Norm$samples[,c(1:9)], launch=TRUE,
          path="5-Glimma", folder="glMDSPlot", html="MDS_plot")
```

 MODIFY CODE at this point to ensure accurate running of the script

#### MDS plot of sample names
```{R}
png("3-QC/MDS Plot_names.png", width = 15, height = 30, units = 'cm', res = 600)
par(mfrow=c(2,1))
plotMDS(y.Norm, top=500, cex=0.8, labels=y.Norm$samples$subject, col=as.numeric(y.Norm$samples$timepoint), main="MDS Plot (SampleID)")
legend("topright", legend=c("PRE", "P0", "P3", "P4"), cex=0.8, col=1:16, pch=16)

plotMDS(y.Norm, top=500, cex=1, pch=21, col="white", bg=as.numeric(y.Norm$samples$timepoint), main="MDS Plot (Symbol)")
legend("topright", legend=c("PRE", "P0", "P3", "P4"), cex=0.8, col=1:16, pch=16)
dev.off()
```
PCH function guidelines: pch=xx; 0= open square; 1= open circle; 15= solid square, 16=solid circle, 21=filled circle/border


## Mean Difference Plots
Library size-adjusted log-fold change between two libraries (the difference) against the average log-expression across those libraries (the mean).
The following command produces an MD plot that compares MD plot before and after normalisation

```{R}
png("3-QC/MDPlots-1LI00.png", width=11, height=20, units='cm', res=200)
par(mfrow=c(2,1))
plotMD(y.Filt, column=1, main="First sample (raw)"); abline(h=0, col="red", lty=2, lwd=2)
plotMD(y.Norm, column=1, main="First sample (TMM-normalised)"); abline(h=0, col="red", lty=2, lwd=2)
dev.off()
```
#Important! Limma design matrix
At this point all transcripts have been mapped to genes, zero counts have been excluded from further analysis, and all remaining counts have been normalised across samples. The following design matrix is the most important step of this script because it is when we actually tell the program what to compare for gene expression changes. (the design below considers the "PRE" timepoint as the baseline (which is correct but only considers the question "what changes between timepoints for all samples whether they are responders or not"))

```{R}
group <- y$samples$timepoint
design <- model.matrix(~0+group)
colnames(design) <- levels(group)
design
```
## Paired design using DuplicateCorrelation
https://support.bioconductor.org/p/52920/  https://support.bioconductor.org/p/59700/
```{R}
Group <- y.Norm$samples$timepoint
design <- model.matrix(~0 + Group)

v <- voom(y.Norm, design)
colnames(design) <- c("P0","P3","P4","PRE")
corfit <- duplicateCorrelation(v, design, block=y.Norm$samples$subject)
v <- voom(y.Norm, design, block=y.Norm$samples$subject, correlation=corfit$consensus)

save(v, file="H:/Transcriptome_R/4-Output/v.rda")
```

## Use the SVD functon in ChAMP to assess which covariates correlate with the top principle components
```{R}
library(ChAMP)

pd <- v$targets[,c("sampleID","subject",
                   "timepoint","id","run_date")]
colnames(pd)[1] <- "Sample_Names"
champ.SVD(beta=v$E, pd=pd, PDFplot=TRUE, Rplot=FALSE, resultsDir="./3-QC/")

glMDSPlot(v, top=500, labels=v$targets$subject,
          groups=v$targets[,c(1:9)], launch=TRUE,
          path="5-Glimma", folder="glMDSPlot", html="MDS_plot_voom")
```
#### MD plot before and after normalisation
```{R}
png("H:/Transcriptome_R/3-QC/MDPlots-sample54.png", width=10, height=20, units='cm', res=300)
par(mfrow=c(3,1))
plotMD(y.Filt, column=54, main="sample 54 (raw)"); abline(h=0, col="red", lty=2, lwd=2)
plotMD(y.Norm, column=54, main="sample 54 (TMM-normalised)"); abline(h=0, col="red", lty=2, lwd=2)
plotMD(v, column=54, main="sample 54 (voom)"); abline(h=0, col="red", lty=2, lwd=2)
dev.off()
```

# Linear model fitting (eBayes or TREAT)
Actual design matrix (i.e. every timepoint - PRE timepoint)
```{R}
lfc <- log2(1.1) # use for TREAT only

fit <- lmFit(v, design, block=y.Norm$samples$subject, correlation=corfit$consensus)
cont.matrix = makeContrasts("P0-PRE" = P0-PRE,
                            "P3-PRE" = P3-PRE,
                            "P4-PRE" = P4-PRE,
                            levels=design)
fit2 = contrasts.fit(fit, cont.matrix)
fit2 = eBayes(fit2)
#fit2 <- treat(fit2, lfc=lfc)
```
#### SA plots
```{R}
png("H:/Transcriptome_R/3-QC/SA Plot.png", width = 30, height = 15, units = 'cm', res = 600)
par(mfrow=c(1,2))
v <- voom(y.Norm, design, plot=TRUE)
plotSA(fit2, main="Final model: Mean−variance trend")
dev.off()
```

## Differential Expression and total number of signifcant genes
```{R}
results <- decideTests(fit2)
summary(results)
vennCounts(results)

ngenes <-  which(results[,1] != 0 | results[,2] != 0 | results[,3] != 0)
ngenes <- length(ngenes)
```
#### Venn Diagram of design matrix and genes in common
```{R}
png("H:/Transcriptome_R/4-Output/Venn Diagram.png", width = 30, height = 15, units = 'cm', res = 600)
par(mfrow=c(1,2))
vennDiagram(results,include=c("both"), circle.col=c("blue","yellow","green","orange"), counts.col=c("blue3"), cex=c(1,0.8,0.8))
vennDiagram(results,include=c("up","down"), circle.col=c("blue","yellow","green","orange"), counts.col=c("red","green3"), cex=c(1,0.7,0.7))
mtext("geneSMART (eBayes)", side = 3, line = -2, outer = TRUE, cex=1.8)
mtext(paste0("(", ngenes," genes)"), side = 3, line = -3.5, outer = TRUE, cex=1.5)
dev.off()
```
#### MD Plots for each matrix component 
Note: these plots are interactive 
```{R}
png("H:/Transcriptome_R/4-Output/MDplot-Average.png", width = 20, height = 30, units = 'cm', res = 600)
par(mfrow=c(3,1))
plotMD(fit2, coef=1, status=results[,"P0-PRE"], values = c(-1, 1), main=colnames(fit2)[1])
plotMD(fit2, coef=2, status=results[,"P3-PRE"], values = c(-1, 1), main=colnames(fit2)[2])
plotMD(fit2, coef=3, status=results[,"P4-PRE"], values = c(-1, 1), main=colnames(fit2)[3])
dev.off()

glMDPlot(fit2, coef=1, status=results[,"P0-PRE"], main=colnames(fit2)[1],
         side.main="SYMBOL", counts=y.Norm$counts, groups=v$targets$timepoint,
         launch=TRUE, path="5-Glimma", folder="glMDPlots", html="P0-PRE")

glMDPlot(fit2, coef=2, status=results[,"P3-PRE"], main=colnames(fit2)[2],
         side.main="SYMBOL", counts=y.Norm$counts, groups=v$targets$timepoint,
         launch=TRUE, path="5-Glimma", folder="glMDPlots", html="P3-PRE")

glMDPlot(fit2, coef=3, status=results[,"P4-PRE"], main=colnames(fit2)[3],
         side.main="SYMBOL", counts=y.Norm$counts, groups=v$targets$timepoint,
         launch=TRUE, path="5-Glimma", folder="glMDPlots", html="P4-PRE")
```

####Volcano Plots
https://bioconductor.org/packages/release/bioc/vignettes/Glimma/inst/doc/Glimma.pdf
```{R}
png("H:/Transcriptome_R/4-Output/Volcano_plots.png", width = 20, height = 30, units = 'cm', res = 600)
par(mfrow=c(3,1))
cutoff = -log10(0.05)
volcanoplot(fit2,coef=1,highlight=25,names=v$genes$SYMBOL,main=colnames(fit2)[1])
volcanoplot(fit2,coef=2,highlight=25,names=v$genes$SYMBOL,main=colnames(fit2)[2])
volcanoplot(fit2,coef=3,highlight=25,names=v$genes$SYMBOL,main=colnames(fit2)[3])
dev.off()


glXYPlot(x=fit2$coef[,"P0-PRE"], y=fit2$lod[,"P0-PRE"], xlab="logFC", ylab="logodds",
         status=results[,"P0-PRE"], main=colnames(fit2)[1], side.main="SYMBOL",
         anno=fit2$genes, counts=y.Norm$counts, groups=v$targets$timepoint,
         path="5-Glimma", folder="glVolcano", html="P0-PRE")

glXYPlot(x=fit2$coef[,"P3-PRE"], y=fit2$lod[,"P3-PRE"], xlab="logFC", ylab="logodds",
         status=results[,"P3-PRE"], main=colnames(fit2)[2], side.main="SYMBOL",
         anno=fit2$genes, counts=y.Norm$counts, groups=v$targets$timepoint,
         path="5-Glimma", folder="glVolcano", html="P3-PRE")

glXYPlot(x=fit2$coef[,"P4-PRE"], y=fit2$lod[,"P4-PRE"], xlab="logFC", ylab="logodds",
         status=results[,"P4-PRE"], main=colnames(fit2)[3], side.main="SYMBOL",
         anno=fit2$genes, counts=y.Norm$counts, groups=v$targets$timepoint,
         path="5-Glimma", folder="glVolcano", html="P4-PRE")
```


## which genes respond to exercise - P0 relative to PRE?
```{R}
r1 <- topTable(fit2, adjust="BH", coef=1, n=Inf)
write.csv(r1, file="H:/Transcriptome_R/4-Output/topTable_P0-PRE.csv")
sig <- r1$adj.P.Val <0.05
cat("No.Sig.Genes.P0-PRE:", length(which(sig==1)))
```
## which genes respond to exercise - P3 relative to PRE?
```{R}
r2 <- topTable(fit2, adjust="BH", coef=2, n=Inf)
write.csv(r2, file="H:/Transcriptome_R/4-Output/topTable_P3-PRE.csv")
sig <- r2$adj.P.Val <0.05
cat("No.Sig.Genes.P3-PRE:", length(which(sig==1)))
```
## which genes respond to exercise - P4 relative to PRE?
```{R}
r3 <- topTable(fit2, adjust="BH", coef=3, n=Inf)
write.csv(r3, file="H:/Transcriptome_R/4-Output/topTable_P4-PRE.csv")
sig <- r3$adj.P.Val <0.05
cat("No.Sig.Genes.P4-PRE:", length(which(sig==1)))

sessionInfo()
```
END