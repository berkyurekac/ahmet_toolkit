# this script was created by Ahmet Can Berkyurek on 10/04/2020, to correlate mRNA and smallRNa expression in bins! 
library(gplots)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library("genefilter")
library(pheatmap)
library(SummarizedExperiment)
library(gdata)
library(dplyr)
library(GenomicAlignments)
library(dplyr)
library(tidyr)
library(scales)
library(reshape2)
library(DESeq2)

# load featureCounts count.txt file. before loading, delete the first two rows for clarity.!!!textwrangler should help...
setwd("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/Figure_6/smallRNA_bins")


########
########
#ribo libraries mRNA
######

countspolyA <- read.delim("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/RNA-seq_ribozero/counts_Lisa-ribozero.txt", header = TRUE)
head(countspolyA)

colnames(countspolyA) <- c("genes", "chr", "start", "end", "starnd", "length", "wt1x", "wt2x", "wt3x", "mj261_1x", "mj261_2x", "mj261_3x")
head(countspolyA)

#rRNA <- read.delim("rRNAs.txt", header = TRUE)
#histone <- read.delim("histone_genes.txt", header = TRUE)

#countspolyA <- subset(countspolyA, !(countspolyA$genes) %in% rRNA$Gene.stable.ID)
#countspolyA <- subset(countspolyA, !(countspolyA$genes) %in% histone$Gene.stable.ID)

# quantile normalisation for every 25 percent

countspolyA <- countspolyA[order(countspolyA$wt1x),]

b <- c(rep("A",9348), rep("B",9348), rep("C",9348), rep("D",9348), rep("E",9348))

countspolyA$label <- b
head(countspolyA)

cluster1 <- countspolyA[countspolyA$label == "A",]
cluster2 <- countspolyA[countspolyA$label == "B",]
cluster3 <- countspolyA[countspolyA$label == "C",]
cluster4 <- countspolyA[countspolyA$label == "D",]
cluster5 <- countspolyA[countspolyA$label == "E",]

#normalize cluster1
cluster1[,7:12] <- cluster1[,7:12]/((cluster1$length)*0.001)

head(cluster1)

cluster1$wt1x <- cluster1$wt1x/(sum(cluster1$wt1x)/1e6)
cluster1$wt2x <- cluster1$wt2x/(sum(cluster1$wt2x)/1e6)
cluster1$wt3x <- cluster1$wt3x/(sum(cluster1$wt3x)/1e6)

cluster1$mj261_1x <- cluster1$mj261_1x/(sum(cluster1$mj261_1x)/1e6)
cluster1$mj261_2x <- cluster1$mj261_2x/(sum(cluster1$mj261_2x)/1e6)
cluster1$mj261_3x <- cluster1$mj261_3x/(sum(cluster1$mj261_3x)/1e6)

head(cluster1)

#TPM normalization Divide the read counts by the length of each gene in kilobases. This gives you reads per kilobase (RPK).
#Count up all the RPK values in a sample and divide this number by 1,000,000. This is your “per million” scaling factor.
#Divide the RPK values by the “per million” scaling factor. This gives you TPM.

countspolyA[,7:12] <- countspolyA[,7:12]/((countspolyA$length)*0.001)

head(countspolyA)

countspolyA$wt1x <- countspolyA$wt1x/(sum(countspolyA$wt1x)/1e6)
countspolyA$wt2x <- countspolyA$wt2x/(sum(countspolyA$wt2x)/1e6)
countspolyA$wt3x <- countspolyA$wt3x/(sum(countspolyA$wt3x)/1e6)

countspolyA$mj261_1x <- countspolyA$mj261_1x/(sum(countspolyA$mj261_1x)/1e6)
countspolyA$mj261_2x <- countspolyA$mj261_2x/(sum(countspolyA$mj261_2x)/1e6)
countspolyA$mj261_3x <- countspolyA$mj261_3x/(sum(countspolyA$mj261_3x)/1e6)

head(countspolyA)


#smallRNA libraries 

countssmall <- read.delim("mj261_smallRNA_total_ce11_antisense.txt", header = TRUE)
head(countssmall)

colnames(countssmall) <- c("genes", "chr", "start", "end", "starnd", "length", "wt1", "wt2", "wt3", "mj261_1", "mj261_2", "mj261_3", "rescue1", "rescue2", "rescue3")
head(countssmall)

countssmall <- subset(countssmall, !(countssmall$genes) %in% rRNA$Gene.stable.ID)
countssmall <- subset(countssmall, !(countssmall$genes) %in% histone$Gene.stable.ID)


#TPM normalization Divide the read counts by the length of each gene in kilobases. This gives you reads per kilobase (RPK).
#Count up all the RPK values in a sample and divide this number by 1,000,000. This is your “per million” scaling factor.
#Divide the RPK values by the “per million” scaling factor. This gives you TPM.

countssmall[,7:15] <- countssmall[,7:15]/((countssmall$length)*0.001)

head(countssmall)

countssmall$wt1 <- countssmall$wt1/(sum(countssmall$wt1)/1e6)
countssmall$wt2 <- countssmall$wt2/(sum(countssmall$wt2)/1e6)
countspolyA$wt3 <- countssmall$wt3/(sum(countssmall$wt3)/1e6)

countssmall$mj261_1 <- countssmall$mj261_1/(sum(countssmall$mj261_1)/1e6)
countssmall$mj261_2 <- countssmall$mj261_2/(sum(countssmall$mj261_2)/1e6)
countssmall$mj261_3 <- countssmall$mj261_3/(sum(countssmall$mj261_3)/1e6)

countssmall$rescue1 <- countssmall$rescue1/(sum(countssmall$rescue1)/1e6)
countssmall$rescue2 <- countssmall$rescue2/(sum(countssmall$rescue2)/1e6)
countssmall$rescue3 <- countssmall$rescue3/(sum(countssmall$rescue3)/1e6)

head(countssmall)


# combine polyA mRNA and small RNA - first cluster with small RNA small to high, then label rows for bining!

polyAmerged <- data.frame(countspolyA[1:12],countssmall[7:15])
head(polyAmerged)

polyAmerged$mRNAwtmean <- rowMeans(polyAmerged[,7:9])
polyAmerged$mRNAmj261mean <- rowMeans(polyAmerged[,10:12])

polyAmerged$smallwtmean <- rowMeans(polyAmerged[,13:15])
polyAmerged$smallmj261mean <- rowMeans(polyAmerged[,16:18])
polyAmerged$smallrescuemean <- rowMeans(polyAmerged[,19:21])


head(polyAmerged)

data <- polyAmerged[,c(1,7:26)]
head(data)
# order data in smallRNA mutant increasing order
data <- data[order(data$smallwtmean),]
head(data)


b <- c(rep("A",2332), rep("B",2332), rep("C",2332), rep("D",2332), rep("E",2332), rep("F",2332), rep("G",2332),rep("H",2332),rep("I",2332),rep("J",2332),rep("K",2332),rep("L",2332),rep("M",2332),rep("N",2332),rep("O",2332),rep("P",2332),rep("R",2332),rep("S",2332),rep("T",2332),rep("U",2347))

b <- as.vector(b)

data$label <- b

head(data)

sensor <- subset(data, data$genes == "piRNASensor")

p <- ggplot(data, aes(x=data$label, y=log10(data$mRNAwtmean+1))) + geom_boxplot(outlier.colour="Grey", outlier.shape=NA, outlier.size=2.5, notch=FALSE, coef = 6) + geom_boxplot(data=sensor, aes(x=sensor$label, y=log10(sensor$mRNAwtmean+1)), color="purple", size=2) + theme_classic() + scale_y_continuous(limits = c(0,2.5))+ xlab("smallRNA log10(TPM+1) bins") + ylab("wild type polyA totalRNA log10(TPM+1)")

ggsave("ribo_wildtype_mRNAsmallRNA.pdf", p, height = 10, width = 15, useDingbats=FALSE)

g <- p + geom_boxplot(data=data, aes(x=data$label, y=log10(data$mRNAmj261mean+1)), color="Red",outlier.shape=NA,linetype="dashed", outlier.size = 1, alpha = 0.5) + geom_boxplot(data=sensor, aes(x=sensor$label, y=log10(sensor$mRNAmj261mean+1)), color="green", size=2)

ggsave("ribo_wildtype_mutant_mRNAsmallRNA.pdf", g, height = 10, width = 15, useDingbats=FALSE)


data <- data[order(data$smallmj261mean),]
head(data)

p2 <- ggplot(data, aes(x=data$label, y=log10(data$mRNAwtmean+1))) + geom_boxplot(outlier.colour="Grey", outlier.shape=NA, outlier.size=2, notch=FALSE, coef = 6) + theme_classic() + scale_y_continuous(limits = c(0,4))+ xlab("smallRNA log10(TPM+1) bins") + ylab("rpb-9(mj261)V polyA totalRNA log10(TPM+1)")


g2 <- p2 +  geom_boxplot(data=data, aes(x=data$label, y=log10(data$mRNAmj261mean+1)), color="Red",outlier.shape=NA,linetype="dashed", outlier.size = 1, alpha = 0.5)

ggsave("polyA_mj261_mRNAsmallRNA.pdf", p2, height = 10, width = 15, useDingbats=FALSE)


#  subset last 4 clusters

sub <- data[data$label == c('R','S','T','U'),]

sub <- sub[order(sub$smallwtmean),]
head(sub)

c <- c(rep("A",117), rep("B",117), rep("C",117), rep("D",117), rep("E",117), rep("F",117), rep("G",117),rep("H",117),rep("I",117),rep("J",117),rep("K",117),rep("L",117),rep("M",117),rep("N",117),rep("O",117),rep("P",117),rep("R",117),rep("S",117),rep("T",117),rep("U",117))

sub$label <- c

head(sub)

sensor <- subset(sub, sub$genes == "piRNASensor")

p <- ggplot(sub, aes(x=sub$label, y=log10(sub$mRNAwtmean+1))) + geom_boxplot(outlier.colour="Grey", outlier.shape=NA, outlier.size=2.5, notch=FALSE, coef = 6) + geom_boxplot(data=sensor, aes(x=sensor$label, y=log10(sensor$mRNAwtmean+1)), color="purple", size=2) + theme_classic() + scale_y_continuous(limits = c(0,4))+ xlab("smallRNA log10(TPM+1) bins") + ylab("wild type polyA totalRNA log10(TPM+1)")

ggsave("ribo_wildtype_mRNAsmallRNA_last4clusters.pdf", p, height = 10, width = 15, useDingbats=FALSE)

g <- p + geom_boxplot(data=sub, aes(x=sub$label, y=log10(sub$mRNAmj261mean+1)), color="Red",outlier.shape=NA,linetype="dashed", outlier.size = 1, alpha = 0.5) + geom_boxplot(data=sensor, aes(x=sensor$label, y=log10(sensor$mRNAmj261mean+1)), color="green", size=2)

ggsave("ribo_wildtype_mutant_mRNAsmallRNA_last4clusters.pdf", g, height = 10, width = 15, useDingbats=FALSE)

#  subset upregulated genes


polyAsig <- read.csv("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/RNA-seq_ribozero/mj261_ribozero_diff_osc.csv", header = TRUE)

polyAsig2 <- subset(polyAsig, polyAsig$threshold == "A")

datasig <- subset(data, data$genes %in% polyAsig2$X)
datasig <- datasig[order(datasig$smallmj261mean),]
head(datasig)

c <- c(rep("A",44), rep("B",44), rep("C",44), rep("D",44), rep("E",44), rep("F",44), rep("G",44),rep("H",44),rep("I",44),rep("J",44),rep("K",44),rep("L",44),rep("M",44),rep("N",44),rep("O",44),rep("P",44),rep("R",44),rep("S",44),rep("T",44),rep("U",49))

c <- as.vector(c)

datasig$label <- c

head(datasig)

p <- ggplot(datasig, aes(x=datasig$label, y=log10(datasig$mRNAwtmean+1))) + geom_boxplot(outlier.colour="Grey", outlier.shape=NA, outlier.size=2, notch=FALSE, coef = 6) + theme_classic() + scale_y_continuous(limits = c(0,2.5))+ xlab("smallRNA log10(TPM+1) bins") + ylab("wild type polyA totalRNA log10(TPM+1)")

ggsave("ribo_wildtype_mRNAsmallRNA.pdf", p, height = 10, width = 15, useDingbats=FALSE)

g <- p + geom_boxplot(data=datasig, aes(x=datasig$label, y=log10(datasig$mRNAmj261mean+1)), color="Red",outlier.shape=NA,linetype="dashed", outlier.size = 1, alpha = 0.5)

ggsave("ribo_wildtype_mutant_mRNAsmallRNA_upregulated.pdf", g, height = 10, width = 15, useDingbats=FALSE)


data <- data[order(data$smallmj261mean),]
head(data)

p2 <- ggplot(data, aes(x=data$label, y=data$mRNAmj261mean)) + geom_boxplot(outlier.colour="Grey", outlier.shape=NA, outlier.size=2, notch=FALSE, coef = 6) + theme_classic() + scale_y_continuous(limits = c(0,4))+ xlab("smallRNA log10(TPM+1) bins") + ylab("rpb-9(mj261)V polyA totalRNA log10(TPM+1)")


g2 <- p2 +  geom_boxplot(data=data, aes(x=data$label, y=data$mRNAwtmean), color="Red",outlier.shape=NA,linetype="dashed", outlier.size = 1, alpha = 0.5)

ggsave("ribo_mj261_mRNAsmallRNA.pdf", p2, height = 10, width = 15, useDingbats=FALSE)


# CLASS I,II,III genes

classI <- read.delim("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/Figure_4/subsets/ribozero/upregulated_classI_ribo.tsv", header = FALSE)
datasig <- subset(data, data$genes %in% classI$V1)
datasig <- datasig[order(datasig$smallwtmean),]
head(datasig)


c <- c(rep("A",7), rep("B",7), rep("C",7), rep("D",7), rep("E",7), rep("F",7), rep("G",7),rep("H",7),rep("I",7),rep("J",7),rep("K",7),rep("L",7),rep("M",7),rep("N",7),rep("0",7),rep("P",3))


d <- c(rep("A",30), rep("B",30), rep("C",30), rep("D",30), rep("E",29))

d <- as.vector(d)


datasig$label <- c

head(datasig)

p <- ggplot(datasig, aes(x=datasig$label, y=log10(datasig$mRNAwtmean+1))) + geom_boxplot(outlier.colour="Grey", outlier.shape=NA, outlier.size=2, notch=FALSE, coef = 6) + theme_classic() + scale_y_continuous(limits = c(0,2.5))+ xlab("smallRNA log10(TPM+1) bins") + ylab("wild type polyA totalRNA log10(TPM+1)")

ggsave("ribo_wildtype_mRNAsmallRNA_classI.pdf", p, height = 10, width = 15, useDingbats=FALSE)

g <- p + geom_boxplot(data=datasig, aes(x=datasig$label, y=log10(datasig$mRNAmj261mean+1)), color="Red",outlier.shape=NA,linetype="dashed", outlier.size = 1, alpha = 0.5)

ggsave("ribo_wildtype_mutant_mRNAsmallRNA_upregulated.pdf", g, height = 10, width = 15, useDingbats=FALSE)


classII <- read.delim("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/Figure_4/subsets/ribozero/upregulated_classII_ribo.tsv", header = FALSE)
datasig <- subset(data, data$genes %in% classII$V1)
datasig <- datasig[order(datasig$smallwtmean),]
head(datasig)


c <- c(rep("A",7), rep("B",7), rep("C",7), rep("D",7), rep("E",7), rep("F",7), rep("G",7),rep("H",7),rep("I",7),rep("J",7),rep("K",7),rep("L",7),rep("M",7),rep("N",7),rep("0",9))


d <- c(rep("A",30), rep("B",30), rep("C",30), rep("D",30), rep("E",29))

d <- as.vector(d)


datasig$label <- c

head(datasig)

p <- ggplot(datasig, aes(x=datasig$label, y=log10(datasig$mRNAwtmean+1))) + geom_boxplot(outlier.colour="Grey", outlier.shape=NA, outlier.size=2, notch=FALSE, coef = 6) + theme_classic() + scale_y_continuous(limits = c(0,2.5))+ xlab("smallRNA log10(TPM+1) bins") + ylab("wild type polyA totalRNA log10(TPM+1)")

ggsave("polyA_wildtype_mRNAsmallRNA_classII.pdf", p, height = 10, width = 15, useDingbats=FALSE)

g <- p + geom_boxplot(data=datasig, aes(x=datasig$label, y=log10(datasig$mRNAmj261mean+1)), color="Red",outlier.shape=NA,linetype="dashed", outlier.size = 1, alpha = 0.5)

ggsave("polyA_wildtype_mutant_mRNAsmallRNA_upregulated.pdf", g, height = 10, width = 15, useDingbats=FALSE)



classIII <- read.delim("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/Figure_4/subsets/ribozero/upregulated_classIII_ribo.tsv", header = FALSE)
datasig <- subset(data, data$genes %in% classIII$V1)
datasig <- datasig[order(datasig$smallwtmean),]
head(datasig)


c <- c(rep("A",7), rep("B",7), rep("C",7), rep("D",7), rep("E",7), rep("F",7), rep("G",7),rep("H",7),rep("I",7),rep("J",7),rep("K",7),rep("L",7),rep("M",7),rep("N",6))


d <- c(rep("A",30), rep("B",30), rep("C",30), rep("D",30), rep("E",29))

d <- as.vector(d)


datasig$label <- c

head(datasig)

p <- ggplot(datasig, aes(x=datasig$label, y=log10(datasig$mRNAwtmean+1))) + geom_boxplot(outlier.colour="Grey", outlier.shape=NA, outlier.size=2, notch=FALSE, coef = 6) + theme_classic() + scale_y_continuous(limits = c(0,2.5))+ xlab("smallRNA log10(TPM+1) bins") + ylab("wild type polyA totalRNA log10(TPM+1)")

ggsave("polyA_wildtype_mRNAsmallRNA_classIII.pdf", p, height = 10, width = 15, useDingbats=FALSE)

g <- p + geom_boxplot(data=datasig, aes(x=datasig$label, y=log10(datasig$mRNAmj261mean+1)), color="Red",outlier.shape=NA,linetype="dashed", outlier.size = 1, alpha = 0.5)

ggsave("polyA_wildtype_mutant_mRNAsmallRNA_upregulated.pdf", g, height = 10, width = 15, useDingbats=FALSE)

