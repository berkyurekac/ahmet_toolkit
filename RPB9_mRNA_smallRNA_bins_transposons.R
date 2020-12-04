# this script was created by Ahmet Can Berkyurek on 14/04/2020, to correlate mRNA and smallRNa expression in bins for transposons. 
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

setwd("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/RNA-seq_polyA/transposons/all_replicates/")

#introduce the countdata
wt1x <- read.delim("sx1316_rep1_combinedAligned.sortedByCoord_mj261polyA_TE_unique.txt", header = FALSE, stringsAsFactors = FALSE)

wt2x <- read.delim("sx1316_rep2_combined_polyA_mj261Aligned.sortedByCoord_mj261polyA_TE_unique.txt", header = FALSE, stringsAsFactors = FALSE)

wt3x <- read.delim("sx1316_rep3_combined_polyA_mj261Aligned.sortedByCoord_mj261polyA_TE_unique.txt", header = FALSE, stringsAsFactors = FALSE)

mj261_1x <- read.delim("sx1902_rep1_combined_polyA_mj261Aligned.sortedByCoord_mj261polyA_TE_unique.txt", header = FALSE, stringsAsFactors = FALSE)

mj261_2x <- read.delim("sx1902_rep2_combined_polyA_mj261Aligned.sortedByCoord_mj261polyA_TE_unique.txt", header = FALSE, stringsAsFactors = FALSE)

mj261_3x <- read.delim("sx1902_rep3_combined_polyA_mj261Aligned.sortedByCoord_mj261polyA_TE_unique.txt", header = FALSE, stringsAsFactors = FALSE)


countdata <- data.frame(wt1x[,c(1:4,7)], wt2x[,7], wt3x[,7], mj261_1x[,7], mj261_2x[,7], mj261_3x[,7])

colnames(countdata) <- c("chr", "start", "end","genes","wt1x", "wt2x", "wt3x", "mj261_1x", "mj261_2x", "mj261_3x")
head(countdata)

countspolyA <- countdata

countspolyA$length <- countspolyA$end-countspolyA$start
head(countspolyA)

countspolyA <- subset(countspolyA, rowMeans(countspolyA[,5:10]) > 1)




#TPM normalization Divide the read counts by the length of each gene in kilobases. This gives you reads per kilobase (RPK).
#Count up all the RPK values in a sample and divide this number by 1,000,000. This is your “per million” scaling factor.
#Divide the RPK values by the “per million” scaling factor. This gives you TPM.

countspolyA[,5:10] <- countspolyA[,5:10]/((countspolyA$length)*0.001)

head(countspolyA)

countspolyA$wt1x <- countspolyA$wt1x/(sum(countspolyA$wt1x)/1e6)
countspolyA$wt2x <- countspolyA$wt2x/(sum(countspolyA$wt2x)/1e6)
countspolyA$wt3x <- countspolyA$wt3x/(sum(countspolyA$wt3x)/1e6)

countspolyA$mj261_1x <- countspolyA$mj261_1x/(sum(countspolyA$mj261_1x)/1e6)
countspolyA$mj261_2x <- countspolyA$mj261_2x/(sum(countspolyA$mj261_2x)/1e6)
countspolyA$mj261_3x <- countspolyA$mj261_3x/(sum(countspolyA$mj261_3x)/1e6)

head(countspolyA)


setwd("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/RNA-seq_polyA/transposons")

#smallRNA libraries 
wt1 <- read.delim("isa_1_mj261__mj261smallRNA_TE_unique.txt", header = FALSE, stringsAsFactors = FALSE)

wt2 <- read.delim("isa_3_mj261__mj261smallRNA_TE_unique.txt", header = FALSE, stringsAsFactors = FALSE)

wt3 <- read.delim("isa_3_mj261__mj261smallRNA_TE_unique.txt", header = FALSE, stringsAsFactors = FALSE)

mj261_1 <- read.delim("isa_4_mj261__mj261smallRNA_TE_unique.txt", header = FALSE, stringsAsFactors = FALSE)

mj261_2 <- read.delim("isa_5_mj261__mj261smallRNA_TE_unique.txt", header = FALSE, stringsAsFactors = FALSE)

mj261_3 <- read.delim("isa_6_mj261__mj261smallRNA_TE_unique.txt", header = FALSE, stringsAsFactors = FALSE)

rescue1 <- read.delim("isa_7_mj261__mj261smallRNA_TE_unique.txt", header = FALSE, stringsAsFactors = FALSE)

rescue2 <- read.delim("isa_8_mj261__mj261smallRNA_TE_unique.txt", header = FALSE, stringsAsFactors = FALSE)

rescue3 <- read.delim("isa_9_mj261__mj261smallRNA_TE_unique.txt", header = FALSE, stringsAsFactors = FALSE)


countdata <- data.frame(wt1[,c(1:4,7)], wt2[,7], wt3[,7], mj261_1[,7], mj261_2[,7], mj261_3[,7],rescue1[,7],rescue2[,7],rescue3[,7])

colnames(countdata) <- c("chr", "start", "end","genes","wt1", "wt2", "wt3", "mj261_1", "mj261_2", "mj261_3", "rescue1", "rescue2", "rescue3")
head(countdata)

countssmall <- countdata

countssmall$length <- countssmall$end-countssmall$start
head(countssmall)




#TPM normalization Divide the read counts by the length of each gene in kilobases. This gives you reads per kilobase (RPK).
#Count up all the RPK values in a sample and divide this number by 1,000,000. This is your “per million” scaling factor.
#Divide the RPK values by the “per million” scaling factor. This gives you TPM.

countssmall[,5:13] <- countssmall[,5:13]/((countssmall$length)*0.001)

head(countssmall)

countssmall$wt1 <- countssmall$wt1/(sum(countssmall$wt1)/1e6)
countssmall$wt2 <- countssmall$wt2/(sum(countssmall$wt2)/1e6)
countssmall$wt3 <- countssmall$wt3/(sum(countssmall$wt3)/1e6)

countssmall$mj261_1 <- countssmall$mj261_1/(sum(countssmall$mj261_1)/1e6)
countssmall$mj261_2 <- countssmall$mj261_2/(sum(countssmall$mj261_2)/1e6)
countssmall$mj261_3 <- countssmall$mj261_3/(sum(countssmall$mj261_3)/1e6)

countssmall$rescue1 <- countssmall$rescue1/(sum(countssmall$rescue1)/1e6)
countssmall$rescue2 <- countssmall$rescue2/(sum(countssmall$rescue2)/1e6)
countssmall$rescue3 <- countssmall$rescue3/(sum(countssmall$rescue3)/1e6)

head(countssmall)

countssmall <- subset(countssmall, countssmall$genes %in% countspolyA$genes)


# combine polyA mRNA and small RNA - first cluster with small RNA small to high, then label rows for bining!

polyAmerged <- data.frame(countspolyA[1:10],countssmall[5:13])
head(polyAmerged)

polyAmerged$mRNAwtmean <- rowMeans(polyAmerged[,5:7])
polyAmerged$mRNAmj261mean <- rowMeans(polyAmerged[,8:10])

polyAmerged$smallwtmean <- rowMeans(polyAmerged[,11:13])
polyAmerged$smallmj261mean <- rowMeans(polyAmerged[,14:16])
polyAmerged$smallrescuemean <- rowMeans(polyAmerged[,17:19])


head(polyAmerged)

data <- polyAmerged[,c(4:24)]
head(data)
# order data in smallRNA mutant increasing order
data <- data[order(data$smallwtmean),]
head(data)


b <- c(rep("A",764), rep("B",764), rep("C",764), rep("D",764), rep("E",764), rep("F",764), rep("G",764),rep("H",764),rep("I",764),rep("J",768))

b <- as.vector(b)

data$label <- b

head(data)

write.csv(data, "RPB9polyA_transposons_smallRNA_bins_log10TPM.csv")

CEMUDR <- subset(data, data$genes == "CEMUDR1;2902771;2909933")
CHAPAEV  <- subset(data, data$genes == "Chapaev-2_CE;2910168;2913994")

sig <- rbind(CEMUDR,CHAPAEV)

p <- ggplot(data, aes(x=data$label, y=log2(data$mRNAwtmean))) + geom_boxplot(outlier.colour="Grey", outlier.shape=NA, outlier.size=2.5, notch=FALSE, coef = 6) + geom_boxplot(data=sig, aes(x=sig$label, y=log2(sig$mRNAwtmean)), color="Black", size=2) + theme_classic() + scale_y_continuous(limits = c(0,1))+ xlab("smallRNA log10(TPM) bins") + ylab("wild type polyA totalRNA log10(TPM)")

ggsave("polyA_wildtype_mRNAsmallRNA.pdf", p, height = 10, width = 15, useDingbats=FALSE)

g <- p + geom_boxplot(data=data, position = "dodge",aes(x=data$label, y=log2(data$mRNAmj261mean)), color="Red",outlier.shape=NA,linetype="dashed", outlier.size = 1, alpha = 0.5) + geom_boxplot(data=sig, aes(x=sig$label, y=log2(sig$mRNAmj261mean)), color="Blue", size=2)

ggsave("polyA_wildtype_mutant_mRNAsmallRNA_transposons.pdf", g, height = 10, width = 15, useDingbats=FALSE)


data <- data[order(data$smallmj261mean),]
head(data)

p2 <- ggplot(data, aes(x=data$label, y=log10(data$mRNAwtmean))) + geom_boxplot(outlier.colour="Grey", outlier.shape=NA, outlier.size=2, notch=FALSE, coef = 6) + geom_boxplot(data=sig, aes(x=sig$label, y=log10(sig$mRNAwtmean)), color="Grey", size=2) + theme_classic() + scale_y_continuous(limits = c(0,4))+ xlab("smallRNA log10(TPM+1) bins") + ylab("rpb-9(mj261)V polyA totalRNA log10(TPM+1)")


g2 <- p2 +  geom_boxplot(data=data, aes(x=data$label, y=log10(data$mRNAmj261mean)), color="Red",outlier.shape=NA,linetype="dashed", outlier.size = 1, alpha = 0.5) + geom_boxplot(data=sig, aes(x=sig$label, y=log10(sig$mRNAmj261mean)), color="Blue", size=2)

ggsave("polyA_mj261_mRNAsmallRNA_transposons.pdf", g2, height = 10, width = 15, useDingbats=FALSE)





#######
# ribo-zero total RNA-seq
######
setwd("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/RNA-seq_ribozero/transposons/")

#introduce the countdata
wt1x <- read.delim("SLX-15940.NEBNext01.HWFK7BBXX.s_5.r_1_ribo_mj261Aligned.sortedByCoord.out_mj261ribo_TE_unique.txt", header = FALSE, stringsAsFactors = FALSE)

wt2x <- read.delim("SLX-15940.NEBNext02.HWFK7BBXX.s_5.r_1.ribo_mj261Aligned.sortedByCoord.out_mj261ribo_TE_unique.txt", header = FALSE, stringsAsFactors = FALSE)

wt3x <- read.delim("SLX-15940.NEBNext03.HWFK7BBXX.s_5.r_1_ribo_mj261Aligned.sortedByCoord.out_mj261ribo_TE_unique.txt", header = FALSE, stringsAsFactors = FALSE)

mj261_1x <- read.delim("SLX-15940.NEBNext10.HWFK7BBXX.s_5.r_1_ribo_mj261Aligned.sortedByCoord.out_mj261ribo_TE_unique.txt", header = FALSE, stringsAsFactors = FALSE)

mj261_2x <- read.delim("SLX-15940.NEBNext11.HWFK7BBXX.s_5.r_1_ribo_mj261Aligned.sortedByCoord.out_mj261ribo_TE_unique.txt", header = FALSE, stringsAsFactors = FALSE)

mj261_3x <- read.delim("SLX-15940.NEBNext12.HWFK7BBXX.s_5.r_1_ribo_mj261Aligned.sortedByCoord.out_mj261ribo_TE_unique.txt", header = FALSE, stringsAsFactors = FALSE)


countdata <- data.frame(wt1x[,c(1:4,7)], wt2x[,7], wt3x[,7], mj261_1x[,7], mj261_2x[,7], mj261_3x[,7])

colnames(countdata) <- c("chr", "start", "end","genes","wt1x", "wt2x", "wt3x", "mj261_1x", "mj261_2x", "mj261_3x")
head(countdata)

countspolyA <- countdata

countspolyA$length <- countspolyA$end-countspolyA$start
head(countspolyA)





#TPM normalization Divide the read counts by the length of each gene in kilobases. This gives you reads per kilobase (RPK).
#Count up all the RPK values in a sample and divide this number by 1,000,000. This is your “per million” scaling factor.
#Divide the RPK values by the “per million” scaling factor. This gives you TPM.

countspolyA[,5:10] <- countspolyA[,5:10]/((countspolyA$length)*0.001)

head(countspolyA)

countspolyA$wt1x <- countspolyA$wt1x/(sum(countspolyA$wt1x)/1e6)
countspolyA$wt2x <- countspolyA$wt2x/(sum(countspolyA$wt2x)/1e6)
countspolyA$wt3x <- countspolyA$wt3x/(sum(countspolyA$wt3x)/1e6)

countspolyA$mj261_1x <- countspolyA$mj261_1x/(sum(countspolyA$mj261_1x)/1e6)
countspolyA$mj261_2x <- countspolyA$mj261_2x/(sum(countspolyA$mj261_2x)/1e6)
countspolyA$mj261_3x <- countspolyA$mj261_3x/(sum(countspolyA$mj261_3x)/1e6)

countspolyA <- subset(countspolyA, rowMeans(countspolyA[,5:10]) > 0.5)

head(countspolyA)


setwd("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/RNA-seq_polyA/transposons")

#smallRNA libraries 
wt1 <- read.delim("isa_1_mj261__mj261smallRNA_TE_unique.txt", header = FALSE, stringsAsFactors = FALSE)

wt2 <- read.delim("isa_3_mj261__mj261smallRNA_TE_unique.txt", header = FALSE, stringsAsFactors = FALSE)

wt3 <- read.delim("isa_3_mj261__mj261smallRNA_TE_unique.txt", header = FALSE, stringsAsFactors = FALSE)

mj261_1 <- read.delim("isa_4_mj261__mj261smallRNA_TE_unique.txt", header = FALSE, stringsAsFactors = FALSE)

mj261_2 <- read.delim("isa_5_mj261__mj261smallRNA_TE_unique.txt", header = FALSE, stringsAsFactors = FALSE)

mj261_3 <- read.delim("isa_6_mj261__mj261smallRNA_TE_unique.txt", header = FALSE, stringsAsFactors = FALSE)

rescue1 <- read.delim("isa_7_mj261__mj261smallRNA_TE_unique.txt", header = FALSE, stringsAsFactors = FALSE)

rescue2 <- read.delim("isa_8_mj261__mj261smallRNA_TE_unique.txt", header = FALSE, stringsAsFactors = FALSE)

rescue3 <- read.delim("isa_9_mj261__mj261smallRNA_TE_unique.txt", header = FALSE, stringsAsFactors = FALSE)


countdata <- data.frame(wt1[,c(1:4,7)], wt2[,7], wt3[,7], mj261_1[,7], mj261_2[,7], mj261_3[,7],rescue1[,7],rescue2[,7],rescue3[,7])

colnames(countdata) <- c("chr", "start", "end","genes","wt1", "wt2", "wt3", "mj261_1", "mj261_2", "mj261_3", "rescue1", "rescue2", "rescue3")
head(countdata)

countssmall <- countdata

countssmall$length <- countssmall$end-countssmall$start
head(countssmall)




#TPM normalization Divide the read counts by the length of each gene in kilobases. This gives you reads per kilobase (RPK).
#Count up all the RPK values in a sample and divide this number by 1,000,000. This is your “per million” scaling factor.
#Divide the RPK values by the “per million” scaling factor. This gives you TPM.

countssmall[,5:13] <- countssmall[,5:13]/((countssmall$length)*0.001)

head(countssmall)

countssmall$wt1 <- countssmall$wt1/(sum(countssmall$wt1)/1e6)
countssmall$wt2 <- countssmall$wt2/(sum(countssmall$wt2)/1e6)
countssmall$wt3 <- countssmall$wt3/(sum(countssmall$wt3)/1e6)

countssmall$mj261_1 <- countssmall$mj261_1/(sum(countssmall$mj261_1)/1e6)
countssmall$mj261_2 <- countssmall$mj261_2/(sum(countssmall$mj261_2)/1e6)
countssmall$mj261_3 <- countssmall$mj261_3/(sum(countssmall$mj261_3)/1e6)

countssmall$rescue1 <- countssmall$rescue1/(sum(countssmall$rescue1)/1e6)
countssmall$rescue2 <- countssmall$rescue2/(sum(countssmall$rescue2)/1e6)
countssmall$rescue3 <- countssmall$rescue3/(sum(countssmall$rescue3)/1e6)

head(countssmall)

countssmall <- subset(countssmall, countssmall$genes %in% countspolyA$genes)


# combine polyA mRNA and small RNA - first cluster with small RNA small to high, then label rows for bining!

polyAmerged <- data.frame(countspolyA[1:10],countssmall[5:13])
head(polyAmerged)

polyAmerged$mRNAwtmean <- rowMeans(polyAmerged[,5:7])
polyAmerged$mRNAmj261mean <- rowMeans(polyAmerged[,8:10])

polyAmerged$smallwtmean <- rowMeans(polyAmerged[,11:13])
polyAmerged$smallmj261mean <- rowMeans(polyAmerged[,14:16])
polyAmerged$smallrescuemean <- rowMeans(polyAmerged[,17:19])


head(polyAmerged)

data <- polyAmerged[,c(4:24)]
head(data)
# order data in smallRNA mutant increasing order
data <- data[order(data$smallwtmean),]
head(data)

b <- c(rep("A",3116), rep("B",3116), rep("C",3116), rep("D",3116), rep("E",3116), rep("F",3116), rep("G",3116),rep("H",3116),rep("I",3116),rep("J",3116),rep("K",3116),rep("L",3116),rep("M",3116),rep("N",3116),rep("O",3116),rep("P",3116),rep("R",3116),rep("S",3116),rep("T",3116),rep("U",3127))

b <- as.vector(b)

data$label <- b




b <- c(rep("A",1677), rep("B",1677), rep("C",1677), rep("D",1677), rep("E",1677), rep("F",1677), rep("G",1677),rep("H",1677),rep("I",1677),rep("J",1678))

b <- as.vector(b)

data$label <- b

head(data)

write.csv(data, "RPB9ribo_transposons_smallRNA_bins_log10TPM.csv")

CEMUDR <- subset(data, data$genes == "CEMUDR1;2902771;2909933")
CHAPAEV  <- subset(data, data$genes == "Chapaev-2_CE;2910168;2913994")

sig <- rbind(CEMUDR,CHAPAEV)

p <- ggplot(data, aes(x=data$label, y=log10(data$mRNAwtmean))) + geom_boxplot(outlier.colour="Grey", outlier.shape=NA, outlier.size=2.5, notch=FALSE, coef = 6) + geom_boxplot(data=sig, aes(x=sig$label, y=log10(sig$mRNAwtmean)), color="Black", size=1) + theme_classic() + scale_y_continuous(limits = c(0,1.5))+ xlab("smallRNA log10(TPM) bins") + ylab("wild type ribozero totalRNA log10(TPM)")

ggsave("polyA_wildtype_mRNAsmallRNA.pdf", p, height = 10, width = 15, useDingbats=FALSE)

g <- p + geom_boxplot(data=data, position = "dodge",aes(x=data$label, y=log10(data$mRNAmj261mean)), color="Red",outlier.shape=NA,linetype="dashed", outlier.size = 1, alpha = 0.5) + geom_boxplot(data=sig, aes(x=sig$label, y=log10(sig$mRNAmj261mean)), color="Blue", size=2)

ggsave("polyA_wildtype_mutant_mRNAsmallRNA_transposons.pdf", g, height = 10, width = 15, useDingbats=FALSE)

