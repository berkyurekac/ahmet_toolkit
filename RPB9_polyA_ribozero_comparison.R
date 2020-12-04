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
setwd("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/RNA_seq_comparison")


polyA <- read.csv("mj261_polyA_normalized_allreplicates.csv", header = TRUE, stringsAsFactors = FALSE)

ribo <- read.csv("mj261_ribozero_normalized.csv", header = TRUE, stringsAsFactors = FALSE)
colnames(ribo) <- c("X", "wt1x","wt2x","wt3x","mut1x","mut2x","mut3x")

head(polyA)
head(ribo)

merged <- merge(polyA, ribo, by.x="X")
head(merged)

merged$polyAwt <- rowMeans(merged[,2:4])
merged$polyAmut <- rowMeans(merged[,5:7])
merged$ribowt <- rowMeans(merged[,8:10])
merged$ribomut <- rowMeans(merged[,11:13])

head(merged)

p <- ggplot(merged, aes(x=log2(merged$polyAwt+1), y=log2(merged$ribowt+1))) + geom_point() + theme_classic()


genes <- read.delim("mj261_polyA_total_ce11_sub.txt", header = TRUE)
write.csv(genes$Geneid, "protein_coding_genes.csv", row.names = FALSE, quote = FALSE)

mergedsub <- merged[merged$X %in% genes$Geneid,]

p2 <- ggplot(mergedsub, aes(x=log2(mergedsub$polyAwt+1), y=log2(mergedsub$ribowt+1))) + geom_point() + geom_abline(slope =  0.9589, intercept = 2.546e-03 ) + theme_classic() + xlab("polyA_wt") + ylab("Ribo_wt")


reg = lm(mergedsub$polyAwt ~mergedsub$ribowt)
summary(reg) 


ggsave("polyAwt_ribowt.pdf", p2, height = 10, width = 10, useDingbats=FALSE)




p3 <- ggplot(mergedsub, aes(x=log2(mergedsub$polyAmut+1), y=log2(mergedsub$ribomut+1))) + geom_point() + geom_abline(slope =  0.8704, intercept = 5.631e-03 ) + theme_classic() + xlab("polyA_mutant") + ylab("Ribo_mutant")


reg = lm(mergedsub$polyAmut ~mergedsub$ribomut)
summary(reg) 


ggsave("polyAmut_ribomut.pdf", p3, height = 10, width = 10, useDingbats=FALSE)











