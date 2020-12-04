# this script was created by Ahmet Can Berkyurek on 15/12/2019, to check for correlation between RPB-9 and PolII ChIP-seq data.
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


setwd("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/cambridge-UK/NGS/ChIP-seq/ollas-RPB9/RPB9_vs_PolII/new")

rpb_rep1 <- read.table("SX3320_vs_N2_log2_ChIP_log2_average.txt", header = FALSE, stringsAsFactors = FALSE)
rpb_rep2 <- read.table("SX3321_vs_N2_log2_ChIP_log2_average.txt", header = FALSE, stringsAsFactors = FALSE)

ama_rep1 <- read.table("SX1316_rep1_log2_ChIP_log2_average.txt", header = FALSE, stringsAsFactors = FALSE)
ama_rep2 <- read.table("SX1316_rep2_log2_ChIP_log2_average.txt", header = FALSE, stringsAsFactors = FALSE)
ama_rep3 <- read.table("SX1316_rep3_log2_ChIP_log2_average.txt", header = FALSE, stringsAsFactors = FALSE)


head(rpb_rep1)

frame <- data.frame(rpb_rep1[,c(1,6)], rpb_rep2$V6, ama_rep1$V6, ama_rep2$V6, ama_rep3$V6)
colnames(frame) <- c("genes", "rpb_rep1", "rpb_rep2", "ama_rep1", "ama_rep2", "ama_rep3")
head(frame)

frame$rpb_mean <- rowMeans(frame[,2:3])
frame$ama_mean <- rowMeans(frame[,4:6])

head(frame)

p <- ggplot(frame, aes(x=frame$rpb_mean, y=frame$ama_mean))  + geom_point(color= "Grey")+ scale_y_continuous(limits = c(-6,6)) + scale_x_continuous(limits = c(-6,6)) + theme_classic() + geom_vline(xintercept = 0, colour = "grey")

ggsave("RPB9_vs_polII.pdf", p, height = 10, width = 10, useDingbats=FALSE)



frame_rpb_bigger_ama <- subset(frame, frame$rpb_mean > 0 & frame$ama_mean < 0)

gonad <- read.delim("oogenesis_genes_reinke.bed", header = FALSE)
sperm <- read.delim("sperm_genes_reinke.bed", header = FALSE)
neuron <- read.delim("neuron_genes_watson.bed", header = FALSE)
Ubiq <- read.delim("Ubiq_genes_names.bed", header = TRUE)

frameubiq <- subset(frame, frame$genes %in% Ubiq$WBGene00002061)

p <- ggplot(frameubiq, aes(x=frameubiq$rpb_mean, y=frameubiq$ama_mean))  + geom_point(color= "Grey") + scale_y_continuous(limits = c(-6,6)) + scale_x_continuous(limits = c(-6,6)) + theme_classic() + geom_vline(xintercept = 0, colour = "grey") 

ggsave("RPB9_vs_polIIUbiqgenes.pdf", p, height = 10, width = 10, useDingbats=FALSE)



framegonad <- subset(frame, frame$genes %in% gonad$V1)

p <- ggplot(framegonad, aes(x=framegonad$rpb_mean, y=framegonad$ama_mean))  + geom_point(color= "Grey") + scale_y_continuous(limits = c(-6,6)) + scale_x_continuous(limits = c(-6,6)) + theme_classic() + geom_vline(xintercept = 0, colour = "grey") 

ggsave("RPB9_vs_polIIgonadgenes.pdf", p, height = 10, width = 10, useDingbats=FALSE)


framesperm <- subset(frame, frame$genes %in% sperm$V1)

p <- ggplot(framesperm, aes(x=framesperm$rpb_mean, y=framesperm$ama_mean))  + geom_point(color= "Grey") + scale_y_continuous(limits = c(-6,6)) + scale_x_continuous(limits = c(-6,6)) + theme_classic() + geom_vline(xintercept = 0, colour = "grey") 

ggsave("RPB9_vs_polIIspermgenes.pdf", p, height = 10, width = 10, useDingbats=FALSE)



frameneuron <- subset(frame, frame$genes %in% neuron$V1)

p <- ggplot(frameneuron, aes(x=frameneuron$rpb_mean, y=frameneuron$ama_mean))  + geom_point(color= "Grey") + scale_y_continuous(limits = c(-6,6)) + scale_x_continuous(limits = c(-6,6)) + theme_classic() + geom_vline(xintercept = 0, colour = "grey") 

ggsave("RPB9_vs_polIIneurongenes.pdf", p, height = 10, width = 10, useDingbats=FALSE)


# 




p <- ggplot(frame, aes(x=frame$rpb_rep1, y=frame$rpb_rep2))  + geom_point(color= "Grey") + scale_y_continuous(limits = c(-6,6)) + scale_x_continuous(limits = c(-6,6)) + theme_classic()

ggsave("RPB9_replicates.pdf", p, height = 10, width = 10, useDingbats=FALSE)

p <- ggplot(frame, aes(x=frame$ama_rep1, y=frame$ama_rep2))  + geom_point(color= "Grey") + scale_y_continuous(limits = c(-6,6)) + scale_x_continuous(limits = c(-6,6)) + theme_classic()

ggsave("AMA1_replicates_1.pdf", p, height = 10, width = 10, useDingbats=FALSE)

p <- ggplot(frame, aes(x=frame$ama_rep2, y=frame$ama_rep3))  + geom_point(color= "Grey") + scale_y_continuous(limits = c(-6,6)) + scale_x_continuous(limits = c(-6,6)) + theme_classic()

ggsave("AMA1_replicates_2.pdf", p, height = 10, width = 10, useDingbats=FALSE)

















