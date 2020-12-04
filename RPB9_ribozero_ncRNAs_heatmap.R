# this script was created by Ahmet Can Berkyurek on 21/04/2020, to analyze ncRNAs in the ribozer oRPB9 RNA-seq data
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
setwd("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/RNA-seq_ribozero/snRNAs/")


########
# introduce ncRNA genes
########


down <- read.delim("snRNAs_diff.txt", header = TRUE)

head(down)

downNC <- subset(down, down$Gene.type =="ncRNA")
downSN <- subset(down, down$Gene.type =="snRNA")
downSNo <- subset(down, down$Gene.type =="snoRNA")

snRNA_all <- read.delim("snRNAs.txt", header = TRUE)

ncRNA_all <- read.delim("ncRNAs.txt", header = T)

snoRNA_all <- read.delim("snoRNAs.txt", header = T)

ribo <- read.csv("mj261_ribozero_normalized.csv", header = TRUE, stringsAsFactors = FALSE)

head(ribo)

ribo$meanwt <- rowMeans(ribo[,2:4])
head(ribo)

ribo[,2:8] <- log2((ribo[,2:8]+1)/(ribo$meanwt+1))
head(ribo)

riboNC <- subset(ribo, ribo$X %in% downNC$Gene.stable.ID)
head(riboNC)
write.csv(riboNC, "downregulated_ncRNAs_ribozero.csv")

riboSN <- subset(ribo, ribo$X %in% downSN$Gene.stable.ID)
write.csv(riboSN, "downregulated_snRNAs_ribozero.csv")

ribosnRNAall <- subset(ribo, ribo$X %in% snRNA_all$Gene.stable.ID)
write.csv(ribosnRNAall, "downregulated_snRNAs_all_ribozero.csv")

riboncRNAall <- subset(ribo, ribo$X %in% ncRNA_all$Gene.stable.ID)
write.csv(ribosnRNAall, "downregulated_ncRNAs_all_ribozero.csv")

ribosnoRNAall <- subset(ribo, ribo$X %in% snoRNA_all$Gene.stable.ID)
write.csv(ribosnRNAall, "downregulated_snoRNAs_all_ribozero.csv")

riboSNo <- subset(ribo, ribo$X %in% downSNo$Gene.stable.ID)
write.csv(riboSNo, "downregulated_snoRNAs_all_ribozero.csv")


# draw heatmaps fpr each subset!
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 299)
col_breaks = c(seq(-5,-3,length=100),  # for red
               seq(-2.99,2.99,length=100),           # for white
               seq(3,5,length=100))

x <- data.frame(col_breaks)
image(rotate(as.matrix(x)), col = my_palette, breaks = col_breaks)


# all genes
rotate <- function(x) t(apply(x, 2, rev))
image(rotate(as.matrix(riboNC[2:7])), col = my_palette, breaks = col_breaks)

rotate <- function(x) t(apply(x, 2, rev))
image(rotate(as.matrix(riboSN[2:7])), col = my_palette, breaks = col_breaks)

rotate <- function(x) t(apply(x, 2, rev))
image(rotate(as.matrix(ribosnRNAall[2:7])), col = my_palette, breaks = col_breaks)

rotate <- function(x) t(apply(x, 2, rev))
image(rotate(as.matrix(riboncRNAall[2:7])), col = my_palette, breaks = col_breaks)

rotate <- function(x) t(apply(x, 2, rev))
image(rotate(as.matrix(ribosnoRNAall[2:7])), col = my_palette, breaks = col_breaks)



rotate <- function(x) t(apply(x, 2, rev))
image(rotate(as.matrix(riboSNo[2:7])), col = my_palette, breaks = col_breaks)

setwd("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/RNA-seq_ribozero/histones/")

histones <- read.delim("histone_genes.txt", header = TRUE)
ribohistone <- subset(ribo, ribo$X %in% histones$Gene.stable.ID)

ribo$meanwt <- rowMeans(ribo[,2:4])
head(ribo)
ribo[,2:8] <- (ribo[,2:8]+1)/(ribo$meanwt+1)
head(ribo)

ribohistone2 <- subset(ribo, ribo$X %in% histones$Gene.stable.ID)


rotate <- function(x) t(apply(x, 2, rev))
image(rotate(as.matrix(ribohistone[2:7])), col = my_palette, breaks = col_breaks)

setwd("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/RNA-seq_ribozero/")

sig <- read.csv("mj261_ribozero_diff_osc.csv", header = T, stringsAsFactors = F)
head(sig)

sig <- as.data.frame(sig)

sig$log2FoldChange <- as.numeric(sig$log2FoldChange)

histones <- read.delim("./histones/histone_genes.txt", header = T)
head(histones)

sighistones <- sig[sig$X %in% histones$Gene.stable.ID,]
write.csv(sighistones, "histone_diff_ribozero.csv")

up <- sig[sig$threshold == "A",]
up <- na.omit(up)
down <- sig[sig$threshold == "B",]
down <- na.omit(down)
Sensor <- sig[sig$genes == "piRNASensor",]


# volcano plot
p <- ggplot(sig, aes(x=log2FoldChange, y=-log10(padj))) +
  geom_point(color = "Grey", size=3) + geom_point(data=up, aes(x=up$log2FoldChange, y=-log10(up$padj)), color = "red", size=4.5) + geom_point(data=Sensor, aes(x=Sensor$log2FoldChange, y=-log10(Sensor$padj)), color = "darkgreen", size=6)  + geom_point(data=down, aes(x=down$log2FoldChange, y=-log10(down$padj)), color = "blue", size=4.5) +
  scale_x_continuous(limits = c(-10,10)) + geom_point(data=sighistones, aes(x=sighistones$log2FoldChange, y=-log10(sighistones$padj)), size=2, color="black") +
  geom_vline(xintercept = 0, colour = "grey") + 
  theme_classic() + xlab("log2FoldChange") + ylab("-log10(padj)")

ggsave("./histones/ribozero_histoneslabeled_volcano.pdf", p, height = 15, width = 10, useDingbats=FALSE)




