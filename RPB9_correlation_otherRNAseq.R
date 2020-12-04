# this script was created by Ahmet Can Berkyurek on 12/12/2018, to correlate differential gene expression analysis with RNA-seq data for pichip candidates.
library(DESeq2)
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
library(data.table)


# load featureCounts count.txt file. before loading, delete the first two rows for clarity.!!!textwrangler should help...
setwd("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/cambridge-UK/NGS/RNAseq/alper/TE")

# first between alleles
# correlation between mj609 and 610
mj609 <- read.csv("lisa_ribozero_allgenes_padj.csv", header = TRUE, stringsAsFactors = FALSE)
rownames(mj609) <- mj609$X
head(mj609)

mj609_sig <- subset(mj609, mj609$threshold != "C")
mj609_sig <- subset(mj609_sig, mj609_sig$sig == "FDR<0.05")
colnames(mj609_sig) <- c("a", "b", "c", "d", "e", "f", "g", "h", "i")
head(mj609_sig)

mj610 <- read.csv("mj613_diff_oscillating_parp1.csv", header = TRUE, stringsAsFactors = FALSE)
rownames(mj610) <- mj610$X
head(mj610)

mj610_sig <- subset(mj610, mj610$sig == "FDR<0.05")
colnames(mj610_sig) <- c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j")
head(mj610_sig)



# introduce hrde-1 and emb-4 data

# combined two data sets

hrde1 <- read.csv("hrde1_diff_osc.csv", header = TRUE, stringsAsFactors = FALSE)
rownames(hrde1) <- hrde1$X
head(hrde1)

hrde1_sig <- subset(hrde1, !(hrde1$threshold) == "C")
head(hrde1_sig)

f <- merge(mj609_sig, hrde1_sig, by.x="a", by.y="X", all = TRUE)
f[is.na(f)] <- 0
f <- data.frame(f)
head(f)

f1 <- data.frame( as.numeric(f$c), as.numeric(f$log2FoldChange))
head(f1)

tnks2corr <- cor(f1$as.numeric.f.c., f1$as.numeric.f.log2FoldChange., method = "spearman", use = "pairwise.complete.obs")
tnks2corr
# -0.02601978

breaks <- c(-10,-6,-4,-2,-1,0,1,2,4,6,10)

g <- ggplot(f1, aes(x=f1$as.numeric.f.c., y=f1$as.numeric.f.log2FoldChange.)) + geom_point(color = "Grey", size=3) + scale_x_continuous(limits = c(-10,10)) + scale_y_continuous(limits = c(-10,10)) + theme_classic() + xlab("log2(tnks-2(mj609))") + ylab("log2(hrde-1(tm1200))")

ggsave("mj261_SX2000_correlation_new.pdf", g, height = 10, width = 10, useDingbats=FALSE)


#subset emb-4 targets
emb4 <- read.csv("emb4_diff_osc.csv", header = TRUE, stringsAsFactors = FALSE)
rownames(emb4) <- emb4$X
head(emb4)

emb4_sig <- subset(emb4, !(emb4$threshold) == "C")
head(emb4_sig)


f <- merge(mj609_sig, emb4_sig, by.x="a", by.y="X", all = TRUE)
f[is.na(f)] <- 0
f <- as.data.frame(f)
head(f)

f1 <- data.frame( as.numeric(f$c), as.numeric(f$log2FoldChange))
head(f1)

tnks2corr <- cor(f1$as.numeric.f.log2FoldChange., f1$as.numeric.f.c., method = "spearman", use = "pairwise.complete.obs")
tnks2corr
# 0.4048734

breaks <- c(-10,-6,-4,-2,-1,0,1,2,4,6,10)

g <- ggplot(f1, aes(x=f1$as.numeric.f.c., y=f1$as.numeric.f.log2FoldChange.)) + geom_point(color = "Grey", size=3) + scale_x_continuous(limits = c(-10,10)) + scale_y_continuous(limits = c(-10,10)) + theme_classic() + xlab("log2(tnks-2(mj609))") + ylab("log2(emb-4)")

ggsave("mj261_SX2929_correlation_new.pdf", g, height = 10, width = 10, useDingbats=FALSE)

#subset set-25 and set-32 targets
set25 <- read.csv("set32_vs_wt_ribozero_allgenes.csv", header = TRUE, stringsAsFactors = FALSE)
rownames(set25) <- set25$X
head(set25)

set25_sig <- subset(set25, abs(set25$log2FoldChange) >= 1)
set25_sig <- subset(set25_sig, set25_sig$padj <= 0.05)
head(set25_sig)


f <- merge(mj609_sig, set25_sig, by.x="a", by.y="X", all = TRUE)
f[is.na(f)] <- 0
f <- as.data.frame(f)
head(f)

f1 <- data.frame( as.numeric(f$c), as.numeric(f$log2FoldChange))
head(f1)

tnks2corr <- cor(f1$as.numeric.f.log2FoldChange., f1$as.numeric.f.c., method = "spearman", use = "pairwise.complete.obs")
tnks2corr
# 0.4048734

breaks <- c(-10,-6,-4,-2,-1,0,1,2,4,6,10)

g <- ggplot(f1, aes(x=f1$as.numeric.f.c., y=f1$as.numeric.f.log2FoldChange.)) + geom_point(color = "Grey", size=3) + scale_x_continuous(limits = c(-10,10)) + scale_y_continuous(limits = c(-10,10)) + theme_classic() + xlab("log2(tnks-2(mj609))") + ylab("log2(emb-4)")

ggsave("mj261_SX2929_correlation_new.pdf", g, height = 10, width = 10, useDingbats=FALSE)
