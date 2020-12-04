# this script was created by Ahmet Can Berkyurek on 01/04/2020, to Transposon analyze differential gene expression analysis for mj261 rpb-9 mutant
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

#polyA libraries
# load featureCounts count.txt file. before loading, delete the first two rows for clarity.!!!textwrangler should help...
# filter too many reads from ribo and polyA libraries. file name is "ribo_highreads.csv"
setwd("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/RNA-seq_polyA/transposons/all_replicates/")

#introduce the countdata from featurecounts (fractional)
wt1 <- read.delim("sx1316_rep1_combinedAligned.sortedByCoord_mj261polyA_TE_unique.txt", header = FALSE, stringsAsFactors = FALSE)

wt2 <- read.delim("sx1316_rep2_combined_polyA_mj261Aligned.sortedByCoord_mj261polyA_TE_unique.txt", header = FALSE, stringsAsFactors = FALSE)

wt3 <- read.delim("sx1316_rep3_combined_polyA_mj261Aligned.sortedByCoord_mj261polyA_TE_unique.txt", header = FALSE, stringsAsFactors = FALSE)

mj261_1 <- read.delim("sx1902_rep1_combined_polyA_mj261Aligned.sortedByCoord_mj261polyA_TE_unique.txt", header = FALSE, stringsAsFactors = FALSE)

mj261_2 <- read.delim("sx1902_rep2_combined_polyA_mj261Aligned.sortedByCoord_mj261polyA_TE_unique.txt", header = FALSE, stringsAsFactors = FALSE)

mj261_3 <- read.delim("sx1902_rep3_combined_polyA_mj261Aligned.sortedByCoord_mj261polyA_TE_unique.txt", header = FALSE, stringsAsFactors = FALSE)


countdata <- data.frame(wt1[,c(4,7)], wt2[,7], wt3[,7], mj261_1[,7], mj261_2[,7], mj261_3[,7])

colnames(countdata) <- c("genes", "wt1", "wt2", "wt3", "mj261_1", "mj261_2", "mj261_3")
head(countdata)

countdata_mj609 <- countdata[,2:7]
rownames(countdata_mj609) <- countdata$genes
head(countdata_mj609)


plotPCAEx = function(object, PCx = 1, PCy = 2, cond="condition", ntop=500, labels = TRUE)
{
  library(RColorBrewer)
  # calculate the variance for each gene
  rv <- rowVars(assay(object))
  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(object)[select,]))
  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  if (!all(cond %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  cond.df <- as.data.frame(colData(object)[, cond, drop=FALSE])
  # add the intgroup factors together to create a new grouping factor
  group <- if (length(cond) > 1) {
    factor(apply( cond.df, 1, paste, collapse=" : "))
  } else {
    colData(object)[[cond]]
  }
  # assembly the data for the plot
  d <- data.frame(PCa=pca$x[,PCx], PCb=pca$x[,PCy], group=group, cond.df, name=colnames(object))
  pc1 <- ggplot(data=d, aes_string(x="PCa", y="PCb", color="group")) + geom_point(size = 1, shape = 15) + 
    scale_colour_brewer(type = "qual",palette = "Dark2")+ 
    xlab(paste0("PC",PCx,": ",round(percentVar[PCx] * 100),"% variance")) +
    ylab(paste0("PC",PCy,": ",round(percentVar[PCy] * 100),"% variance")) +
    coord_fixed()
  pc2 <- pc1 + geom_point()
  pc3 <- pc2 + theme_classic()
  #  Finally add the labels, using ggrepel so that they dont write over each other or the points  
  if (labels)
  {
    library("ggrepel")  
    pc3 + geom_text_repel(aes(label = name),
                          color = "gray20",
                          data = d,
                          force = 10)
  }else{ plot(pc3) }
}




# Convert to matrix and remove all negative values and convert them to integers.
countdata_mj609 <- as.matrix(countdata_mj609)


(condition <- factor(c(rep("wt", 3), rep("mj261", 3))))
(coldata <- data.frame(row.names=colnames(countdata_mj609), condition))

# convert to Deseq2 matrix
dds <- DESeqDataSetFromMatrix(countData=countdata_mj609, colData=coldata, design=~condition)
dds

# keep the comparison wt. 
dds$condition <- relevel(dds$condition, ref = "wt")
dds <- DESeq(dds, betaPrior = FALSE)

rld <- rlog(dds)

plotPCAEx(rld, 1,2,"condition",1000)
plotPCAEx(rld, 2,3,"condition",1000)

barplot(counts(dds, normalized = F), las = 2)

par(mfrow = c(2,2))
plot(log2(normCounts[,4] +1), log2(normCounts[,5] +1),pch = ".")
plot(log2(normCounts[,4] +1), log2(normCounts[,6] +1),pch = ".")
plot(log2(normCounts[,5] +1), log2(normCounts[,6] +1),pch = ".")
par(mfrow = c(1,1))


hist(assay(rld))
## print the results with log2fc and padj
res <- results(dds)
dds_normalized <- counts(dds, normalized=T)
dds_normalizeddf <- data.frame(counts(dds, normalized=T))
dds_normalizeddf$genes <- rownames(dds_normalizeddf)
res <- res[order(-res$log2FoldChange),]
head(res)
write.csv(dds_normalized, file="mj261_polyA_transposons_normalized_allreplicates.csv") 

resdata = as.data.frame(dplyr::mutate(as.data.frame(res), sig=ifelse(res$padj<0.01, "FDR<0.01", "Not Sig")), row.names=rownames(res))
resdata <- resdata %>% mutate(threshold = ifelse(log2FoldChange >=  1 & padj < 0.01,"A", ifelse(log2FoldChange <= -1 & padj < 0.01 , "B", "C")))
resdata <- resdata[order(-resdata$log2FoldChange),]
rownames(resdata) <- rownames(res)
resdata$genes <- rownames(resdata)
head(resdata)
write.csv(resdata, "mj261_polyA_transposons_diff.csv")

up <- subset(resdata, resdata$threshold == "A")
down <- subset(resdata, resdata$threshold == "B")


# volcano plot
p <- ggplot(resdata, aes(x=resdata$log2FoldChange, y=-log10(resdata$padj))) + 
  geom_point(color = "Grey", size=3) + geom_point(data=up, aes(x=up$log2FoldChange, y=-log10(up$padj)), color = "red", size=4.5) + geom_point(data=down, aes(x=down$log2FoldChange, y=-log10(down$padj)), color = "blue", size=4.5) +
  scale_x_continuous(limits = c(-10,10)) + 
  geom_vline(xintercept = 0, colour = "grey") + 
  theme_classic() + xlab("log2FoldChange") + ylab("-log10(padj)")

ggsave("mj261_polyA_transposons_volcano_allreplicates.pdf", p, height = 10, width = 10, useDingbats=FALSE)


# scatter plots for wt vs mut, use deseq2 normalized values

data <- read.csv("mj261_polyA_transposons_normalized_allreplicates.csv", header = TRUE, stringsAsFactors = FALSE)
head(data)


data$meanwt <- rowMeans(data[,2:4])
data$meanwt2 <- log10(rowMeans(data[,2:4])+1)

data$meanmj261 <- rowMeans(data[,5:7])
data$meanmj261_2 <- log10(rowMeans(data[,5:7])+1)


head(data)

datasig <- subset(data, data$X %in% sig5$X)

reg = lm(data$meanwt ~ data$meanmj261)
summary(reg)  #####0.402297    0.006389   62.97   <2e-16 ***

data2 <- subset(data, data$meanwt2 >= 1)
data3 <- subset(data2, data2$meanwt2 > (data2$meanmj261_2)*2)


g <- ggplot(data, aes(x=log2(data$meanwt+1), y=log2(data$meanmj261+1)))  + geom_point(stat = "identity") + geom_point(data=datasig, aes(x=log2(datasig$meanwt+1), y=log2(datasig$meanmj261+1)), color= "Blue", size=4)  + theme_classic() + geom_abline(slope = 0.9513 , intercept = 0.0001684 ) + xlab("wild type log2(TPM+1)") + ylab("rpb-9(mj261)V log2(TPM+1)")

ggsave("mj261_polyA_transposons_scatter_allreplicates.pdf", g, height = 10, width = 10, useDingbats=FALSE)

g2 <- ggplot(data, aes(x=log2(data$meanwt+1), y=log2(data$meanmj261+1))) + geom_point(data=datasig, aes(x=log2(datasig$meanwt+1), y=log2(datasig$meanmj261+1)), color= "Blue", size=4) + geom_point(stat = "identity")   + theme_classic() + geom_abline(slope = 0.9871 , intercept = 0.0006019 ) + xlab("wild type log2(TPM+1)") + ylab("rpb-9(mj261)V log2(TPM+1)")

ggsave("mj261_polyA_transposons_scatter_rep1.pdf", g2, height = 10, width = 10, useDingbats=FALSE)


######
#show diff genes in heatmap!filter the ones that overlap with diff protein coding genes.
######
setwd("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/RNA-seq_polyA/transposons/")

data <- read.csv("mj261_polyA_transposons_normalized.csv", header = TRUE, stringsAsFactors = FALSE)
head(data)

sig <- read.csv("mj261_polyA_transposons_diff.csv", header = TRUE, stringsAsFactors = TRUE)
head(sig)
sig <- as.data.frame(sig)

sig2 <- sig[abs(sig$log2FoldChange) >= 1 & sig$padj < 0.01, ]
sig3 <- na.omit(sig2)


sig4 <- sig3[c(13,17,22:35) , ]


data$meanwt <- rowMeans(data[,2:3])

datanorm <- data.frame(data$X, data[2:5]/(data$meanwt+1))
head(datanorm)
datanormsubset <- datanorm[datanorm$data.X %in% sig4$X,]
head(datanormsubset)

datanormsubset[,2:5] <- log2(datanormsubset[,2:5])

rownames(datanormsubset) <- datanormsubset$data.X


setwd("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/RNA-seq_polyA/transposons/all_replicates/")

data <- read.csv("mj261_polyA_transposons_normalized_allreplicates.csv", header = TRUE, stringsAsFactors = FALSE)
head(data)

sig <- read.csv("mj261_polyA_transposons_diff_allreplicates.csv", header = TRUE, stringsAsFactors = TRUE)
head(sig)
sig <- as.data.frame(sig)

sig2 <- sig[abs(sig$log2FoldChange) >= 1 & sig$padj < 0.01, ]
sig3 <- na.omit(sig2)


sig4 <- sig3[c(4,23,28,58:80) , ]



data$meanwt <- rowMeans(data[,2:4])

datanorm <- data.frame(data$X, data[2:7]/(data$meanwt+1))
head(datanorm)
datanormsubset <- datanorm[datanorm$data.X %in% sig4$X,]
head(datanormsubset)

datanormsubset[,2:7] <- log2(datanormsubset[,2:7])

rownames(datanormsubset) <- datanormsubset$data.X




sig5 <- sig3[c(4,23,28), ]



sig5 <- sig3[c(4,23,28), ]


data$meanwt <- rowMeans(data[,2:4])

datanorm <- data.frame(data$X, data[2:7]/(data$meanwt+1))
head(datanorm)
datanormsubset <- datanorm[datanorm$data.X %in% sig4$X,]
head(datanormsubset)

datanormsubset[,2:7] <- log2(datanormsubset[,2:7])

rownames(datanormsubset) <- datanormsubset$data.X


melted <- melt(datanormsubset)
head(melted)

# use ggplot2, this brings better vector graphics. geom_tile function is for the heatmap!

p <- ggplot(melted, aes(x = variable, y = data.X, fill = log2(value))) + geom_tile(color = "gray") + scale_fill_gradientn(colours = c("white","grey","black","red","red4"), 
                                                                                                              values = rescale(c(-6,-3,0,3,6)), guide = "colourbar") + 
  theme(axis.text.y = element_text(color="black", 
                                   
                                   size=2))


my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 299)
col_breaks = c(seq(-3,-1.5,length=100),  # for red
               seq(-1.49,1.49,length=100),           # for white
               seq(1.5,3,length=100))

heatmap.2(as.matrix(datanormsubset[2:7]),      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          col=my_palette,       # use on color palette defined earlier
          breaks=col_breaks,    # enable color transition at specified limits
          Colv = FALSE,
          Rowv = FALSE,
          dendrogram="none")     # only draw a row dendrogram 


#polyA libraries
# load featureCounts count.txt file. before loading, delete the first two rows for clarity.!!!textwrangler should help...
# filter too many reads from ribo and polyA libraries. file name is "ribo_highreads.csv"
setwd("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/RNA-seq_ribozero/transposons")

#introduce the countdata from featurecounts (fractional)
wt1 <- read.delim("SLX-15940.NEBNext01.HWFK7BBXX.s_5.r_1_ribo_mj261Aligned.sortedByCoord.out_mj261ribo_TE_unique.txt", header = FALSE, stringsAsFactors = FALSE)

wt2 <- read.delim("SLX-15940.NEBNext02.HWFK7BBXX.s_5.r_1.ribo_mj261Aligned.sortedByCoord.out_mj261ribo_TE_unique.txt", header = FALSE, stringsAsFactors = FALSE)

wt3 <- read.delim("SLX-15940.NEBNext03.HWFK7BBXX.s_5.r_1_ribo_mj261Aligned.sortedByCoord.out_mj261ribo_TE_unique.txt", header = FALSE, stringsAsFactors = FALSE)

mj261_1 <- read.delim("SLX-15940.NEBNext10.HWFK7BBXX.s_5.r_1_ribo_mj261Aligned.sortedByCoord.out_mj261ribo_TE_unique.txt", header = FALSE, stringsAsFactors = FALSE)

mj261_2 <- read.delim("SLX-15940.NEBNext11.HWFK7BBXX.s_5.r_1_ribo_mj261Aligned.sortedByCoord.out_mj261ribo_TE_unique.txt", header = FALSE, stringsAsFactors = FALSE)

mj261_3 <- read.delim("SLX-15940.NEBNext12.HWFK7BBXX.s_5.r_1_ribo_mj261Aligned.sortedByCoord.out_mj261ribo_TE_unique.txt", header = FALSE, stringsAsFactors = FALSE)


countdata <- data.frame(wt1[,c(4,7)], wt2[,7], wt3[,7], mj261_1[,7], mj261_2[,7], mj261_3[,7])

colnames(countdata) <- c("genes", "wt1", "wt2", "wt3", "mj261_1", "mj261_2", "mj261_3")
head(countdata)

countdata_mj609 <- countdata[,2:7]
rownames(countdata_mj609) <- countdata$genes
head(countdata_mj609)


plotPCAEx = function(object, PCx = 1, PCy = 2, cond="condition", ntop=500, labels = TRUE)
{
  library(RColorBrewer)
  # calculate the variance for each gene
  rv <- rowVars(assay(object))
  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(object)[select,]))
  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  if (!all(cond %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  cond.df <- as.data.frame(colData(object)[, cond, drop=FALSE])
  # add the intgroup factors together to create a new grouping factor
  group <- if (length(cond) > 1) {
    factor(apply( cond.df, 1, paste, collapse=" : "))
  } else {
    colData(object)[[cond]]
  }
  # assembly the data for the plot
  d <- data.frame(PCa=pca$x[,PCx], PCb=pca$x[,PCy], group=group, cond.df, name=colnames(object))
  pc1 <- ggplot(data=d, aes_string(x="PCa", y="PCb", color="group")) + geom_point(size = 1, shape = 15) + 
    scale_colour_brewer(type = "qual",palette = "Dark2")+ 
    xlab(paste0("PC",PCx,": ",round(percentVar[PCx] * 100),"% variance")) +
    ylab(paste0("PC",PCy,": ",round(percentVar[PCy] * 100),"% variance")) +
    coord_fixed()
  pc2 <- pc1 + geom_point()
  pc3 <- pc2 + theme_classic()
  #  Finally add the labels, using ggrepel so that they dont write over each other or the points  
  if (labels)
  {
    library("ggrepel")  
    pc3 + geom_text_repel(aes(label = name),
                          color = "gray20",
                          data = d,
                          force = 10)
  }else{ plot(pc3) }
}




# Convert to matrix and remove all negative values and convert them to integers.
countdata_mj609 <- as.matrix(countdata_mj609)


(condition <- factor(c(rep("wt", 3), rep("mj261", 3))))
(coldata <- data.frame(row.names=colnames(countdata_mj609), condition))

# convert to Deseq2 matrix
dds <- DESeqDataSetFromMatrix(countData=countdata_mj609, colData=coldata, design=~condition)
dds

# keep the comparison wt. 
dds$condition <- relevel(dds$condition, ref = "wt")
dds <- DESeq(dds, betaPrior = FALSE)

rld <- rlog(dds)

plotPCAEx(rld, 1,2,"condition",1000)
plotPCAEx(rld, 2,3,"condition",1000)

barplot(counts(dds, normalized = F), las = 2)

par(mfrow = c(2,2))
plot(log2(normCounts[,4] +1), log2(normCounts[,5] +1),pch = ".")
plot(log2(normCounts[,4] +1), log2(normCounts[,6] +1),pch = ".")
plot(log2(normCounts[,5] +1), log2(normCounts[,6] +1),pch = ".")
par(mfrow = c(1,1))


hist(assay(rld))
## print the results with log2fc and padj
res <- results(dds)
dds_normalized <- counts(dds, normalized=T)
dds_normalizeddf <- data.frame(counts(dds, normalized=T))
dds_normalizeddf$genes <- rownames(dds_normalizeddf)
res <- res[order(-res$log2FoldChange),]
head(res)
write.csv(dds_normalized, file="mj261_ribo_transposons_normalized.csv") 

resdata = as.data.frame(dplyr::mutate(as.data.frame(res), sig=ifelse(res$padj<0.01, "FDR<0.01", "Not Sig")), row.names=rownames(res))
resdata <- resdata %>% mutate(threshold = ifelse(log2FoldChange >=  1 & padj < 0.01,"A", ifelse(log2FoldChange <= -1 & padj < 0.01 , "B", "C")))
resdata <- resdata[order(-resdata$log2FoldChange),]
rownames(resdata) <- rownames(res)
resdata$genes <- rownames(resdata)
head(resdata)
write.csv(resdata, "mj261_ribo_transposons_diff.csv")

up <- subset(resdata, resdata$threshold == "A")
down <- subset(resdata, resdata$threshold == "B")


# volcano plot
p <- ggplot(resdata, aes(x=resdata$log2FoldChange, y=-log10(resdata$padj))) + 
  geom_point(color = "Grey", size=3) + geom_point(data=up, aes(x=up$log2FoldChange, y=-log10(up$padj)), color = "red", size=4.5) + geom_point(data=down, aes(x=down$log2FoldChange, y=-log10(down$padj)), color = "blue", size=4.5) +
  scale_x_continuous(limits = c(-10,10)) + 
  geom_vline(xintercept = 0, colour = "grey") + 
  theme_classic() + xlab("log2FoldChange") + ylab("-log10(padj)")

ggsave("mj261_ribo_transposons_volcano.pdf", p, height = 10, width = 10, useDingbats=FALSE)


# scatter plots for wt vs mut, use deseq2 normalized values

data <- read.csv("mj261_ribo_transposons_normalized.csv", header = TRUE, stringsAsFactors = FALSE)
head(data)


data$meanwt <- rowMeans(data[,2:4])
data$meanwt2 <- log10(rowMeans(data[,2:4])+1)

data$meanmj261 <- rowMeans(data[,5:7])
data$meanmj261_2 <- log10(rowMeans(data[,5:7])+1)


head(data)

datasig <- subset(data, data$X %in% sig5$X)

reg = lm(data$meanwt ~ data$meanmj261)
summary(reg)  #####0.402297    0.006389   62.97   <2e-16 ***

data2 <- subset(data, data$meanwt2 >= 1)
data3 <- subset(data2, data2$meanwt2 > (data2$meanmj261_2)*2)


g <- ggplot(data, aes(x=log2(data$meanwt+1), y=log2(data$meanmj261+1)))   + geom_point(data=datasig, aes(x=log2(datasig$meanwt+1), y=log2(datasig$meanmj261+1)), color= "Blue", size=4) + geom_point(stat = "identity") + theme_classic() + geom_abline(slope = 0.7527 , intercept = 0.0001684 ) + xlab("wild type log2(TPM+1)") + ylab("rpb-9(mj261)V log2(TPM+1)")

ggsave("mj261_ribo_transposons_scatter.pdf", g, height = 10, width = 10, useDingbats=FALSE)

g2 <- ggplot(data, aes(x=log2(data$meanwt+1), y=log2(data$meanmj261+1))) + geom_point(data=datasig, aes(x=log2(datasig$meanwt+1), y=log2(datasig$meanmj261+1)), color= "Blue", size=4) + geom_point(stat = "identity")   + theme_classic() + geom_abline(slope = 0.9871 , intercept = 0.0006019 ) + xlab("wild type log2(TPM+1)") + ylab("rpb-9(mj261)V log2(TPM+1)")

ggsave("mj261_polyA_transposons_scatter_rep1.pdf", g2, height = 10, width = 10, useDingbats=FALSE)


######
#show diff genes in heatmap!filter the ones that overlap with diff protein coding genes.
######

data <- read.csv("mj261_ribo_transposons_normalized.csv", header = TRUE, stringsAsFactors = FALSE)
head(data)

sig <- read.csv("mj261_ribo_transposons_diff.csv", header = TRUE, stringsAsFactors = TRUE)
head(sig)
sig <- as.data.frame(sig)

sig2 <- sig[abs(sig$log2FoldChange) >= 1 & sig$padj < 0.01, ]
sig3 <- na.omit(sig2)


sig4 <- sig3[c(27,28,44:63) , ]
sig5 <- sig3[c(27,28) , ]


data$meanwt <- rowMeans(data[,2:4])

datanorm <- data.frame(data$X, data[2:7]/(data$meanwt+1))
head(datanorm)
datanormsubset <- datanorm[datanorm$data.X %in% sig4$X,]
head(datanormsubset)

datanormsubset[,2:7] <- log2(datanormsubset[,2:7])

rownames(datanormsubset) <- datanormsubset$data.X


setwd("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/RNA-seq_polyA/transposons/all_replicates/")

data <- read.csv("mj261_polyA_transposons_normalized_allreplicates.csv", header = TRUE, stringsAsFactors = FALSE)
head(data)

sig <- read.csv("mj261_polyA_transposons_diff_allreplicates.csv", header = TRUE, stringsAsFactors = TRUE)
head(sig)
sig <- as.data.frame(sig)

sig2 <- sig[abs(sig$log2FoldChange) >= 1 & sig$padj < 0.01, ]
sig3 <- na.omit(sig2)


sig4 <- sig3[c(4,23,28,58:80) , ]



data$meanwt <- rowMeans(data[,2:4])

datanorm <- data.frame(data$X, data[2:7]/(data$meanwt+1))
head(datanorm)
datanormsubset <- datanorm[datanorm$data.X %in% sig4$X,]
head(datanormsubset)

datanormsubset[,2:7] <- log2(datanormsubset[,2:7])

rownames(datanormsubset) <- datanormsubset$data.X




sig5 <- sig3[c(4,23,28), ]



sig5 <- sig3[c(4,23,28), ]


data$meanwt <- rowMeans(data[,2:4])

datanorm <- data.frame(data$X, data[2:7]/(data$meanwt+1))
head(datanorm)
datanormsubset <- datanorm[datanorm$data.X %in% sig4$X,]
head(datanormsubset)

datanormsubset[,2:7] <- log2(datanormsubset[,2:7])

rownames(datanormsubset) <- datanormsubset$data.X


melted <- melt(datanormsubset)
head(melted)

# use ggplot2, this brings better vector graphics. geom_tile function is for the heatmap!

p <- ggplot(melted, aes(x = variable, y = data.X, fill = log2(value))) + geom_tile(color = "gray") + scale_fill_gradientn(colours = c("white","grey","black","red","red4"), 
                                                                                                                          values = rescale(c(-6,-3,0,3,6)), guide = "colourbar") + 
  theme(axis.text.y = element_text(color="black", 
                                   
                                   size=2))


my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 299)
col_breaks = c(seq(-3,-1.5,length=100),  # for red
               seq(-1.49,1.49,length=100),           # for white
               seq(1.5,3,length=100))

heatmap.2(as.matrix(datanormsubset[2:7]),      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          col=my_palette,       # use on color palette defined earlier
          breaks=col_breaks,    # enable color transition at specified limits
          Colv = FALSE,
          Rowv = FALSE,
          dendrogram="none")     # only draw a row dendrogram 




