# this script was created by Ahmet Can Berkyurek on 13/07/18, to analyze AMA1 enrichment on transposable elements. 
library(GenomicAlignments)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library("genefilter")
library(pheatmap)
library(SummarizedExperiment)
library(gdata)
library(dplyr)
library(tidyr)
library(scales)
library(reshape2)
# load featureCounts count.txt file. before loading, delete the first two rows for clarity.!!!textwrangler should help...
setwd("/Volumes/miska/Ahmet Berkyurek/sequencing/NGS/ChIP-seq/AMA-1/ce11")

countdata <- read.table("/Volumes/miska/Ahmet Berkyurek/sequencing/NGS/ChIP-seq/AMA-1/ce11/AMA1_ChIP_TEcounts.txt", header=TRUE, row.names=1)
head(countdata)
#All counts
# Remove first five columns (chr, start, end, strand, length)
countdata <- countdata[ ,5:ncol(countdata)]

# Remove .bam or .sam from filenames
colnames(countdata) <- gsub("\\.[sb]am$", "", colnames(countdata))

# Convert to matrix and remove all negative values and convert them to integers.
countdata <- as.matrix(countdata)
head(countdata)

countdataframe <- data.frame(countdata)
head(countdataframe)
colnames(countdataframe) <- c('Gene_length', 'wt_1', 'wt_2', 'wt_3', 'mj261_1','mj261_2', 'mj261_3')

## determine library size factor
## wt replicates
wt1 <- countBam("/Volumes/miska/Ahmet\ Berkyurek/sequencing/NGS/ChIP-seq/AMA-1/ce11/wt/SX1316_ChIP_ce11.sorted.bam")

wt1length <- wt1$records
head(wt1length)


wt2 <- countBam("/Volumes/miska/Ahmet\ Berkyurek/sequencing/NGS/ChIP-seq/AMA-1/ce11/wt/SX1316_ChIP_ce11_rep2.sorted.bam")

wt2length <- wt2$records
head(wt2length)

wt3 <- countBam("/Volumes/miska/Ahmet\ Berkyurek/sequencing/NGS/ChIP-seq/AMA-1/ce11/wt/SX1316_ChIP_ce11_rep3.sorted.bam")

wt3length <- wt3$records
head(wt3length)

## mj261 replicates
mj261_1 <- countBam("/Volumes/miska/Ahmet\ Berkyurek/sequencing/NGS/ChIP-seq/AMA-1/ce11/mj261/SX1984_ChIP_ce11.sorted.bam")

mj261_1length <- mj261_1$records
head(mj261_1length)


mj261_2 <- countBam("/Volumes/miska/Ahmet\ Berkyurek/sequencing/NGS/ChIP-seq/AMA-1/ce11/mj261/SX1984_ChIP_ce11_rep2.sorted.bam")

mj261_2length <- mj261_2$records
head(mj261_2length)

mj261_3 <- countBam("/Volumes/miska/Ahmet\ Berkyurek/sequencing/NGS/ChIP-seq/AMA-1/ce11/mj261/SX1984_ChIP_ce11_rep3.sorted.bam")

mj261_3length <- mj261_3$records
head(mj261_3length)

## determine individual factors by dividing with the average

averagelibrarysize <- mean(wt1length,wt2length,wt3length,mj261_1length,mj261_2length,mj261_3length)
print(averagelibrarysize)

wt1factor  = wt1length/averagelibrarysize
print(wt1factor)

wt2factor = wt2length/averagelibrarysize
print(wt2factor)

wt3factor = wt3length/averagelibrarysize
print(wt3factor)

mj261_1factor = mj261_1length/averagelibrarysize
print(mj261_1factor)

mj261_2factor = mj261_2length/averagelibrarysize
print(mj261_2factor)

mj261_3factor = mj261_3length/averagelibrarysize
print(mj261_3factor)

## divide each count column with individual size factors to overcome the sequencing bias!

countdataframenrm <- data.frame(countdataframe$wt_1/wt1factor,countdataframe$wt_2/wt2factor,countdataframe$wt_3/wt3factor,countdataframe$mj261_1/mj261_1factor,countdataframe$mj261_2/mj261_2factor,countdataframe$mj261_3/mj261_3factor)
rownames(countdataframenrm) = rownames(countdataframe)
colnames(countdataframenrm) = c('wt1','wt2','wt3','mj261_1','mj261_2','mj261_3')
head(countdataframenrm)


## divide each column with average wt counts


countdataframenrm <- countdataframenrm[1:6+0.1]
countdataframenrm$meanwt <- rowMeans(countdataframenrm[4:6])
head(countdataframenrm)

countdataframenrmfinal <- countdataframenrm[1:6]/countdataframenrm$meanwt
head(countdataframenrmfinal)
countdataframenrmfinal <- countdataframenrmfinal[order(-countdataframenrmfinal$wt2),]

countdataframenrmfinal$genes <- rownames(countdataframenrmfinal)                                                  
head(countdataframenrmfinal)


melted <- countdataframenrmfinal %>% gather(sample, FE, -genes)
head(melted)

# use ggplot2, this brings better vector graphics. geom_tile function is for the heatmap!

p <- ggplot(melted, aes(x = sample, y = genes, fill = FE)) + geom_tile(color = "gray") + scale_fill_gradientn(colours = c("white","grey","black","red","red4"), 
                                                                                                                      values = rescale(c(0,2.5,5,7.5,10)), guide = "colourbar") + 
  theme(axis.text.y = element_text(color="black", 
                                   
                                   size=2))


ggsave("snRNAs_heatmap_newnames_all.pdf",p, width = 6, height = 12)


# select the first 25 rows!

top <- countdataframenrmfinal[1:25,]
write.csv(top$genes, 'top25.csv')
topmelted <- top %>% gather(sample2, FE2,-genes)

p <- ggplot(topmelted, aes(x = sample2, y = genes, fill = FE2)) + geom_tile(color = "gray") + scale_fill_gradientn(colours = c("white","grey","black","red","red4"), 
                                                                                                              values = rescale(c(0,2.5,5,7.5,10)), guide = "colourbar") + 
  theme(axis.text.y = element_text(color="black", 
                                   
                                   size=12))


ggsave("snRNAs_heatmap_newnames_all.pdf",p, width = 6, height = 12)



