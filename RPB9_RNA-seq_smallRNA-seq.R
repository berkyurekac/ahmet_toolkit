# this script was created by Ahmet Can Berkyurek on 10/04/2020, to analyze differential gene expression analysis for RPB- total and small RNA
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


# introduce small RNA raw data and genes from RNA-seq/ChIPseq correlation
setwd("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/Figure_6")

smallRNA <- read.csv("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/RNA-seq_correlation/raw_counts_smallRNAss.csv", header = TRUE, stringsAsFactors = FALSE)
head(smallRNA)

smallRNA$wt1 <- (smallRNA$wt1/sum(smallRNA$wt1))*1e6
smallRNA$wt2 <- (smallRNA$wt2/sum(smallRNA$wt2))*1e6
smallRNA$wt3 <- (smallRNA$wt3/sum(smallRNA$wt3))*1e6

smallRNA$mutant1 <- (smallRNA$mutant1/sum(smallRNA$mutant1))*1e6
smallRNA$mutant2 <- (smallRNA$mutant2/sum(smallRNA$mutant2))*1e6
smallRNA$mutant3 <- (smallRNA$mutant3/sum(smallRNA$mutant3))*1e6

smallRNA$rescue1 <- (smallRNA$rescue1/sum(smallRNA$rescue1))*1e6
smallRNA$rescue2 <- (smallRNA$rescue2/sum(smallRNA$rescue2))*1e6
smallRNA$rescue3 <- (smallRNA$rescue3/sum(smallRNA$rescue3))*1e6

head(smallRNA)


smallRNA$meanwt <- rowMeans(smallRNA[,2:4])

smallRNA[,2:11] <- log2((smallRNA[,2:11]+1)/(smallRNA$meanwt+1))
rownames(smallRNA) <- smallRNA$X
head(smallRNA)


#polyA stringent

upclassI <- read.delim("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/Figure_4/subsets/upclassI.tsv", header = FALSE, stringsAsFactors = FALSE)

upclassII <- read.delim("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/Figure_4/subsets/upclassII.tsv", header = FALSE, stringsAsFactors = FALSE)

upclassIII <- read.delim("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/Figure_4/subsets/upclassIII.tsv", header = FALSE, stringsAsFactors = FALSE)


#subset with smallRNA 

upclassIsmallRNA <- subset(smallRNA, smallRNA$X %in% upclassI$V1)

upclassIIsmallRNA <- subset(smallRNA, smallRNA$X %in% upclassII$V1)

upclassIIIsmallRNA <- subset(smallRNA, smallRNA$X %in% upclassIII$V1)

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 299)
col_breaks = c(seq(-5,-1.5,length=100),  # for red
               seq(-1.49,1.49,length=100),           # for white
               seq(1.5,5,length=100))

heatmap.2(as.matrix(smallRNA[2:10]),      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          col=my_palette,       # use on color palette defined earlier
          breaks=col_breaks,    # enable color transition at specified limits
          Colv = FALSE,
          Rowv = FALSE,
          dendrogram="none")     # only draw a row dendrogram 




#polyA all_replicates upclassI_new.tsv

upclassI <- read.delim("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/Figure_4/subsets/polyA_allreplicates/upclassI.tsv", header = FALSE, stringsAsFactors = FALSE)

upclassInew <- read.delim("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/Figure_4/subsets/polyA_allreplicates/upclassI_new.tsv", header = FALSE, stringsAsFactors = FALSE)

upclassII <- read.delim("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/Figure_4/subsets/polyA_allreplicates/upclassII.tsv", header = FALSE, stringsAsFactors = FALSE)

upclassIII <- read.delim("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/Figure_4/subsets/polyA_allreplicates/upclassIII.tsv", header = FALSE, stringsAsFactors = FALSE)

sperm <- read.delim("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/RNA-seq_correlation/sperm_genes_reinke.bed", header = FALSE)

#subset with smallRNA 

upclassIsmallRNA <- subset(smallRNA, smallRNA$X %in% upclassI$V1)

upclassInewsmallRNA <- subset(smallRNA, smallRNA$X %in% upclassInew$V1)


upclassIIsmallRNA <- subset(smallRNA, smallRNA$X %in% upclassII$V1)

upclassIIIsmallRNA <- subset(smallRNA, smallRNA$X %in% upclassIII$V1)

spermsmallRNA <- subset(smallRNA, smallRNA$X %in% sperm$V1)

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 299)
col_breaks = c(seq(-5,-1.5,length=100),  # for red
               seq(-1.49,1.49,length=100),           # for white
               seq(1.5,5,length=100))

heatmap.2(as.matrix(upclassInewsmallRNA[2:10]),      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          col=my_palette,       # use on color palette defined earlier
          breaks=col_breaks,    # enable color transition at specified limits
          Colv = FALSE,
          Rowv = TRUE,
          dendrogram="none")     # only draw a row dendrogram 




