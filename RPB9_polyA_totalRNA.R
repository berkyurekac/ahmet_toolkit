# this script was created by Ahmet Can Berkyurek on 01/04/2020, to analyze differential gene expression analysis for mj261 rpb-9 mutant
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
setwd("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/RNA-seq_polyA/all_replicates")

#introduce the countdata from featurecounts (fractional)
countdata <- read.table("mj261_polyA_total_ce11.txt", header=TRUE, row.names=1)

# Remove first five columns (chr, start, end, strand, length)
countdata <- countdata[ ,6:ncol(countdata)]
countdata <- data.frame(countdata)
head(countdata)

colnames(countdata) <- c("wt1", "wt2", "wt3", "mj261_1", "mj261_2", "mj261_3")
head(countdata)


countdata_mj609 <- countdata[,c(1,2,3,4,5,6)]
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

countdata_mj609 <- round(countdata_mj609)
head(countdata_mj609)
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
write.csv(dds_normalized, file="mj261_polyA_normalized_allreplicates.csv") 

resdata = as.data.frame(dplyr::mutate(as.data.frame(res), sig=ifelse(res$padj<0.01, "FDR<0.01", "Not Sig")), row.names=rownames(res))
resdata <- resdata %>% mutate(threshold = ifelse(log2FoldChange >= 1 & padj < 0.01,"A", ifelse(log2FoldChange <= -1 & padj < 0.01 , "B", "C")))
resdata <- resdata[order(-resdata$log2FoldChange),]
rownames(resdata) <- rownames(res)
resdata$genes <- rownames(resdata)
head(resdata)

osc <- read.csv("/Volumes/miska/Ahmet Berkyurek/sequencing/NGS/RNA-seq/Lisa_ribozero/mapped/all_genes/oscilatinggenes.tsv")
osc <- data.frame(osc)
head(osc)

diff_ext <- subset(resdata, !(rownames(resdata) %in% osc$WBGene00000002))
head(diff_ext)
write.csv(diff_ext, "mj261_polyA_diff_osc_allreplicates.csv")

non <- subset(diff_ext, diff_ext$threshold == "C")
write.csv(non$genes, "non-regulated.csv")


up <- subset(diff_ext, diff_ext$threshold == "A")
down <- subset(diff_ext, diff_ext$threshold == "B")
Sensor <- subset(diff_ext, diff_ext$genes == "piRNASensor")
piRNAgenes <- read.delim("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/Final Figures/EMBO/revision/bioinformatics/piRNA_pathway_genes/piRNApathwaygenes.tsv", header = F)
piRNAgenesset <- subset(diff_ext, diff_ext$genes %in% piRNAgenes$V1)
write.csv(piRNAgenesset, "piRNAgenes_diff.csv")

setwd("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/Final Figures/EMBO/revision/bioinformatics/piRNA_pathway_genes/")

# volcano plot
p <- ggplot(diff_ext, aes(x=diff_ext$log2FoldChange, y=-log10(diff_ext$padj))) + 
  geom_point(color = "Grey", size=3) + geom_point(data=up, aes(x=up$log2FoldChange, y=-log10(up$padj)), color = "red", size=4.5)+ geom_point(data=piRNAgenesset, aes(x=piRNAgenesset$log2FoldChange, y=-log10(piRNAgenesset$padj)), color = "yellow", size=1.5) + geom_point(data=Sensor, aes(x=Sensor$log2FoldChange, y=-log10(Sensor$padj)), color = "darkgreen", size=6)  + geom_point(data=down, aes(x=down$log2FoldChange, y=-log10(down$padj)), color = "blue", size=4.5) +
  scale_x_continuous(limits = c(-10,10)) + 
  geom_vline(xintercept = 0, colour = "grey") + 
  theme_classic() + xlab("log2FoldChange") + ylab("-log10(padj)")

ggsave("mj261_polyA_osc_volcano_allreplicates_piRNAgenes.pdf", p, height = 15, width = 10, useDingbats=FALSE)


#####
# check for specific gene targets - germline -sperm - soma- neuron - hrde1 targets
######
setwd("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/RNA-seq_correlation")


shared <- read.delim("Fig2bsharedgermline geneset.txt", header = TRUE, stringsAsFactors = FALSE)
head(shared)

write.csv(shared$primerpair, "oocyte_sperm_commongenes.csv")

spermoocytecommon <- read.delim("sperm_oocyte_commongenes.tsv", header = FALSE)

spermreinke <- read.csv("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/RNA-seq_correlation/sperm_genes_reinke.bed", header = FALSE, stringsAsFactors = FALSE)

spermortiz <- read.csv("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/RNA-seq_correlation/sperm_genes_ortiz.bed", header = FALSE, stringsAsFactors = FALSE)

oocytereinke <- read.csv("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/RNA-seq_correlation/oogenesis_genes_reinke.bed", header = FALSE, stringsAsFactors = FALSE)

oocyteortiz <- read.csv("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/RNA-seq_correlation/oocyte_genes_ortiz.bed", header = FALSE, stringsAsFactors = FALSE)

somareinke <- read.csv("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/RNA-seq_correlation/soma_genes_reinke.bed", header = FALSE, stringsAsFactors = FALSE)

neuron <- read.csv("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/RNA-seq_correlation/neuron_genes_watson.bed", header = FALSE, stringsAsFactors = FALSE)

hrde1targets <- read.csv("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/RNA-seq_correlation/hrde1_allgenes_osc.csv", header = FALSE)

CSR1targets <- read.csv("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/RNA-seq_correlation/CSR1_IP_smallRNAtargets_conine.tsv", header = FALSE, stringsAsFactors = FALSE)

ALG <- read.csv("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/RNA-seq_correlation/ALG3_4_targets_genes_conine.tsv", header = FALSE, stringsAsFactors = FALSE)

strome <- read.csv("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/RNA-seq_correlation/Strome_2017_PNAS_genes.csv", header = TRUE, stringsAsFactors = FALSE)

stromesperm <- subset(strome, strome$sperm.genes == "1")

stromesoma <- subset(strome, strome$soma.specific.genes == "1")

stromeubiq <- subset(strome, strome$ubiquitous.genes == "1")

fabiangermline <- read.delim("Fabiangermline.tsv", header = FALSE)
fabiangermline1 <- subset(fabiangermline, !(fabiangermline$V1) %in% spermreinke$V1)
fabiangermline2 <- subset(fabiangermline1, !(fabiangermline1$V1) %in% oocytereinke$V1)


mj609spermreinke <- subset(diff_ext, diff_ext$genes %in% spermreinke$V1)

spermsig <- subset(up, up$genes %in% spermreinke$V1)

mj609spermortiz <- subset(diff_ext, diff_ext$genes %in% spermortiz$V1)

mj609oocytereinke <- subset(diff_ext, diff_ext$genes %in% oocytereinke$V1)

mj609oocyeortiz <- subset(diff_ext, diff_ext$genes %in% oocyteortiz$V1)

mj609soma <- subset(diff_ext, diff_ext$genes %in% somareinke$V1)

mj609neuron <- subset(diff_ext, diff_ext$genes %in% neuron$V1)

mj609CSR1targets <- subset(diff_ext, diff_ext$genes %in% CSR1targets$V1)

mj609hrde1targets <- subset(diff_ext, diff_ext$genes %in% hrde1targets$V1)

mj609ALG <- subset(diff_ext, diff_ext$genes %in% ALG$V1)

mj609stromesperm <- subset(diff_ext, diff_ext$genes %in% stromesperm$Wormbase.ID)

mj609stromesoma <- subset(diff_ext, diff_ext$genes %in% stromesoma$Wormbase.ID)

mj609stromaubiq <- subset(diff_ext, diff_ext$genes %in% stromeubiq$Wormbase.ID)

Nosoma <- subset(diff_ext, !(diff_ext$genes) %in% mj609soma$genes)

mj609fabiangermline <- subset(diff_ext, diff_ext$genes %in% fabiangermline2$V1)



#####   
# label diff genes with tissue soecific names, then box plot!
##### 

sets <- diff_ext %>%
  mutate(Label = case_when(diff_ext$genes %in% mj609spermreinke$genes ~ "A",
                           diff_ext$genes %in% mj609oocytereinke$genes ~ "B",
                           diff_ext$genes %in% mj609soma$genes ~ "C",
                           TRUE ~ "M"))


p <- ggplot(sets, aes(x=sets$Label, y=sets$log2FoldChange)) + geom_boxplot( color="Black", size=1) + theme_classic() + scale_x_discrete(limits=c("A", "B", "C","M")) + scale_y_continuous(limits = c(-5,5))

ggsave("mj261_polyA_genesets_boxplot_allreplicates.pdf", p, height = 10, width = 10, useDingbats=FALSE)


sets2 <- diff_ext %>%
  mutate(Label = case_when(diff_ext$genes %in% mj609stromesperm$genes ~ "A",
                           diff_ext$genes %in% mj609stromesoma$genes ~ "B",
                           diff_ext$genes %in% mj609stromaubiq$genes ~ "C",
                           TRUE ~ "M"))

p <- ggplot(sets2, aes(x=sets2$Label, y=sets2$log2FoldChange)) + geom_boxplot( color="Black", size=1) + theme_classic() + scale_x_discrete(limits=c("A", "B", "C","M")) + scale_y_continuous(limits = c(-5,5))

ggsave("mj261_polyA_genesets_boxplot.pdf", p, height = 10, width = 10, useDingbats=FALSE)

#####
# for polyA, plot all genes substracted with soma, only soma, and all genes! diff_ext$genes %in% mj609soma$genes ~ "A",
######

sets <- diff_ext %>%
  mutate(Label = case_when(diff_ext$genes %in% mj609fabiangermline$genes  ~ "C",
                           # diff_ext$genes %in% Nosoma$genes ~ "B",
                           # diff_ext$genes %in% mj609fabiangermline$genes  ~ "C",
                           FALSE ~ "M"))


p <- ggplot(sets, aes(x=sets$Label, y=sets$log2FoldChange)) + geom_boxplot( color="Black", size=1) + theme_classic() + scale_x_discrete(limits=c("A","B","C")) + scale_y_continuous(limits = c(-5,5))

ggsave("X5.pdf", p, height = 10, width = 10, useDingbats=FALSE)

#####
# for polyA, then label sperm, oocyte, all germline
######

sets <- diff_ext %>%
  mutate(Label = case_when(diff_ext$genes %in% mj609spermreinke$genes ~ "A",
                           diff_ext$genes %in% mj609oocytereinke$genes ~ "B",
                           diff_ext$genes %in% mj609fabiangermline$genes  ~ "C",
                           TRUE ~ "M"))


p <- ggplot(sets, aes(x=sets$Label, y=sets$log2FoldChange)) + geom_boxplot( color="Black", size=1) + theme_classic() + scale_x_discrete(limits=c("A", "B","M")) + scale_y_continuous(limits = c(-5,5))

ggsave("mj261_ribo_genesets_boxplot_germline_sperm_oocyte_therest_allreplicates.pdf", p, height = 10, width = 10, useDingbats=FALSE)



#ribozero libraries
# load featureCounts count.txt file. before loading, delete the first two rows for clarity.!!!textwrangler should help...
# filter too many reads from ribo and polyA libraries. file name is "ribo_highreads.csv"
setwd("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/RNA-seq_ribozero")

#introduce the countdata from featurecounts (fractional)
countdata <- read.table("counts_Lisa-ribozero.txt", header=TRUE, row.names=1)

# Remove first five columns (chr, start, end, strand, length)
countdata <- countdata[ ,6:ncol(countdata)]
countdata <- data.frame(countdata)
head(countdata)

colnames(countdata) <- c("wt1", "wt2", "wt3", "mj261_1", "mj261_2", "mj261_3")
head(countdata)

countdata_mj609 <- countdata
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

countdata_mj609 <- round(countdata_mj609)
head(countdata_mj609)
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
write.csv(dds_normalized, file="mj261_ribozero_normalized.csv") 

resdata = as.data.frame(dplyr::mutate(as.data.frame(res), sig=ifelse(res$padj<0.01, "FDR<0.01", "Not Sig")), row.names=rownames(res))
resdata <- resdata %>% mutate(threshold = ifelse(log2FoldChange >= 1 & padj < 0.01,"A", ifelse(log2FoldChange <= -1 & padj < 0.01 , "B", "C")))
resdata <- resdata[order(-resdata$log2FoldChange),]
rownames(resdata) <- rownames(res)
resdata$genes <- rownames(resdata)
head(resdata)

osc <- read.csv("/Volumes/miska/Ahmet Berkyurek/sequencing/NGS/RNA-seq/Lisa_ribozero/mapped/all_genes/oscilatinggenes.tsv")
osc <- data.frame(osc)
head(osc)

diff_ext <- subset(resdata, !(rownames(resdata) %in% osc$WBGene00000002))
head(diff_ext)
write.csv(diff_ext, "mj261_ribozero_diff_osc.csv")

up <- subset(diff_ext, diff_ext$threshold == "A")
down <- subset(diff_ext, diff_ext$threshold == "B")
Sensor <- subset(diff_ext, diff_ext$genes == "piRNASensor")


# volcano plot
p <- ggplot(diff_ext, aes(x=diff_ext$log2FoldChange, y=-log10(diff_ext$padj))) + 
  geom_point(color = "Grey", size=3) + geom_point(data=up, aes(x=up$log2FoldChange, y=-log10(up$padj)), color = "red", size=4.5) + geom_point(data=Sensor, aes(x=Sensor$log2FoldChange, y=-log10(Sensor$padj)), color = "darkgreen", size=6)  + geom_point(data=down, aes(x=down$log2FoldChange, y=-log10(down$padj)), color = "blue", size=4.5) +
  scale_x_continuous(limits = c(-10,10)) + 
  geom_vline(xintercept = 0, colour = "grey") + 
  theme_classic() + xlab("log2FoldChange") + ylab("-log10(padj)")

ggsave("mj261_ribozero_osc_volcano.pdf", p, height = 15, width = 10, useDingbats=FALSE)


#####
# check for specific gene targets - germline -sperm - soma- neuron - hrde1 targets
######
setwd("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/RNA-seq_correlation")

shared <- read.delim("Fig2bsharedgermline geneset.txt", header = TRUE, stringsAsFactors = FALSE)
head(shared)

write.csv(shared$primerpair, "oocyte_sperm_commongenes.csv")

spermoocytecommon <- read.delim("sperm_oocyte_commongenes.tsv", header = FALSE)

spermreinke <- read.csv("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/RNA-seq_correlation/sperm_genes_reinke.bed", header = FALSE, stringsAsFactors = FALSE)

spermortiz <- read.csv("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/RNA-seq_correlation/sperm_genes_ortiz.bed", header = FALSE, stringsAsFactors = FALSE)

oocytereinke <- read.csv("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/RNA-seq_correlation/oogenesis_genes_reinke.bed", header = FALSE, stringsAsFactors = FALSE)

oocyteortiz <- read.csv("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/RNA-seq_correlation/oocyte_genes_ortiz.bed", header = FALSE, stringsAsFactors = FALSE)

somareinke <- read.csv("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/RNA-seq_correlation/soma_genes_reinke.bed", header = FALSE, stringsAsFactors = FALSE)

neuron <- read.csv("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/RNA-seq_correlation/neuron_genes_watson.bed", header = FALSE, stringsAsFactors = FALSE)

hrde1targets <- read.csv("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/RNA-seq_correlation/hrde1_allgenes_osc.csv", header = FALSE)

CSR1targets <- read.csv("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/RNA-seq_correlation/CSR1_IP_smallRNAtargets_conine.tsv", header = FALSE, stringsAsFactors = FALSE)

ALG <- read.csv("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/RNA-seq_correlation/ALG3_4_targets_genes_conine.tsv", header = FALSE, stringsAsFactors = FALSE)

strome <- read.csv("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/RNA-seq_correlation/Strome_2017_PNAS_genes.csv", header = TRUE, stringsAsFactors = FALSE)

stromesperm <- subset(strome, strome$sperm.genes == "1")

stromesoma <- subset(strome, strome$soma.specific.genes == "1")

stromeubiq <- subset(strome, strome$ubiquitous.genes == "1")


piRNAtargets <- read.csv("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/RNA-seq_correlation/piRNAtargets.txt", header = FALSE, stringsAsFactors = FALSE)



mj609spermreinke <- subset(diff_ext, diff_ext$genes %in% spermreinke$V1)

mj609spermortiz <- subset(diff_ext, diff_ext$genes %in% spermortiz$V1)

mj609oocytereinke <- subset(diff_ext, diff_ext$genes %in% oocytereinke$V1)

mj609oocyeortiz <- subset(diff_ext, diff_ext$genes %in% oocyteortiz$V1)

mj609soma <- subset(diff_ext, diff_ext$genes %in% somareinke$V1)

mj609neuron <- subset(diff_ext, diff_ext$genes %in% neuron$V1)

mj609CSR1targets <- subset(diff_ext, diff_ext$genes %in% CSR1targets$V1)

mj609hrde1targets <- subset(diff_ext, diff_ext$genes %in% hrde1targets$V1)

mj609ALG <- subset(diff_ext, diff_ext$genes %in% ALG$V1)

mj609stromesperm <- subset(diff_ext, diff_ext$genes %in% stromesperm$Wormbase.ID)

mj609stromesoma <- subset(diff_ext, diff_ext$genes %in% stromesoma$Wormbase.ID)

mj609stromaubiq <- subset(diff_ext, diff_ext$genes %in% stromeubiq$Wormbase.ID)

mj609piRNatargets <- subset(diff_ext, diff_ext$genes %in% piRNAtargets$V1)


Nosoma <- subset(diff_ext, !(diff_ext$genes) %in% mj609soma$genes)

mj609fabiangermline <- subset(diff_ext, diff_ext$genes %in% fabiangermline2$V1)


#####   
# label diff genes with tissue soecific names, then box plot!
##### 

sets <- diff_ext %>%
  mutate(Label = case_when(diff_ext$genes %in% mj609spermreinke$genes ~ "A",
                           diff_ext$genes %in% mj609oocytereinke$genes ~ "B",
                           #diff_ext$genes %in% mj609soma$genes ~ "C",
                           #diff_ext$genes %in% mj609neuron$genes ~ "D",
                           #diff_ext$genes %in% mj609piRNatargets$genes ~ "E",
                           TRUE ~ "M"))


p <- ggplot(sets, aes(x=sets$Label, y=sets$log2FoldChange)) + geom_boxplot( color="Black", size=1) + theme_classic() + scale_x_discrete(limits=c("A", "B", "M")) + scale_y_continuous(limits = c(-5,5))

ggsave("mj261_ribozero_genesets_boxplot2.pdf", p, height = 10, width = 10, useDingbats=FALSE)


sets <- diff_ext %>%
  mutate(Label = case_when(diff_ext$genes %in% mj609spermreinke$genes ~ "A",
                           diff_ext$genes %in% mj609oocytereinke$genes ~ "B",
                           diff_ext$genes %in% mj609fabiangermline$genes  ~ "C",
                           FALSE ~ "M"))


p <- ggplot(sets2, aes(x=sets2$Label, y=sets2$log2FoldChange)) + geom_boxplot( color="Black", size=1) + theme_classic() + scale_x_discrete(limits=c("A", "B", "C","M")) + scale_y_continuous(limits = c(-5,5))

ggsave("mj261_polyA_genesets_boxplot.pdf", p, height = 10, width = 10, useDingbats=FALSE)


#####   
# scatter plot for diff genes
##### 


setwd("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/RNA-seq_polyA/")

polyA <- read.csv("mj261_polyA_diff_osc.csv", header = TRUE)
polyAup <- subset(polyA, polyA$threshold == "A")
polyAdown <- subset(polyA, polyA$threshold == "B")

polyAnorm <- read.csv("mj261_polyA_normalized.csv", header = TRUE)

head(polyAnorm)

polyAnorm$meanwt <- rowMeans(polyAnorm[,2:4])
polyAnorm$meanmj261 <- rowMeans(polyAnorm[,5:6])

polyAupsig <- subset(polyAnorm, polyAnorm$X %in% polyAup$X)
polyAdownsig <- subset(polyAnorm, polyAnorm$X %in% polyAdown$X)
sensor <- subset(polyAnorm, polyAnorm$X == "piRNASensor")


reg = lm(polyAnorm$meanwt ~ polyAnorm$meanmj261)
summary(reg)

g <- ggplot(polyAnorm, aes(x=log2(polyAnorm$meanwt+1), y=log2(polyAnorm$meanmj261+1))) + geom_point(data=polyAupsig, aes(x=log2(polyAupsig$meanwt+1), y=log2(polyAupsig$meanmj261+1)), color= "Red", size=4) + geom_point(data=sensor, aes(x=log2(sensor$meanwt+1), y=log2(sensor$meanmj261+1)), color= "Darkgreen", size=4)   + geom_point(data=polyAdownsig, aes(x=log2(polyAdownsig$meanwt+1), y=log2(polyAdownsig$meanmj261+1)), color= "Blue", size=4) + geom_point(stat = "identity")  + theme_classic() + geom_abline(slope = 0.6541 , intercept = 0.016 ) + xlab("wild type log2(TPM+1)") + ylab("rpb-9(mj261)V log2(TPM+1)")

ggsave("mj261_polyA_transposons_scatter_allreplicates.pdf", g, height = 10, width = 10, useDingbats=FALSE)


#####   
# scatter plot for diff genes
##### 


setwd("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/RNA-seq_polyA/all_replicates")

polyA <- read.csv("mj261_polyA_diff_osc_allreplicates.csv", header = TRUE)
polyAup <- subset(polyA, polyA$threshold == "A")
polyAdown <- subset(polyA, polyA$threshold == "B")

polyAnorm <- read.csv("mj261_polyA_normalized_allreplicates.csv", header = TRUE)

head(polyAnorm)

polyAnorm$meanwt <- rowMeans(polyAnorm[,2:4])
polyAnorm$meanmj261 <- rowMeans(polyAnorm[,5:7])

polyAupsig <- subset(polyAnorm, polyAnorm$X %in% polyAup$X)
polyAdownsig <- subset(polyAnorm, polyAnorm$X %in% polyAdown$X)
sensor <- subset(polyAnorm, polyAnorm$X == "piRNASensor")


reg = lm(polyAnorm$meanwt ~ polyAnorm$meanmj261)
summary(reg)

g <- ggplot(polyAnorm, aes(x=log2(polyAnorm$meanwt+1), y=log2(polyAnorm$meanmj261+1))) + geom_point(data=polyAupsig, aes(x=log2(polyAupsig$meanwt+1), y=log2(polyAupsig$meanmj261+1)), color= "Red", size=4) + geom_point(data=sensor, aes(x=log2(sensor$meanwt+1), y=log2(sensor$meanmj261+1)), color= "Darkgreen", size=4)   + geom_point(data=polyAdownsig, aes(x=log2(polyAdownsig$meanwt+1), y=log2(polyAdownsig$meanmj261+1)), color= "Blue", size=4) + geom_point(stat = "identity")  + theme_classic() + geom_abline(slope = 0.7527 , intercept = 0.016 ) + xlab("wild type log2(TPM+1)") + ylab("rpb-9(mj261)V log2(TPM+1)")

ggsave("mj261_polyA_scatter_allreplicates.pdf", g, height = 10, width = 10, useDingbats=FALSE)




#####   
# scatter plot for diff genes
##### 


setwd("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/RNA-seq_ribozero/")

polyA <- read.csv("mj261_ribozero_diff_osc.csv", header = TRUE)
polyAup <- subset(polyA, polyA$threshold == "A")
polyAdown <- subset(polyA, polyA$threshold == "B")

polyAnorm <- read.csv("mj261_ribozero_normalized.csv", header = TRUE)

head(polyAnorm)

polyAnorm$meanwt <- rowMeans(polyAnorm[,2:4])
polyAnorm$meanmj261 <- rowMeans(polyAnorm[,5:7])

polyAupsig <- subset(polyAnorm, polyAnorm$X %in% polyAup$X)
polyAdownsig <- subset(polyAnorm, polyAnorm$X %in% polyAdown$X)
sensor <- subset(polyAnorm, polyAnorm$X == "piRNASensor")


reg = lm(polyAnorm$meanwt ~ polyAnorm$meanmj261)
summary(reg)

g <- ggplot(polyAnorm, aes(x=log2(polyAnorm$meanwt+1), y=log2(polyAnorm$meanmj261+1))) + geom_point(data=polyAupsig, aes(x=log2(polyAupsig$meanwt+1), y=log2(polyAupsig$meanmj261+1)), color= "Red", size=4) + geom_point(data=sensor, aes(x=log2(sensor$meanwt+1), y=log2(sensor$meanmj261+1)), color= "Darkgreen", size=7.5)   + geom_point(data=polyAdownsig, aes(x=log2(polyAdownsig$meanwt+1), y=log2(polyAdownsig$meanmj261+1)), color= "Blue", size=4) + geom_point(stat = "identity")  + theme_classic() + geom_abline(slope = 0.9983 , intercept = 3.516e-04 ) + xlab("wild type log2(TPM+1)") + ylab("rpb-9(mj261)V log2(TPM+1)")

ggsave("mj261_ribo_scatter.pdf", g, height = 10, width = 10, useDingbats=FALSE)














