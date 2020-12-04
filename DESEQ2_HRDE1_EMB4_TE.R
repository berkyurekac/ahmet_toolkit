# this script was created by Ahmet Can Berkyurek on 03/05/2019, to analyze differential TE expression analysis with RNA-seq data for hrde-1 and emb-4. 
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


# load featureCounts count.txt file. before loading, delete the first two rows for clarity.!!!textwrangler should help...
setwd("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/cambridge-UK/NGS/RNAseq/alper/TE")

wt1 <- read.delim("ERR1530686Aligned.TEs_wt_expression.txt", header = FALSE, stringsAsFactors = FALSE)
wt2 <- read.delim("ERR1530688Aligned.TEs_wt_expression.txt", header = FALSE, stringsAsFactors = FALSE)
wt3 <- read.delim("ERR1530690Aligned.TEs_wt_expression.txt", header = FALSE, stringsAsFactors = FALSE)
hrde1 <- read.delim("ERR1530693Aligned.TEs_wt_expression.txt", header = FALSE, stringsAsFactors = FALSE)
hrde2 <- read.delim("ERR1530695Aligned.TEs_wt_expression.txt", header = FALSE, stringsAsFactors = FALSE)
hrde3 <- read.delim("ERR1530697Aligned.TEs_wt_expression.txt", header = FALSE, stringsAsFactors = FALSE)
emb1 <- read.delim("ERR1530700Aligned.TEs_wt_expression.txt", header = FALSE, stringsAsFactors = FALSE)
emb2 <- read.delim("ERR1530702Aligned.TEs_wt_expression.txt", header = FALSE, stringsAsFactors = FALSE)
emb3 <- read.delim("ERR1530704Aligned.TEs_wt_expression.txt", header = FALSE, stringsAsFactors = FALSE)


frame <- data.frame(wt1[,c(4,7)],wt2$V7, wt3$V7, hrde1$V7, hrde2$V7, hrde3$V7, emb1$V7, emb2$V7, emb3$V7)
head(frame)
colnames(frame) <- c("genes","wt-1", "wt-2", "wt-3","hrde1-1", "hrde1-2", "hrde1-3", "emb4-1", "emb4-2", "emb4-3")
rownames(frame) <- frame$genes
head(frame)

# Remove first five columns (chr, start, end, strand, length)

countdata_hrde1 <- frame[,c(2,3,4,5,6,7)]
head(countdata_hrde1)

#keep <- rowSums((countdata)) >= 10
#dds <- countdata[keep,]
#head(dds)

# Convert to matrix and remove all negative values and convert them to integers.

countdata_hrde1 <- round(countdata_hrde1)
head(countdata_hrde1)
(condition <- factor(c(rep("wt", 3), rep("hrde1", 3))))
(coldata <- data.frame(row.names=colnames(countdata_hrde1), condition))

# convert to Deseq2 matrix
dds <- DESeqDataSetFromMatrix(countData=countdata_hrde1, colData=coldata, design=~condition)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds

# keep the comparison wt. 
dds$condition <- relevel(dds$condition, ref = "wt")
dds <- DESeq(dds, betaPrior = FALSE)


## print the results with log2fc and padj
res <- results(dds)

## Order by adjusted p-value
res <- res[order(-res$log2FoldChange),]
head(res)
write.csv(res, "mj609_diff_TE.csv")


## Merge with normalized count data
plotMA(res,ylim=c(-10,10),main="Ribozero Total RNA")


resdata = as.data.frame(dplyr::mutate(as.data.frame(res), sig=ifelse(res$padj<0.05, "FDR<0.05", "Not Sig")), row.names=rownames(res))
resdata <- resdata %>% mutate(threshold = ifelse(log2FoldChange >= 1 & padj < 0.05,"A", ifelse(log2FoldChange< -1 & padj < 0.05 , "B", "C")))
resdata <- resdata[order(-resdata$log2FoldChange),]
rownames(resdata) <- rownames(res)
resdata$genes <- rownames(resdata)
head(resdata)
write.csv(resdata, "hrde1_diff_TE.csv")

osc <- read.csv("/Volumes/miska/Ahmet Berkyurek/sequencing/NGS/RNA-seq/Lisa_ribozero/mapped/all_genes/oscilatinggenes.tsv")
osc <- data.frame(osc)
head(osc)

diff_ext <- subset(resdata, !(rownames(resdata) %in% osc$WBGene00000002))
head(diff_ext)
write.csv(diff_ext, "emb4_diff_osc.csv")

intestine <- read.delim("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/cambridge-UK/data/PARP-1/ChIP-seq/tissue/intestine_genes_names.bed", header = FALSE)
head(intestine)

int <- subset(diff_ext, rownames(diff_ext) %in% intestine$V1 )
head(int)

hrde1_subset <- subset(diff_ext, rownames(resdata) %in% hrde1_frame$V1)
head(hrde1_subset)
write.csv(hrde1_subset, "ribozero_Lisa_hrde1_parp1.csv")

gonad_subset <- subset(diff_ext, rownames(diff_ext) %in% gonad$V1)
head(gonad_subset)

a <- c("piRNASensor")
a <- as.data.frame(a)
colnames(a) <- c("name")
a

b <- subset(diff_ext, rownames(diff_ext) %in% a$name)
b

c <- subset(diff_ext, rownames(diff_ext) %in% hrde1_frame$V1)
c

d <-  subset(diff_ext, rownames(diff_ext) %in% gonad$V1)
d


#volcanoplot
p <- ggplot(resdata, aes(x=log2FoldChange, y=-log10(padj))) +
  geom_point(aes(colour = threshold), size=3) +
  geom_point(data = b, color = "Darkgreen", size=3) + 
  geom_point(data = d, color = "Purple", size = 2) +
  geom_point(data = c, color = "Yellow", size=1) + scale_colour_manual(values = c("A"= "red", "B"="blue",  "C"= "black")) + scale_x_continuous(limits = c(-8,8)) + theme_classic()

p <- ggplot(diff_ext, aes(x=log2FoldChange, y=-log10(padj))) +
  geom_point(aes(colour = threshold), size=3) +
  geom_point(data = int, color = "Darkgreen", size=1) + scale_colour_manual(values = c("A"= "red", "B"="blue",  "C"= "black")) + scale_x_continuous(limits = c(-8,8)) + theme_classic()

ggsave("mj613_intestine_RNAseq.pdf", p, height = 10, width = 10, useDingbats=FALSE)

# correlation between replicates!
## Examine plot of p-values)
hist(res$pvalue, breaks=50, col="grey")

# PCA of the replicates!
(mycols <- brewer.pal(8, "Dark2")[1:length(unique(condition))])


rld_pca <- function (rld, intgroup = "condition", ntop = 500, colors=NULL, legendpos="bottomleft", main="PCA Biplot", textcx=1, ...) {
  require(genefilter)
  require(calibrate)
  require(RColorBrewer)
  rv = rowVars(assay(rld))
  select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(rld)[select, ]))
  fac = factor(apply(as.data.frame(colData(rld)[, intgroup, drop = FALSE]), 1, paste, collapse = " : "))
  if (is.null(colors)) {
    if (nlevels(fac) >= 3) {
      colors = brewer.pal(nlevels(fac), "Paired")
    }   else {
      colors = c("black", "red")
    }
  }
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab <- paste0("PC1 (",as.character(pc2var),"%)")
  plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], pch=21, xlab=pc1lab, ylab=pc2lab, main=main, ...)
  with(as.data.frame(pca$x), textxy(PC1, PC2, labs=rownames(as.data.frame(pca$x)), cex=textcx))
  legend(legendpos, legend=levels(fac), col=colors, pch=20)
  #     rldyplot(PC2 ~ PC1, groups = fac, data = as.data.frame(pca$rld),
  #            pch = 16, cerld = 2, aspect = "iso", col = colours, main = draw.key(key = list(rect = list(col = colours),
  #                                                                                         terldt = list(levels(fac)), rep = FALSE)))
}
png("qc-pca.png", 1000, 1000, pointsize=20)
rld_pca(rld, colors=mycols, intgroup="condition", xlim=c(-50, 50),ylim=c(-50,50))
dev.off()

# determine library size factor: 
# I will run the manual analysis on 25 genes from lisa's AMA-1 ChIP-data. 
## wt replicates


wt1length <-  56111530
head(wt1length)



wt2length <-  89769958
head(wt2length)


wt3length <- 75693370
head(wt3length)

## mj261 replicates

mj261_1length <- 79040357 
head(mj261_1length)



mj261_2length <- 69150324
head(mj261_2length)


mj261_3length <- 39509223
head(mj261_3length)

## determine individual factors by dividing with the average

averagelibrarysize <- 68212460.3
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

countdataframenrm <- data.frame(dds$wt1/wt1factor,dds$wt2/wt2factor,dds$wt3/wt3factor,dds$mutant1/mj261_1factor,dds$mutant2/mj261_2factor,dds$mutant3/mj261_3factor)
head(countdataframenrm)
rownames(countdataframenrm) = rownames(dds)
colnames(countdataframenrm) = c("wt1", "wt2", "wt3", "mutant1", "mutant2", "mutant3")
head(countdataframenrm)

## divide each column with average wt counts


countdataframenrm <- countdataframenrm[1:6]+0.1
countdataframenrm$meanwt <- rowMeans(countdataframenrm[1:3])
head(countdataframenrm)

countdataframenrmfinal <- countdataframenrm[1:6]/countdataframenrm$meanwt
head(countdataframenrmfinal)
countdataframenrmfinal <- countdataframenrmfinal[order(-countdataframenrmfinal$mutant2),]

countdataframenrmfinal$genes <- rownames(countdataframenrmfinal)                                                  
head(countdataframenrmfinal)

countdataframenrmfinal <- log(countdataframenrmfinal[1:6], 2)
countdataframenrmfinal$genes <- rownames(countdataframenrmfinal)
colnames(countdataframenrmfinal) <- c("wt1", "wt2", "wt3", "mutant1", "mutant2", "mutant3","genes")
head(countdataframenrmfinal)


#subset with 99 gene list from Lisa's ChIP-seq
upregulated <- read.csv("/Volumes/miska/Ahmet Berkyurek/sequencing/NGS/RNA-seq/RBP9/5independent/combined/allsmallRNAs/hrde1targets.tsv")
upregulatedgenes <- data.frame(upregulated$WBGene00000089)
upregulatedgenes2 <- data.frame(upregulatedgenes)
head(upregulatedgenes2)

hrde1 <- read.csv("/Volumes/miska/Ahmet Berkyurek/sequencing/NGS/RNA-seq/Lisa_ribozero/mapped/all_genes/99genelist.csv")
upregulatedgenes <- data.frame(hrde1$X)
hrde1genes <- data.frame(upregulatedgenes)
head(hrde1genes)

upregulated_subset <- subset(countdataframenrmfinal, rownames(countdataframenrmfinal) %in% genes58$`genes68fromtotalRNaseqsubset$genes`)
head(upregulated_subset)
upregulated_subset <- upregulated_subset[order(upregulated_subset[,'genes']), ]

write.csv(upregulated_subset, "99genelist.csv")

upregulated_subset <- upregulated_subset[order(-upregulated_subset$mutant2),]



upregulatedsubsetl2fc1 <- subset(upregulated_subset, mutant2 > 0)
head(upregulatedsubsetl2fc1)
write.csv(upregulatedsubsetl2fc1, "68gebes_Lisa.csv")
genesupregulated68 <- as.data.frame(upregulatedsubsetl2fc1$genes)
head(genesupregulated68)

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 299)
col_breaks = c(seq(-5,-2.0,length=100),  # for red
               seq(-1.99,1.99,length=100),           # for white
               seq(2.0,5,length=100))

heatmap.2(as.matrix(upregulated_subset[1:6]),      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          col=my_palette,       # use on color palette defined earlier
          breaks=col_breaks,    # enable color transition at specified limits
          Colv = FALSE,
          Rowv=FALSE,
          dendrogram="none")     # only draw a row dendrogram 



# subset with top25 genes from the ChIP-seq data
top25 <- read.csv("top25.csv")
top25frame <- as.data.frame(top25)
rownames(top25frame) <- top25frame$x
rownames(top25frame)
top25frame_subset <- subset(countdataframenrmfinal, rownames(countdataframenrmfinal) %in% top25frame$x)

melted <- top25frame_subset %>% gather(sample, FE, -genes)
head(melted)

row_order <- order(melted$FE, decreasing = TRUE)

melted <- melted[row_order, , drop = FALSE]
head(melted)

# subset with 99 genes from Lisa's AMA-1 ChIP-seq data. 
upregulated <- read.csv("/Volumes/miska/Ahmet Berkyurek/sequencing/NGS/RNA-seq/RBP9/5independent/combined/allsmallRNAs/99genelist.tsv")
upregulatedgenes <- data.frame(upregulated$WBGene00000089)
upregulatedgenes2 <- data.frame(upregulatedgenes)
head(upregulatedgenes2)

upregulatednew <- read.xls("/Volumes/miska/Ahmet Berkyurek/sequencing/NGS/RNA-seq/Lisa_ribozero/mapped/all_genes/mj261_upregulated.xlsx")
upregulatedgenesnew <- data.frame(upregulatednew$WBGene00008010)
upregulatedgenes2new <- data.frame(upregulatedgenesnew)
head(upregulatedgenes2new)


upregulated_subset <- subset(countdataframenrmfinal, rownames(countdataframenrmfinal) %in% upregulatedgenes2$upregulated.WBGene00000089)
head(upregulated_subset)

upregulated_subset <- upregulated_subset[order(-upregulated_subset$mutant2),]
list1 <- write.csv(upregulated_subset[1:132,], "list1.csv")
list2 <- write.csv(upregulated_subset[273:443,], "list2.csv")
list3 <- write.csv(upregulated_subset, "list3.csv")

## divide each column with average wt counts


countdataframenrm <- countdataframenrm[1:9]+0.1
countdataframenrm$meanwt <- rowMeans(countdataframenrm[1:3])
head(countdataframenrm)

countdataframenrmfinal <- countdataframenrm[1:9]/countdataframenrm$meanwt
head(countdataframenrmfinal)

countdataframenrmfinal <- log(countdataframenrmfinal[1:9], 2)
countdataframenrmfinal$genes <- rownames(countdataframenrmfinal)
colnames(countdataframenrmfinal) <- c("wt1", "wt2", "wt3", "mutant1", "mutant2", "mutant3", "rescue1", "rescue2", "rescue3", 'genes')
head(countdataframenrmfinal)


my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 299)
col_breaks = c(seq(-5,-2.0,length=100),  # for red
               seq(-1.99,1.99,length=100),           # for white
               seq(2.0,5,length=100))

heatmap.2(as.matrix(upregulated_subset[1:9]),      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          col=my_palette,       # use on color palette defined earlier
          breaks=col_breaks,    # enable color transition at specified limits
          Colv = FALSE,
          dendrogram="row")     # only draw a row dendrogram 


# Note for Ahmet: reshaoe the hetmaps! p <- p + facet_wrap(~rowMeanClass +experiment , scales = "free", ncol = 2)




# use ggplot2, this brings better vector graphics. geom_tile function is for the heatmap!

p <- ggplot(melted, aes(x = sample, y = genes, fill = FE)) + geom_tile(color = "gray") + scale_fill_gradientn(colours = c("white","grey","black","red","red4"), 
                                                                                                              values = rescale(c(0,2.5,5,7.5,10)), guide = "colourbar") + 
  theme(axis.text.y = element_text(color="black", 
                                   
                                   size=6))


ggsave("snRNAs_heatmap_newnames_all.pdf",p, width = 6, height = 12)


# select the first 25 rows!

top <- countdataframenrmfinal[1:25,]
topmelted <- top %>% gather(sample2, FE2,-genes)

p <- ggplot(topmelted, aes(x = sample2, y = genes, fill = FE2)) + geom_tile(color = "gray") + scale_fill_gradientn(colours = c("white","grey","black","red","red4"), 
                                                                                                                   values = rescale(c(0,2.5,5,7.5,10)), guide = "colourbar") + 
  theme(axis.text.y = element_text(color="black", 
                                   
                                   size=12))


ggsave("snRNAs_heatmap_newnames_all.pdf",p, width = 6, height = 12)

# filter the differentially expressed genes ith osciallting genes! 
# intorduce the osciallting genes
diff <- read.csv("/Volumes/miska/Ahmet Berkyurek/sequencing/NGS/RNA-seq/Lisa_ribozero/mapped/all_genes/lisa_ribozero_allgenes_padj.csv")
diff <- data.frame(diff)
rownames(diff) <- diff$X
head(diff)

# introduce the oscilating genes, then subset!

osc <- read.csv("/Volumes/miska/Ahmet Berkyurek/sequencing/NGS/RNA-seq/Lisa_ribozero/mapped/all_genes/oscilatinggenes.tsv")
osc <- data.frame(osc)
head(osc)

# subset the data frame with oscilating genes
diff_ext <- subset(diff, !(rownames(diff) %in% osc$WBGene00000002))
head(diff_ext)
write.csv(diff_ext, "diff_oscillating.csv")

all <- read.csv("lisa_ribozero_allgenes_padj_oscillatingfiltered.csv")
all <- data.frame(all)
rownames(all) <- all$X.1
head(all)

germlinesubset <- subset(all, rownames(all) %in% germline_frame$V1)
head(germlinesubset)
write.csv(germlinesubset, "ribozero_Lisa_germline.csv")

somasubset <- subset(all, rownames(all) %in% soma_frame$V1)
head(somasubset)
write.csv(somasubset, "ribozero_Lisa_soma.csv")

hrde1subset <- subset(all, rownames(all) %in% hrde1_frame$V1)
head(hrde1subset)
write.csv(hrde1subset, "ribozero_Lisa_hrde1.csv")



