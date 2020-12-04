# this script was created by Ahmet Can Berkyurek on 21/03/2019 to check piRNA abundance in mj261 and rescue. 
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
setwd("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/lisa/rpb-p paper figures/piRNAs/miRNAsfiltered")

#introduce the countdata from featurecounts (fractional)
countdata <- read.table("mj261_5dep_21U.txt", header=TRUE, row.names = 1)
head(countdata)
countdata2 <- read.table("desp1_prg1.txt", header=TRUE, row.names = 1)
piRNAs <- read.delim("ce11_piRNAs.txt")
head(piRNAs)
# Remove first five columns (chr, start, end, strand, length)
countdata <- countdata[ ,6:ncol(countdata)]
head(countdata)

countdata2 <- countdata2[ ,6:ncol(countdata2)]
head(countdata2)

# Convert to matrix and remove all negative values and convert them to integers.

dds <- data.frame(countdata[,1:9])
colnames(dds) <- c("wt1", "wt2","wt3", "mutant1", "mutant2", "mutant3","rescue1", "rescue2", "rescue3")
head(dds)

dds <- data.frame(countdata[,1:9],countdata2$X3.SX1888.dep__5dependentAligned.sortedByCoord.out.bam)
colnames(dds) <- c("wt1", "wt2","wt3", "mutant1", "mutant2", "mutant3","rescue1", "rescue2", "rescue3","prg1")
head(dds)

keep <- rowSums(dds) >= 0
dds2 <- dds[keep,]
head(dds2)
colnames(dds2) <- c("wt1", "wt2","wt3", "mutant1", "mutant2", "mutant3","rescue1", "rescue2", "rescue3")
head(dds2)

ddsnormalized <- data.frame(dds2$wt1/(3681387), dds2$wt2/(7966827), dds2$wt3/(2621807),dds2$mutant1/(12959479),dds2$mutant2/(11826203),dds2$mutant3/(10736443),dds2$rescue1/(16037323),dds2$rescue2/(11245124),dds2$rescue3/(13053543),dds2$prg1/(18173039))
ddsnormalized <- ddsnormalized*1000000
head(ddsnormalized)
ddsnormalizedlog <- log(ddsnormalized+1,2)
head(ddsnormalizedlog)

colnames(ddsnormalizedlog) <- c("wt1", "wt2","wt3", "mutant1", "mutant2", "mutant3","rescue1", "rescue2", "rescue3","prg1")
rownames(ddsnormalizedlog) <- rownames(dds2)
ddsnormalizedlog$genes <- rownames(dds2)
head(ddsnormalizedlog)

ddsnormalizedlog2 <- subset(ddsnormalizedlog, ddsnormalizedlog$genes %in% piRNAs$Gene.stable.ID)
head(ddsnormalizedlog2)



ddsnormalizedlogmerged <- data.frame(rowMeans(ddsnormalizedlog2[,1:3]), rowMeans(ddsnormalizedlog2[,4:6]), rowMeans(ddsnormalizedlog2[,7:9]), ddsnormalizedlog2$prg1, ddsnormalizedlog2$genes)
colnames(ddsnormalizedlogmerged) <- c("wt" , "mutant", "rescue", "prg1", "genes")
rownames(ddsnormalizedlogmerged) <- rownames(ddsnormalizedlog2)
head(ddsnormalizedlogmerged)

ddsnormalizedlogmerged <- ddsnormalizedlogmerged[order(ddsnormalizedlogmerged$genes),]
head(ddsnormalizedlogmerged)




b <- c(rep("A",160), rep("B",160), rep("C",160), rep("D",160), rep("E",160), rep("F",160), rep("G",160),rep("H",160),rep("I",160),rep("J",160),rep("K",160),rep("L",160),rep("M",160),rep("N",160),rep("O",160),rep("P",160),rep("R",160),rep("S",160),rep("T",160),rep("U",199))
b <- data.frame(b)
head(b)

ddsnormalizedlogmerged$code <- b$b

a1 <- melt(ddsnormalizedlogmerged)
head(a1)


x <- (subset(ddsnormalizedlogmerged,ddsnormalizedlogmerged$mutant > 0))
y <- (subset(ddsnormalizedlogmerged,ddsnormalizedlogmerged$wt > 0))
z <- (subset(ddsnormalizedlogmerged,ddsnormalizedlogmerged$rescue > 0))
j <- (subset(ddsnormalizedlogmerged,ddsnormalizedlogmerged$prg1 > 0))


conplot <- ggplot(ddsnormalizedlogmerged, aes(x=code, y=wt)) + geom_boxplot(color = "Black", alpha=0.3, fatten=4) + theme_classic() 
conplot2 <- conplot + geom_boxplot(data=ddsnormalizedlogmerged, aes(x=code, y=mutant),color="Blue", alpha=0.1, fatten=4) 
conplot3 <- conplot2 + geom_boxplot(data=ddsnormalizedlogmerged, aes(x=code, y=rescue),color="Red", alpha=0.1, fatten=4)
conplot4 <- conplot3 + geom_boxplot(data=ddsnormalizedlogmerged, aes(x=code, y=prg1),color="Green", alpha=0.1, fatten=4)

plot <- ggplot(ddsnormalizedlogmerged, aes(x=code, y=wt)) + geom_point(color = "Black", stat = "identity") + theme_classic() 
plot2 <- plot + geom_point(data=ddsnormalizedlogmerged, aes(x=code, y=mutant),color="Blue", stat = "identity") + geom_point(data=ddsnormalizedlogmerged, aes(x=code, y=rescue),color="Red", stat = "identity") + geom_point(data=ddsnormalizedlogmerged, aes(x=code, y=prg1),color="Green", stat = "identity")

ggsave("piRNAs_percentile3.pdf", plot2, height = 15, width = 40, useDingbats=FALSE)

sig <- ks.test(ddsnormalizedlogmerged$wt, ddsnormalizedlogmerged$mutant, alternative = c("two.sided"), exact = NULL)


# subset chapaev transposon piRNAs, all replicates or mean!
setwd("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/piRNAs")

data <- read.csv("rpb9(mj261)vswt_allpiRNAs.csv", header = TRUE)
data <- data[order(data$genes),]
data$prg1 <- ddsnormalizedlogmerged$prg1
head(data)

chapaev <- read.delim("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/piRNAs/Capaev_piRNAs.txt", header = TRUE)

cemudr <-  read.delim("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/piRNAs/Cemudr_piRNAs.txt", header = TRUE)

cahapev1 <- subset(ddsnormalizedlog2, ddsnormalizedlog2$genes %in% chapaev$Gene.stable.ID)
WBGene00174441 <- subset(cahapev1, cahapev1$genes == "WBGene00174441")
WBGene00172213 <-  subset(cahapev1, cahapev1$genes == "WBGene00172213")

wt <- c(0.6259938 ,0.5869635,   0)
mut <- c(0, 0   ,   0)
rescue <- c(0.6991478, 0.8482749 ,0.2984503)
prg1 <- c(0,0,0)

WBGene00172213_all <- data.frame(wt, mut, rescue, prg1)


cemudr1 <- subset(ddsnormalizedlog2, ddsnormalizedlog2$genes %in% cemudr$Gene.stable.ID)
WBGene00045722 <- subset(cemudr1, cemudr1$genes == "WBGene00045722")

wt <- c(0.3466866 ,0.461069, 0.4661482)
mut <- c(0.1072374, 0.1171072    ,   0)
rescue <- c(0.642647, 0.4390483, 0.5456179)
prg1 <- c(0,0,0)

WBGene00045722_all <- data.frame(wt, mut, rescue, prg1)

a2 <- melt(WBGene00172213_all)
head(a2)


plot <- ggplot(a2, aes(x=variable, y=value)) + geom_boxplot() + scale_x_discrete(limits = c("wt", "mut","rescue")) + theme_classic() + ylab("log2(reads per million +1)")

ggsave("WBGene00172213_chapaev_allreplicates.pdf", plot, height = 5, width = 5, useDingbats=FALSE)


# subset chapaev transposon piRNAs, all replicates or mean!

setwd("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/piRNAs/")


keep <- rowSums(dds) >= 0
dds2 <- dds[keep,]
head(dds2)
colnames(dds2) <- c("wt1", "wt2","wt3", "mutant1", "mutant2", "mutant3","rescue1", "rescue2", "rescue3")
head(dds2)

factor <- (3681387+7966827+2621807+12959479+11826203+10736443+16037323+11245124+13053543)/9


dds2$wt1 <- dds2$wt1/(3681387/factor)

dds2$wt2 <- dds2$wt2/(7966827/factor)

dds2$wt3 <- dds$wt3/(2621807/factor)

dds2$mutant1 <- dds2$mutant1/(12959479/factor)

dds2$mutant2 <- dds2$mutant2/(11826203/factor)

dds2$mutant3 <- dds2$mutant3/(10736443/factor)

dds2$rescue1 <- dds2$rescue1/(16037323/factor)

dds2$rescue2 <- dds2$rescue2/(11245124/factor)

dds2$rescue3 <- dds2$rescue3/(13053543/factor)

head(dds2)


WBGene00172213 <-  subset(dds2, rownames(dds2) == "WBGene00172213")

wt <- c(5.44047 , 5.027968,   0)
mut <- c(0, 0   ,   0)
rescue <- c(6.244332, 8.014864, 2.301499)

WBGene00172213_all <- data.frame(wt, mut, rescue)


a2 <- melt(WBGene00172213_all)
head(a2)


plot <- ggplot(a2, aes(x=variable, y=value)) + geom_boxplot() + scale_x_discrete(limits = c("wt", "mut","rescue")) + theme_classic() + ylab("Normalized Reads") + scale_y_continuous(limits=c(0,10), breaks=c(1,2,3,4,5,6,7,8,9,10))

ggsave("WBGene00172213_chapaev_allreplicates3.pdf", plot, height = 5, width = 5, useDingbats=FALSE)




WBGene00045722 <-  subset(dds2, rownames(dds2) == "WBGene00045722")

wt <- c(2.720235 ,3.770976 ,3.819594)
mut <- c(0.7727346, 0.7170428  ,     0)
rescue <- c(5.619899 ,3.562162, 4.602997)

WBGene00045722_all <- data.frame(wt, mut, rescue)


a2 <- melt(WBGene00045722_all)
head(a2)


plot <- ggplot(a2, aes(x=variable, y=value)) + geom_boxplot() + scale_x_discrete(limits = c("wt", "mut","rescue")) + theme_classic() + ylab("Normalized Reads") + scale_y_continuous(limits=c(0,10), breaks=c(1,2,3,4,5,6,7,8,9,10))

ggsave("WBGene00045722_cemudr_allreplicates3.pdf", plot, height = 5, width = 5, useDingbats=FALSE)


########
# repeat the analysis on cluster I and cliuster II of chrIV, nd maybe motif-independent piRNAs!
######
#load featureCounts count.txt file. before loading, delete the first two rows for clarity.!!!textwrangler should help...
setwd("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/piRNAs/clusters")

#introduce the countdata from featurecounts (fractional)
countdata <- read.table("mj261_5dep_21U.txt", header=TRUE, row.names = 1)
head(countdata)

countdata <- read.table("./new/mj261_5dep_new.txt", header=TRUE, row.names = 1)
head(countdata)

# Remove first five columns (chr, start, end, strand, length)
countdata <- countdata[ ,6:ncol(countdata)]
head(countdata)



# Convert to matrix and remove all negative values and convert them to integers.

dds <- data.frame(countdata[,1:9])
colnames(dds) <- c("wt1", "wt2","wt3", "mutant1", "mutant2", "mutant3","rescue1", "rescue2", "rescue3")
head(dds)

keep <- rowSums(dds) >= 0
dds2 <- dds[keep,]
head(dds2)
colnames(dds2) <- c("wt1", "wt2","wt3", "mutant1", "mutant2", "mutant3","rescue1", "rescue2", "rescue3")
head(dds2)

ddsnormalized <- data.frame(dds2$wt1/(3681387), dds2$wt2/(7966827), dds2$wt3/(2621807),dds2$mutant1/(12959479),dds2$mutant2/(11826203),dds2$mutant3/(10736443),dds2$rescue1/(16037323),dds2$rescue2/(11245124),dds2$rescue3/(13053543))
ddsnormalized <- ddsnormalized*1000000
head(ddsnormalized)
ddsnormalizedlog <- log10(ddsnormalized+1)
head(ddsnormalizedlog)

ddsnormalized <- data.frame(dds2$wt1/(2634941), dds2$wt3/(1813824),dds2$mutant1/(9000499),dds2$mutant2/(8124378),dds2$mutant3/(7396277),dds2$rescue1/(10924014),dds2$rescue2/(7661417),dds2$rescue3/(3782503))
ddsnormalized <- ddsnormalized*1000000
head(ddsnormalized)
ddsnormalizedlog <- log10(ddsnormalized+1)
head(ddsnormalizedlog)


colnames(ddsnormalizedlog) <- c("wt1","wt3", "mutant1", "mutant2", "mutant3","rescue1", "rescue2", "rescue3")
rownames(ddsnormalizedlog) <- rownames(dds2)
ddsnormalizedlog$genes <- rownames(dds2)
head(ddsnormalizedlog)

piRNAs <- read.delim("piRNAs.txt", header = F)


ddsnormalizedlog2 <- subset(ddsnormalizedlog, ddsnormalizedlog$genes %in% piRNAs$V1)
head(ddsnormalizedlog2)

ddsnormalizedlog2$meanwt <- rowMeans(ddsnormalizedlog2[,1:2])

ddsnormalizedlog2$meanmut <- rowMeans(ddsnormalizedlog2[,3:5])

ddsnormalizedlog2$meanrescue <- rowMeans(ddsnormalizedlog2[,6:8])

head(ddsnormalizedlog2)

ddsnormalizedlog2 <- ddsnormalizedlog2[rowSums(ddsnormalizedlog2[,1:8]) > 0,]


a <- melt(ddsnormalizedlog2)
head(a)

p <- ggplot(a, aes(x=variable, y=value)) + geom_boxplot() + scale_x_discrete(limits=c("meanwt", "meanmut", "meanrescue")) + theme_classic()

ggsave("./new/piRNAs_all_new.pdf", p, height = 5, width = 5, useDingbats=FALSE)




########
#introduce clusterI and cluster II piRNAs, repeat boxplots
########

coord <- read.delim("piRNAlocation.bed", header = F)
head(coord) 

clusterI <- subset(coord, coord$V2 > 4000000 & coord$V2 < 7500000)
write.csv(clusterI, "piRNAs_clusterI.csv")
clusterInames <- read.delim("clusterI_names.tsv", header = F)

clusterII <-  subset(coord, coord$V2 > 13000000 & coord$V2 < 17500000)
write.csv(clusterII, "piRNAs_clusterII.csv")
clusterIInames <- read.delim("clusterII_names.tsv", header = F)


clusterIlist <- subset(ddsnormalizedlog2, ddsnormalizedlog2$genes %in% clusterInames$V1)
b <- melt(clusterIlist)
head(b)

p <- ggplot(b, aes(x=variable, y=value)) + geom_boxplot() + scale_x_discrete(limits=c("meanwt", "meanmut", "meanrescue")) + theme_classic()

ggsave("./new/piRNAs_clusterI_new.pdf", p, height = 5, width = 5, useDingbats=FALSE)


clusterIIlist <- subset(ddsnormalizedlog2, ddsnormalizedlog2$genes %in% clusterIInames$V1)
c <- melt(clusterIIlist)
head(c)

p <- ggplot(c, aes(x=variable, y=value)) + geom_boxplot() + scale_x_discrete(limits=c("meanwt", "meanmut", "meanrescue")) + theme_classic()

ggsave("./new/piRNAs_clusterII_new.pdf", p, height = 5, width = 5, useDingbats=FALSE)


other <- subset(coord, coord$V1 != "chrIV")
write.csv(other, "./new/motif_independent_piRNAs.csv")

motifind <- read.delim("./new/motif_independent.tsv", header = F)

motifindlist <- subset(ddsnormalizedlog2, ddsnormalizedlog2$genes %in% motifind$V1)
d <- melt(motifindlist)
head(d)

p <- ggplot(d, aes(x=variable, y=value)) + geom_boxplot() + scale_x_discrete(limits=c("meanwt", "meanmut", "meanrescue")) + theme_classic()

ggsave("./new/piRNAs_motif_indep_new.pdf", p, height = 5, width = 5, useDingbats=FALSE)





b <- c(rep("A",160), rep("B",160), rep("C",160), rep("D",160), rep("E",160), rep("F",160), rep("G",160),rep("H",160),rep("I",160),rep("J",160),rep("K",160),rep("L",160),rep("M",160),rep("N",160),rep("O",160),rep("P",160),rep("R",160),rep("S",160),rep("T",160),rep("U",199))
b <- data.frame(b)
head(b)

ddsnormalizedlogmerged$code <- b$b

a1 <- melt(ddsnormalizedlogmerged)
head(a1)


x <- (subset(ddsnormalizedlogmerged,ddsnormalizedlogmerged$mutant > 0))
y <- (subset(ddsnormalizedlogmerged,ddsnormalizedlogmerged$wt > 0))
z <- (subset(ddsnormalizedlogmerged,ddsnormalizedlogmerged$rescue > 0))
j <- (subset(ddsnormalizedlogmerged,ddsnormalizedlogmerged$prg1 > 0))


conplot <- ggplot(ddsnormalizedlogmerged, aes(x=code, y=wt)) + geom_boxplot(color = "Black", alpha=0.3, fatten=4) + theme_classic() 
conplot2 <- conplot + geom_boxplot(data=ddsnormalizedlogmerged, aes(x=code, y=mutant),color="Blue", alpha=0.1, fatten=4) 
conplot3 <- conplot2 + geom_boxplot(data=ddsnormalizedlogmerged, aes(x=code, y=rescue),color="Red", alpha=0.1, fatten=4)
conplot4 <- conplot3 + geom_boxplot(data=ddsnormalizedlogmerged, aes(x=code, y=prg1),color="Green", alpha=0.1, fatten=4)

plot <- ggplot(ddsnormalizedlogmerged, aes(x=code, y=wt)) + geom_point(color = "Black", stat = "identity") + theme_classic() 
plot2 <- plot + geom_point(data=ddsnormalizedlogmerged, aes(x=code, y=mutant),color="Blue", stat = "identity") + geom_point(data=ddsnormalizedlogmerged, aes(x=code, y=rescue),color="Red", stat = "identity") + geom_point(data=ddsnormalizedlogmerged, aes(x=code, y=prg1),color="Green", stat = "identity")

ggsave("piRNAs_percentile3.pdf", plot2, height = 15, width = 40, useDingbats=FALSE)

sig <- ks.test(ddsnormalizedlogmerged$wt, ddsnormalizedlogmerged$mutant, alternative = c("two.sided"), exact = NULL)





