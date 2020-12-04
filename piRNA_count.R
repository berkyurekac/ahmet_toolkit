# this script was created by Ahmet Can Berkyurek on 10/10/2018, to analyze piRNAs in small RNA-seq data.
# I want to analyze mj261 5` dependent libraries for piRNAs on chrIV and the ones on other chromosomes. 

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
setwd("/Volumes/miska/Ahmet Berkyurek/sequencing/NGS/RNA-seq/RBP9/5dependent")

#introduce miRNA genes filtered from ENSEMBL 92 gtf file. this has been prepared with grep unix comand on bash. 
piRNAs <- read.delim("Caenorhabditis_elegans.WBcel235.piRNAs.bed", header = FALSE)
piRNAs <- data.frame(piRNAs)
rownames(piRNAs) <- piRNAs$V4
head(piRNAs)
piRNAlist <- data.frame(piRNAs$V4)
rownames(piRNAlist) <- rownames(piRNAs)
head(piRNAlist)

#subset for only chrIV -- > this will give the list for the piRNAs on chrIV only!
piRNAlistchrIV <- subset(piRNAs, piRNAs$V1 == "chrIV")
head(piRNAlistchrIV)
piRNAlistchrIVonly <- data.frame(rownames(piRNAlistchrIV))
head(piRNAlistchrIVonly)
rownames(piRNAlistchrIVonly) <- piRNAlistchrIVonly$rownames.piRNAlistchrIV.
colnames(piRNAlistchrIVonly) <- "genes"
head(piRNAlistchrIVonly)

#subset for those that don't have the chrIV
piRNAnonchrIV <- subset(piRNAs, piRNAs$V1 != "chrIV")
head(piRNAnonchrIV)
piRNAlistnonchrIVonly <- data.frame(rownames(piRNAnonchrIV))
rownames(piRNAlistnonchrIVonly) <- piRNAlistnonchrIVonly$rownames.piRNAnonchrIV.
colnames(piRNAlistnonchrIVonly) <- "genes"
head(piRNAlistnonchrIVonly)

#introduce the countdata from featurecounts (fractional) for input and IP.
input <- read.table("mj261_5dep_21U.txt", header=TRUE, row.names = 1)
input <- input[ ,6:ncol(input)]
input <- data.frame(input)
head(input)
colnames(input) <- c("wt1", "wt2","wt3","mj261_1", "mj261_2", "mj261_3", "mj261_rescue1", "mj261_rescue2", "mj261_rescue3" )
input2 <- data.frame(input$wt1, input$wt3, input$mj261_1, input$mj261_2, input$mj261_3, input$mj261_rescue1, input$mj261_rescue2, input$mj261_rescue3)
colnames(input2) <- c("wt1","wt3","mj261_1", "mj261_2", "mj261_3", "mj261_rescue1", "mj261_rescue2", "mj261_rescue3" )
rownames(input2) <- rownames(input) 
head(input2)

dds<- subset(input2, rownames(input) %in% rownames(piRNAlist))
head(dds)

# normalize to library size factor
wt1 <- 2634941/1000000
wt3 <- 1813824/1000000
mj261_1 <- 9000499/1000000	
mj261_2 <- 8124378/1000000	
mj261_3 <- 7396277/1000000
rescue1 <- 10924014/1000000
rescue2 <- 7661417/1000000
rescue3 <- 3782503/1000000


countdataframenrm <- data.frame(dds$wt1/wt1, dds$wt3/wt3, dds$mj261_1/mj261_1, dds$mj261_2/mj261_2, dds$mj261_3/mj261_3, dds$mj261_rescue1/rescue1, dds$mj261_rescue2/rescue2, dds$mj261_rescue3/rescue3)
head(countdataframenrm)
rownames(countdataframenrm) = rownames(dds)
colnames(countdataframenrm) = c("wt1","wt3","mj261_1", "mj261_2", "mj261_3", "mj261_rescue1", "mj261_rescue2", "mj261_rescue3" )
head(countdataframenrm)


countdataframenrm <- subset(countdataframenrm, rowSums(countdataframenrm) >= 5)

countdataframenrm <- log(countdataframenrm[1:8]+1,2)
countdataframenrm$genes <- rownames(countdataframenrm)
head(countdataframenrm)
tail(countdataframenrm)

melted <- melt(countdataframenrm, id.vars = "genes")

g <- ggplot(melted, aes(variable,value)) + geom_boxplot() + theme(text = element_text(size = 12)) + labs(x = "Experiment", y = "log2(read+1/million)") + theme_classic() + stat_summary(colour = "Darkblue",fun.y=mean, geom="point", shape=23, size=4)
ggsave("allpiRNAs_allreplicates.pdf", g, height = 10, width = 10)

# repeat the same for only chrIV
dds2<- subset(input2, rownames(input) %in% rownames(piRNAlistchrIVonly))
head(dds2)

countdataframenrm <- data.frame(dds2$wt1/wt1, dds2$wt3/wt3, dds2$mj261_1/mj261_1, dds2$mj261_2/mj261_2, dds2$mj261_3/mj261_3, dds2$mj261_rescue1/rescue1, dds2$mj261_rescue2/rescue2, dds2$mj261_rescue3/rescue3)
head(countdataframenrm)
rownames(countdataframenrm) = rownames(dds2)
colnames(countdataframenrm) = c("wt1","wt3","mj261_1", "mj261_2", "mj261_3", "mj261_rescue1", "mj261_rescue2", "mj261_rescue3" )
head(countdataframenrm)

countdataframenrm <- subset(countdataframenrm, rowSums(countdataframenrm) >= 5)


countdataframenrm <- log(countdataframenrm[1:8]+1,2)
countdataframenrm$genes <- rownames(countdataframenrm)
head(countdataframenrm)
tail(countdataframenrm)

melted <- melt(countdataframenrm, id.vars = "genes")

g <- ggplot(melted, aes(variable,value)) + geom_boxplot() + theme(text = element_text(size = 12)) + labs(x = "Experiment", y = "log2(read+1/million)") + theme_classic() + stat_summary(colour = "Darkblue",fun.y=mean, geom="point", shape=23, size=4)
ggsave("chrIVpiRNAs_allreplicates.pdf", g, height = 10, width = 10)

# repeat the same for piRNAs outisde chrIV
dds3 <- subset(input2, rownames(input) %in% rownames(piRNAlistnonchrIVonly))
head(dds3)

countdataframenrm <- data.frame(dds3$wt1/wt1, dds3$wt3/wt3, dds3$mj261_1/mj261_1, dds3$mj261_2/mj261_2, dds3$mj261_3/mj261_3, dds3$mj261_rescue1/rescue1, dds3$mj261_rescue2/rescue2, dds3$mj261_rescue3/rescue3)
head(countdataframenrm)
rownames(countdataframenrm) = rownames(dds3)
colnames(countdataframenrm) = c("wt1","wt3","mj261_1", "mj261_2", "mj261_3", "mj261_rescue1", "mj261_rescue2", "mj261_rescue3" )
head(countdataframenrm)
countdataframenrm <- log(countdataframenrm[1:8]+1,2)
countdataframenrm$meanwt <- rowMeans(countdataframenrm[,1:2])
countdataframenrm$meanmut <- rowMeans(countdataframenrm[,3:5])
countdataframenrm$meanrescue <- rowMeans(countdataframenrm[,6:8])
countdataframenrm$genes <- rownames(countdataframenrm)
head(countdataframenrm)
tail(countdataframenrm)

melted <- melt(countdataframenrm, id.vars = "genes")

g <- ggplot(melted, aes(variable,value)) + geom_boxplot() + theme(text = element_text(size = 12)) + labs(x = "Experiment", y = "log2(read+1/million)") + scale_y_continuous(limits=c(0, 0.1)) + theme_classic() + scale_x_discrete(limits = c("meanwt", "meanmut", "meanrescue"))+ stat_summary(colour = "Darkblue",fun.y=mean, geom="point", shape=23, size=4)
ggsave("nonchrIVpiRNAs_combined_2.pdf", g, height = 10, width = 10)

# take the average of the replicates, repeat the graph. 

head(countdataframenrm)
countdataframenrm2 <- data.frame(rowMeans(countdataframenrm[,1:2]), rowMeans(countdataframenrm[,3:5]), rowMeans(countdataframenrm[,6:8]), countdataframenrm$genes)
rownames(countdataframenrm2) = rownames(dds)
colnames(countdataframenrm2) = c("wt","mj261","rescue", "genes" )
head(countdataframenrm2)

countdataframenrm2 <- countdataframenrm2[order(-countdataframenrm2$wt),]
head(countdataframenrm2)

# run Kolmogorov Smirnov Test

ks.test(countdataframenrm2$wt,countdataframenrm2$mj261,alternative = c("less"), exact = FALSE)

write.csv(countdataframenrm2, "allpiRNAs.csv")

countdataframenrm3 <- countdataframenrm2[order(-countdataframenrm2$wt),]
head(countdataframenrm3)

countdataframenrm4 <- subset(countdataframenrm3, countdataframenrm3$wt > countdataframenrm3$mj261 )
write.csv(countdataframenrm4, "piRNAchnages.csv")


melted <- melt(countdataframenrm2, id.vars = "genes")

g <- ggplot(melted, aes(variable,value)) + geom_boxplot() + theme(text = element_text(size = 12)) + labs(x = "Experiment", y = "log2(read+1/million)") + theme_classic() + stat_summary(colour = "Darkblue",fun.y=mean, geom="point", shape=23, size=4)
ggsave("cpombined_allpiRNAs.pdf", g, height = 10, width = 10)


